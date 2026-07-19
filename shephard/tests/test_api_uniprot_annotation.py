"""
SHEPHARD: tests for the UniProt REST annotation API in shephard.apis.uniprot

Every test in this file runs offline. The single network seam
(uniprot._http_get_json) is monkeypatched with canned UniProt responses, so
the parsing, batching, annotation and error-handling logic is all exercised
without touching the network. The HTTP layer itself is tested separately by
monkeypatching urllib.request.urlopen.

Tests that hit the live UniProt API live in test_api_uniprot_live.py and are
skipped unless SHEPHARD_TEST_NETWORK=1 is set.

Holehouse Lab - Washington University in St. Louis
"""

import io
import json
import urllib.error
import urllib.request

import pytest

from shephard.proteome import Proteome
from shephard.apis import uniprot
from shephard.exceptions import APIException, InterfaceException

# 60-residue test sequence
SEQ = 'ACDEFGHIKLMNPQRSTVWY' * 3


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _proteome(accessions=('P04637',), seq=SEQ):
    return Proteome([{'sequence': seq, 'name': f'protein {a}',
                      'unique_ID': a, 'attributes': {}} for a in accessions])


def _feature(feature_type, start, end, description='', evidences=None,
             feature_id=None, cross_references=None):
    """
    Builds a UniProt-shaped feature dictionary.
    """
    feature = {'type': feature_type,
               'location': {'start': {'value': start, 'modifier': 'EXACT'},
                            'end': {'value': end, 'modifier': 'EXACT'}},
               'description': description}

    if evidences is not None:
        feature['evidences'] = evidences
    if feature_id is not None:
        feature['featureId'] = feature_id
    if cross_references is not None:
        feature['featureCrossReferences'] = cross_references

    return feature


def _record(accession='P04637', features=None, sequence=SEQ, pdb=None,
            secondary_accessions=None):
    """
    Builds a UniProt-shaped record.
    """
    record = {'primaryAccession': accession,
              'features': features if features is not None else [],
              'sequence': {'value': sequence, 'length': len(sequence)}}

    if pdb is not None:
        record['uniProtKBCrossReferences'] = pdb
    if secondary_accessions is not None:
        record['secondaryAccessions'] = secondary_accessions

    return record


def _pdb_xref(pdb_id, chains, method='X-ray', resolution='2.10 A'):
    return {'database': 'PDB', 'id': pdb_id,
            'properties': [{'key': 'Method', 'value': method},
                           {'key': 'Resolution', 'value': resolution},
                           {'key': 'Chains', 'value': chains}]}


class FakeAPI:
    """
    Stands in for uniprot._http_get_json. Records every URL requested so tests
    can assert on batching behaviour, and returns whichever of the canned
    records were asked for.
    """

    def __init__(self, records):
        self.records = {r['primaryAccession']: r for r in records}

        # index secondary accessions too, mirroring what UniProt does
        for r in records:
            for secondary in r.get('secondaryAccessions', []):
                self.records.setdefault(secondary, r)

        self.urls = []

    def __call__(self, url, timeout=30, n_retries=3, _sleep=None):
        self.urls.append(url)

        # pull the requested accessions back out of the URL
        query = url.split('accessions=')[1].split('&')[0]
        requested = [a for a in query.replace('%2C', ',').split(',') if a]

        results = []
        seen = set()
        for accession in requested:
            record = self.records.get(accession)
            if record is not None and record['primaryAccession'] not in seen:
                seen.add(record['primaryAccession'])
                results.append(record)

        return {'results': results}

    @property
    def n_calls(self):
        return len(self.urls)


@pytest.fixture
def fake_api(monkeypatch):
    """
    Installs a FakeAPI over the network seam and hands back an installer so
    each test can define the records it wants.
    """

    def _install(records):
        api = FakeAPI(records)
        monkeypatch.setattr(uniprot, '_http_get_json', api)
        return api

    return _install


# ---------------------------------------------------------------------------
# accession validation and group discovery
# ---------------------------------------------------------------------------

def test_valid_accessions():
    for accession in ['P04637', 'Q9UJX3', 'A0A023GPI8', 'P04637-2', 'O95905']:
        assert uniprot.is_valid_uniprot_accession(accession) is True


def test_invalid_accessions():
    for accession in ['', 'NOTANACC', 'p04637', 'P0463', 'MY_PROTEIN_1',
                      '12345', None, 12345, ['P04637']]:
        assert uniprot.is_valid_uniprot_accession(accession) is False


def test_annotation_groups_are_discoverable():
    groups = uniprot.uniprot_annotation_groups()

    assert groups == sorted(groups)
    assert 'domains' in groups
    assert 'modified_residues' in groups
    assert 'disulfide_bonds' in groups

    # every advertised group must be usable
    for group in groups:
        assert group in uniprot.UNIPROT_ANNOTATION_GROUPS


def test_unknown_annotation_group_raises():
    with pytest.raises(APIException):
        uniprot._resolve_annotation_groups({'not_a_real_group': True})


# ---------------------------------------------------------------------------
# field construction
# ---------------------------------------------------------------------------

def test_build_fields_requests_only_what_is_needed():
    fields = uniprot._build_fields(['domains'], False, False)

    assert fields == ['accession', 'ft_domain']

    # sequence only when we intend to verify it
    assert 'sequence' in uniprot._build_fields(['domains'], False, True)

    # PDB cross-references only when structures were requested
    assert 'xref_pdb' in uniprot._build_fields(['domains'], True, False)
    assert 'xref_pdb' not in fields


def test_build_fields_deduplicates_shared_fields():
    # transmembrane and molecule_processing both pull several fields; nothing
    # should ever appear twice
    fields = uniprot._build_fields(['transmembrane', 'molecule_processing',
                                    'domains'], True, True)
    assert len(fields) == len(set(fields))


# ---------------------------------------------------------------------------
# low level feature parsing
# ---------------------------------------------------------------------------

def test_feature_positions():
    assert uniprot._feature_positions(_feature('Domain', 5, 10)) == (5, 10)

    # string positions are tolerated
    assert uniprot._feature_positions(
        {'location': {'start': {'value': '5'}, 'end': {'value': '10'}}}) == (5, 10)


def test_feature_positions_unusable():
    # unknown position
    assert uniprot._feature_positions(
        {'location': {'start': {'value': None}, 'end': {'value': 10}}}) is None

    # no location at all
    assert uniprot._feature_positions({'type': 'Domain'}) is None

    # malformed location
    assert uniprot._feature_positions({'location': 'nonsense'}) is None
    assert uniprot._feature_positions({'location': {'start': 5, 'end': 10}}) is None

    # end before start
    assert uniprot._feature_positions(_feature('Domain', 20, 5)) is None

    # non-numeric
    assert uniprot._feature_positions(
        {'location': {'start': {'value': 'x'}, 'end': {'value': 10}}}) is None


def test_format_evidence():
    feature = _feature('Domain', 1, 10, evidences=[
        {'evidenceCode': 'ECO:0000269', 'source': 'PubMed', 'id': '25732823'},
        {'evidenceCode': 'ECO:0000255'}])

    assert uniprot._format_evidence(feature) == 'ECO:0000269|PubMed:25732823,ECO:0000255'

    # no evidence is an empty string, not a crash
    assert uniprot._format_evidence(_feature('Domain', 1, 10)) == ''
    assert uniprot._format_evidence({'evidences': ['nonsense']}) == ''


def test_format_cross_references():
    feature = _feature('Binding site', 1, 1,
                       cross_references=[{'database': 'ChEBI', 'id': 'CHEBI:29105'}])

    assert uniprot._format_cross_references(feature) == 'ChEBI:CHEBI:29105'
    assert uniprot._format_cross_references(_feature('Domain', 1, 10)) == ''


def test_annotation_name_slugging():
    assert uniprot._annotation_name('Beta strand', 'uniprot_') == 'uniprot_beta_strand'
    assert uniprot._annotation_name('Cross-link', 'uniprot_') == 'uniprot_cross_link'
    assert uniprot._annotation_name('Zinc finger', 'x_') == 'x_zinc_finger'
    assert uniprot._annotation_name('Modified residue', '') == 'modified_residue'


def test_parse_pdb_chain_ranges():
    assert uniprot._parse_pdb_chain_ranges('A/C=324-358') == [('A/C', 324, 358)]
    assert uniprot._parse_pdb_chain_ranges('A=1-100, B=200-300') == [('A', 1, 100), ('B', 200, 300)]

    # malformed entries are skipped rather than fatal
    assert uniprot._parse_pdb_chain_ranges('garbage') == []
    assert uniprot._parse_pdb_chain_ranges('A=notarange') == []
    assert uniprot._parse_pdb_chain_ranges('A=300-100') == []
    assert uniprot._parse_pdb_chain_ranges('') == []
    assert uniprot._parse_pdb_chain_ranges(None) == []

    # a good range survives alongside a bad one
    assert uniprot._parse_pdb_chain_ranges('A=1-10, junk') == [('A', 1, 10)]


# ---------------------------------------------------------------------------
# single protein annotation
# ---------------------------------------------------------------------------

def test_annotate_protein_domains_and_attributes(fake_api):
    fake_api([_record(features=[
        _feature('Domain', 5, 20, description='Protein kinase',
                 evidences=[{'evidenceCode': 'ECO:0000255', 'source': 'PROSITE-ProRule', 'id': 'PRU00159'}],
                 feature_id='PRO_0000012345')])])

    P = _proteome()
    protein = P.protein('P04637')

    summary = uniprot.annotate_protein_with_uniprot(protein)

    assert summary['domains'] == 1
    assert summary['sites'] == 0
    assert summary['skipped'] == 0

    d = protein.get_domains_by_type('uniprot_domain')[0]
    assert d.start == 5
    assert d.end == 20
    assert d.attribute('source') == 'uniprot'
    assert d.attribute('uniprot_accession') == 'P04637'
    assert d.attribute('uniprot_feature_type') == 'Domain'
    assert d.attribute('description') == 'Protein kinase'
    assert d.attribute('evidence') == 'ECO:0000255|PROSITE-ProRule:PRU00159'
    assert d.attribute('uniprot_feature_id') == 'PRO_0000012345'


def test_annotate_protein_single_residue_becomes_site(fake_api):
    fake_api([_record(features=[
        _feature('Modified residue', 12, 12, description='Phosphoserine'),
        _feature('Binding site', 30, 34, description='ATP')])])

    P = _proteome()
    protein = P.protein('P04637')

    summary = uniprot.annotate_protein_with_uniprot(protein)

    # single residue -> Site
    sites = protein.get_sites_by_type('uniprot_modified_residue', return_list=True)
    assert len(sites) == 1
    assert sites[0].position == 12
    assert sites[0].attribute('description') == 'Phosphoserine'

    # multi-residue feature of a site-ish type -> Domain
    domains = protein.get_domains_by_type('uniprot_binding_site')
    assert len(domains) == 1
    assert (domains[0].start, domains[0].end) == (30, 34)

    assert summary['sites'] == 1
    assert summary['domains'] == 1


def test_annotate_protein_region_of_one_residue_stays_a_domain(fake_api):
    # a region-mode feature must stay a Domain even when it covers one residue
    fake_api([_record(features=[_feature('Region', 7, 7, description='Tiny region')])])

    P = _proteome()
    protein = P.protein('P04637')
    uniprot.annotate_protein_with_uniprot(protein)

    assert len(protein.get_domains_by_type('uniprot_region')) == 1
    assert len(protein.sites) == 0


def test_annotate_protein_disulfide_becomes_paired_sites(fake_api):
    fake_api([_record(features=[_feature('Disulfide bond', 10, 40)])])

    P = _proteome()
    protein = P.protein('P04637')

    summary = uniprot.annotate_protein_with_uniprot(protein)

    sites = protein.get_sites_by_type('uniprot_disulfide_bond', return_list=True)
    assert summary['sites'] == 2
    assert sorted([s.position for s in sites]) == [10, 40]

    # each partner records where the other one is - a bond is not a region
    by_position = {s.position: s for s in sites}
    assert by_position[10].attribute('partner_position') == '40'
    assert by_position[40].attribute('partner_position') == '10'

    # and no Domain was created spanning the two
    assert len(protein.domains) == 0


def test_annotate_protein_interchain_crosslink_is_single_site(fake_api):
    # a cross-link to another chain is reported at one position
    fake_api([_record(features=[_feature('Cross-link', 15, 15,
                                         description='Glycyl lysine isopeptide')])])

    P = _proteome()
    protein = P.protein('P04637')

    summary = uniprot.annotate_protein_with_uniprot(protein)

    assert summary['sites'] == 1
    site = protein.get_sites_by_type('uniprot_cross_link', return_list=True)[0]
    assert site.position == 15
    assert 'partner_position' not in site.attributes


def test_annotate_protein_only_requested_groups_are_applied(fake_api):
    fake_api([_record(features=[_feature('Domain', 5, 20),
                                _feature('Natural variant', 8, 8),
                                _feature('Mutagenesis', 9, 9)])])

    P = _proteome()
    protein = P.protein('P04637')

    # natural_variants and mutagenesis are off by default
    uniprot.annotate_protein_with_uniprot(protein)

    assert len(protein.get_domains_by_type('uniprot_domain')) == 1
    assert len(protein.sites) == 0


def test_annotate_protein_opting_into_variants(fake_api):
    fake_api([_record(features=[_feature('Natural variant', 8, 8, description='R -> H')])])

    P = _proteome()
    protein = P.protein('P04637')

    uniprot.annotate_protein_with_uniprot(protein, natural_variants=True)

    sites = protein.get_sites_by_type('uniprot_natural_variant', return_list=True)
    assert len(sites) == 1
    assert sites[0].attribute('description') == 'R -> H'


def test_annotate_protein_custom_prefix(fake_api):
    fake_api([_record(features=[_feature('Domain', 5, 20)])])

    P = _proteome()
    protein = P.protein('P04637')

    uniprot.annotate_protein_with_uniprot(protein, prefix='up_')

    assert len(protein.get_domains_by_type('up_domain')) == 1


def test_annotate_protein_overlapping_features_do_not_clash(fake_api):
    # UniProt regularly reports distinct features that share coordinates -
    # autoname must keep them all
    fake_api([_record(features=[
        _feature('Region', 1, 30, description='Interaction with A'),
        _feature('Region', 1, 30, description='Interaction with B')])])

    P = _proteome()
    protein = P.protein('P04637')

    summary = uniprot.annotate_protein_with_uniprot(protein)

    assert summary['domains'] == 2
    descriptions = sorted([d.attribute('description')
                           for d in protein.get_domains_by_type('uniprot_region')])
    assert descriptions == ['Interaction with A', 'Interaction with B']


# ---------------------------------------------------------------------------
# structures and tracks
# ---------------------------------------------------------------------------

def test_annotate_protein_experimental_structures(fake_api):
    fake_api([_record(pdb=[_pdb_xref('1ABC', 'A/B=10-25', method='X-ray', resolution='1.90 A'),
                           _pdb_xref('2XYZ', 'A=40-55', method='NMR', resolution='-')])])

    P = _proteome()
    protein = P.protein('P04637')

    summary = uniprot.annotate_protein_with_uniprot(protein,
                                                    experimental_structures=True)

    structures = protein.get_domains_by_type('uniprot_experimental_structure')
    assert summary['domains'] == 2
    assert len(structures) == 2

    by_pdb = {d.attribute('pdb_id'): d for d in structures}

    assert (by_pdb['1ABC'].start, by_pdb['1ABC'].end) == (10, 25)
    assert by_pdb['1ABC'].attribute('method') == 'X-ray'
    assert by_pdb['1ABC'].attribute('resolution') == '1.90 A'
    assert by_pdb['1ABC'].attribute('chains') == 'A/B'
    assert by_pdb['1ABC'].attribute('source') == 'uniprot'

    assert by_pdb['2XYZ'].attribute('method') == 'NMR'


def test_annotate_protein_secondary_structure_track(fake_api):
    fake_api([_record(features=[_feature('Helix', 1, 5),
                                _feature('Beta strand', 10, 12),
                                _feature('Turn', 20, 20)])])

    P = _proteome()
    protein = P.protein('P04637')

    summary = uniprot.annotate_protein_with_uniprot(protein,
                                                    secondary_structure_track=True)

    symbols = protein.track('uniprot_secondary_structure').symbols

    assert summary['tracks'] == 1
    assert len(symbols) == len(SEQ)
    assert ''.join(symbols[0:5]) == 'HHHHH'
    assert ''.join(symbols[9:12]) == 'EEE'
    assert symbols[19] == 'T'
    assert symbols[5] == '-'

    # the track alone must not create secondary structure Domains
    assert len(protein.domains) == 0


def test_annotate_protein_secondary_structure_as_domains(fake_api):
    fake_api([_record(features=[_feature('Helix', 1, 5)])])

    P = _proteome()
    protein = P.protein('P04637')

    uniprot.annotate_protein_with_uniprot(protein, secondary_structure=True)

    assert len(protein.get_domains_by_type('uniprot_helix')) == 1


def test_annotate_protein_structure_coverage_track(fake_api):
    fake_api([_record(pdb=[_pdb_xref('1ABC', 'A=5-10')])])

    P = _proteome()
    protein = P.protein('P04637')

    summary = uniprot.annotate_protein_with_uniprot(protein,
                                                    structure_coverage_track=True)

    values = protein.track('uniprot_structure_coverage').values

    assert summary['tracks'] == 1
    assert len(values) == len(SEQ)
    assert values[0:4] == [0.0, 0.0, 0.0, 0.0]
    assert values[4:10] == [1.0]*6
    assert values[10] == 0.0

    # coverage track alone does not add structure Domains
    assert len(protein.domains) == 0


def test_structure_coverage_track_clamps_to_protein(fake_api):
    # a structure that runs off the end of the local sequence must still give
    # a track of exactly the right length
    fake_api([_record(pdb=[_pdb_xref('1ABC', f'A=50-{len(SEQ) + 40}')])])

    P = _proteome()
    protein = P.protein('P04637')

    uniprot.annotate_protein_with_uniprot(protein, structure_coverage_track=True)

    values = protein.track('uniprot_structure_coverage').values
    assert len(values) == len(SEQ)
    assert values[-1] == 1.0


# ---------------------------------------------------------------------------
# error handling - single protein
# ---------------------------------------------------------------------------

def test_out_of_range_feature_raises_when_safe(fake_api):
    fake_api([_record(features=[_feature('Domain', 5, len(SEQ) + 50)])])

    P = _proteome()
    protein = P.protein('P04637')

    with pytest.raises(APIException) as excinfo:
        uniprot.annotate_protein_with_uniprot(protein)

    # the message must explain WHY, since a sequence version mismatch is the
    # usual cause and is the thing the user needs to go and fix
    assert 'falls outside the protein' in str(excinfo.value)
    assert 'UniProt sequence differ' in str(excinfo.value)


def test_out_of_range_feature_skipped_when_not_safe(fake_api):
    fake_api([_record(features=[_feature('Domain', 5, len(SEQ) + 50),
                                _feature('Domain', 1, 10)])])

    P = _proteome()
    protein = P.protein('P04637')

    summary = uniprot.annotate_protein_with_uniprot(protein, safe=False, verbose=False)

    assert summary['skipped'] == 1
    assert summary['domains'] == 1
    assert len(protein.domains) == 1

    # the surviving domain must be the in-range one
    assert (protein.domains[0].start, protein.domains[0].end) == (1, 10)


def test_out_of_range_site_is_caught(fake_api):
    fake_api([_record(features=[_feature('Modified residue', len(SEQ) + 10, len(SEQ) + 10)])])

    P = _proteome()

    with pytest.raises(APIException) as excinfo:
        uniprot.annotate_protein_with_uniprot(P.protein('P04637'))

    assert 'falls outside the protein' in str(excinfo.value)
    assert len(P.protein('P04637').sites) == 0


def test_feature_without_position_is_skipped(fake_api):
    unusable = {'type': 'Domain', 'description': 'somewhere',
                'location': {'start': {'value': None}, 'end': {'value': None}}}

    fake_api([_record(features=[unusable])])

    P = _proteome()
    protein = P.protein('P04637')

    summary = uniprot.annotate_protein_with_uniprot(protein, safe=False, verbose=False)

    assert summary['skipped'] == 1
    assert summary['domains'] == 0

    # and with safe=True it is an error rather than a silent omission
    with pytest.raises(APIException):
        uniprot.annotate_protein_with_uniprot(P.protein('P04637'))


def test_sequence_mismatch_raises_when_safe(fake_api):
    fake_api([_record(sequence='MMMMMMMMMM', features=[_feature('Domain', 1, 5)])])

    P = _proteome()
    protein = P.protein('P04637')

    with pytest.raises(APIException):
        uniprot.annotate_protein_with_uniprot(protein)


def test_sequence_mismatch_skips_when_not_safe(fake_api):
    fake_api([_record(sequence='MMMMMMMMMM', features=[_feature('Domain', 1, 5)])])

    P = _proteome()
    protein = P.protein('P04637')

    summary = uniprot.annotate_protein_with_uniprot(protein, safe=False, verbose=False)

    assert summary['error'] == 'sequence mismatch'
    assert len(protein.domains) == 0


def test_sequence_verification_can_be_disabled(fake_api):
    fake_api([_record(sequence='MMMMMMMMMM', features=[_feature('Domain', 1, 5)])])

    P = _proteome()
    protein = P.protein('P04637')

    summary = uniprot.annotate_protein_with_uniprot(protein, verify_sequence=False)

    assert summary['domains'] == 1


def test_missing_record_raises_when_safe(fake_api):
    fake_api([])          # UniProt returns nothing for this accession

    P = _proteome()

    with pytest.raises(APIException):
        uniprot.annotate_protein_with_uniprot(P.protein('P04637'))


def test_missing_record_reported_when_not_safe(fake_api):
    fake_api([])

    P = _proteome()
    summary = uniprot.annotate_protein_with_uniprot(P.protein('P04637'),
                                                    safe=False, verbose=False)

    assert summary['error'] == 'no record returned'


def test_invalid_accession_handling(fake_api):
    api = fake_api([_record()])

    P = _proteome(accessions=['MY_PROTEIN_1'])
    protein = P.protein('MY_PROTEIN_1')

    with pytest.raises(APIException):
        uniprot.annotate_protein_with_uniprot(protein)

    summary = uniprot.annotate_protein_with_uniprot(protein, safe=False, verbose=False)
    assert summary['error'] == 'invalid accession'

    # a malformed accession must never reach the API
    assert api.n_calls == 0


def test_explicit_accession_overrides_unique_id(fake_api):
    api = fake_api([_record(accession='P04637', features=[_feature('Domain', 1, 10)])])

    P = _proteome(accessions=['my_local_name'])
    protein = P.protein('my_local_name')

    summary = uniprot.annotate_protein_with_uniprot(protein, accession='P04637')

    assert summary['domains'] == 1
    assert 'P04637' in api.urls[0]


def test_requesting_no_annotations_raises(fake_api):
    fake_api([_record()])
    P = _proteome()

    with pytest.raises(APIException):
        uniprot.annotate_protein_with_uniprot(
            P.protein('P04637'),
            domains=False, regions=False, motifs=False, repeats=False,
            zinc_fingers=False, coiled_coils=False, compositional_bias=False,
            dna_binding=False, transmembrane=False, binding_sites=False,
            active_sites=False, other_sites=False, modified_residues=False,
            glycosylation=False, lipidation=False, disulfide_bonds=False,
            cross_links=False)


def test_non_protein_input_raises(fake_api):
    fake_api([_record()])

    with pytest.raises(InterfaceException):
        uniprot.annotate_protein_with_uniprot('not a protein')


# ---------------------------------------------------------------------------
# proteome (batch) annotation
# ---------------------------------------------------------------------------

def test_annotate_proteome_batches_requests(fake_api):
    accessions = ['P04637', 'P38398', 'Q9UJX3', 'O95905', 'P0DTC2']
    api = fake_api([_record(accession=a, features=[_feature('Domain', 1, 10)])
                    for a in accessions])

    P = _proteome(accessions=accessions)

    summary = uniprot.annotate_proteome_with_uniprot(P, chunk_size=2)

    # 5 accessions in chunks of 2 is 3 calls, NOT 5
    assert api.n_calls == 3

    assert len(summary) == 5
    for accession in accessions:
        assert summary[accession]['domains'] == 1
        assert len(P.protein(accession).domains) == 1


def test_annotate_proteome_single_call_by_default(fake_api):
    accessions = ['P04637', 'P38398', 'Q9UJX3']
    api = fake_api([_record(accession=a) for a in accessions])

    P = _proteome(accessions=accessions)
    uniprot.annotate_proteome_with_uniprot(P)

    assert api.n_calls == 1


def test_annotate_proteome_deduplicates_accessions(fake_api):
    api = fake_api([_record(accession='P04637', features=[_feature('Domain', 1, 10)])])

    # two local proteins that both map to one UniProt accession - e.g. a
    # sequence that appears twice in a FASTA under different local names
    P = Proteome([{'sequence': SEQ, 'name': 'a', 'unique_ID': 'local_a', 'attributes': {}},
                  {'sequence': SEQ, 'name': 'b', 'unique_ID': 'local_b', 'attributes': {}}])

    summary = uniprot.annotate_proteome_with_uniprot(
        P, chunk_size=1, accession_from_protein=lambda p: 'P04637')

    # one accession, so one call even though chunk_size is 1 and there are
    # two proteins
    assert api.n_calls == 1

    # ... and both proteins were annotated from that single record
    assert summary['local_a']['domains'] == 1
    assert summary['local_b']['domains'] == 1
    assert len(P.protein('local_a').domains) == 1
    assert len(P.protein('local_b').domains) == 1


def test_annotate_proteome_resolves_secondary_accessions(fake_api):
    # the proteome uses an old accession that UniProt now lists as secondary
    fake_api([_record(accession='P04637', secondary_accessions=['Q15086'],
                      features=[_feature('Domain', 1, 10)])])

    P = _proteome(accessions=['Q15086'])

    summary = uniprot.annotate_proteome_with_uniprot(P)

    assert summary['Q15086']['domains'] == 1


def test_annotate_proteome_accession_from_protein(fake_api):
    api = fake_api([_record(accession='P04637', features=[_feature('Domain', 1, 10)])])

    P = Proteome([{'sequence': SEQ, 'name': 'n', 'unique_ID': 'local_1',
                   'attributes': {'accession': 'P04637'}}])

    summary = uniprot.annotate_proteome_with_uniprot(
        P, accession_from_protein=lambda p: p.attribute('accession'))

    assert summary['local_1']['domains'] == 1
    assert 'P04637' in api.urls[0]


def test_annotate_proteome_mixed_failures_when_not_safe(fake_api):
    fake_api([_record(accession='P04637', features=[_feature('Domain', 1, 10)]),
              _record(accession='P38398', sequence='MMMM')])

    P = _proteome(accessions=['P04637', 'P38398', 'Q9UJX3'])
    P.add_protein(SEQ, 'local', 'NOT_AN_ACCESSION')

    summary = uniprot.annotate_proteome_with_uniprot(P, safe=False, verbose=False)

    assert summary['P04637']['domains'] == 1
    assert summary['P38398']['error'] == 'sequence mismatch'
    assert summary['Q9UJX3']['error'] == 'no record returned'
    assert summary['NOT_AN_ACCESSION']['error'] == 'invalid accession'

    # the good protein was still annotated despite its neighbours failing
    assert len(P.protein('P04637').domains) == 1


def test_annotate_proteome_missing_record_raises_when_safe(fake_api):
    fake_api([_record(accession='P04637')])

    P = _proteome(accessions=['P04637', 'Q9UJX3'])

    with pytest.raises(APIException):
        uniprot.annotate_proteome_with_uniprot(P)


def test_annotate_proteome_non_proteome_input_raises(fake_api):
    fake_api([])

    with pytest.raises(InterfaceException):
        uniprot.annotate_proteome_with_uniprot('not a proteome')


def test_proteome_and_protein_paths_agree(fake_api):
    """
    The single-protein and batch functions must produce identical annotations
    from the same record - that is the whole point of sharing the engine.
    """
    features = [_feature('Domain', 5, 20, description='d'),
                _feature('Modified residue', 12, 12, description='m'),
                _feature('Disulfide bond', 3, 30),
                _feature('Region', 1, 40, description='r')]

    fake_api([_record(accession='P04637', features=features,
                      pdb=[_pdb_xref('1ABC', 'A=10-25')])])

    single = _proteome()
    uniprot.annotate_protein_with_uniprot(single.protein('P04637'),
                                          experimental_structures=True,
                                          secondary_structure_track=True)

    batch = _proteome()
    uniprot.annotate_proteome_with_uniprot(batch,
                                           experimental_structures=True,
                                           secondary_structure_track=True)

    p_single = single.protein('P04637')
    p_batch = batch.protein('P04637')

    assert sorted([(d.start, d.end, d.domain_type) for d in p_single.domains]) == \
           sorted([(d.start, d.end, d.domain_type) for d in p_batch.domains])

    assert sorted([(s.position, s.site_type) for s in p_single.sites]) == \
           sorted([(s.position, s.site_type) for s in p_batch.sites])

    assert p_single.track_names == p_batch.track_names


def test_annotate_proteome_requests_expected_fields(fake_api):
    api = fake_api([_record(accession='P04637')])

    P = _proteome()
    uniprot.annotate_proteome_with_uniprot(P, domains=True, regions=False,
                                           motifs=False, repeats=False,
                                           zinc_fingers=False, coiled_coils=False,
                                           compositional_bias=False,
                                           dna_binding=False, transmembrane=False,
                                           binding_sites=False, active_sites=False,
                                           other_sites=False, modified_residues=False,
                                           glycosylation=False, lipidation=False,
                                           disulfide_bonds=False, cross_links=False)

    url = api.urls[0]
    assert 'ft_domain' in url

    # nothing else should have been requested
    assert 'ft_variant' not in url
    assert 'xref_pdb' not in url


def test_secondary_structure_track_requests_structure_fields(fake_api):
    api = fake_api([_record(accession='P04637')])

    P = _proteome()
    uniprot.annotate_proteome_with_uniprot(P, secondary_structure_track=True)

    url = api.urls[0]
    for field in ['ft_helix', 'ft_strand', 'ft_turn']:
        assert field in url


# ---------------------------------------------------------------------------
# batching internals
# ---------------------------------------------------------------------------

def test_fetch_records_rejects_bad_chunk_size(fake_api):
    fake_api([_record()])

    with pytest.raises(APIException):
        uniprot._fetch_uniprot_records(['P04637'], ['accession'], chunk_size=0)


def test_fetch_records_rejects_non_object_response(monkeypatch):
    monkeypatch.setattr(uniprot, '_http_get_json',
                        lambda url, timeout=30, n_retries=3: ['not', 'a', 'dict'])

    with pytest.raises(APIException):
        uniprot._fetch_uniprot_records(['P04637'], ['accession'])


def test_fetch_records_tolerates_record_without_accession(monkeypatch):
    monkeypatch.setattr(uniprot, '_http_get_json',
                        lambda url, timeout=30, n_retries=3: {'results': [{'features': []}]})

    assert uniprot._fetch_uniprot_records(['P04637'], ['accession']) == {}


# ---------------------------------------------------------------------------
# HTTP layer
# ---------------------------------------------------------------------------

class FakeResponse(io.BytesIO):
    """
    Minimal stand-in for the object urlopen returns.
    """

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
        return False


def _http_error(code, body=b'', headers=None):
    return urllib.error.HTTPError('http://example.com', code, 'error',
                                  headers if headers is not None else {},
                                  io.BytesIO(body))


def test_http_get_json_success(monkeypatch):
    payload = {'results': [{'primaryAccession': 'P04637'}]}
    monkeypatch.setattr(urllib.request, 'urlopen',
                        lambda request, timeout=None: FakeResponse(json.dumps(payload).encode()))

    assert uniprot._http_get_json('http://example.com') == payload


def test_http_get_json_sets_user_agent(monkeypatch):
    captured = {}

    def _urlopen(request, timeout=None):
        captured['headers'] = request.headers
        return FakeResponse(b'{}')

    monkeypatch.setattr(urllib.request, 'urlopen', _urlopen)
    uniprot._http_get_json('http://example.com')

    # urllib capitalises header names
    assert 'SHEPHARD' in captured['headers']['User-agent']


def test_http_get_json_retries_transient_failure(monkeypatch):
    attempts = {'n': 0}

    def _urlopen(request, timeout=None):
        attempts['n'] = attempts['n'] + 1
        if attempts['n'] < 3:
            raise _http_error(503)
        return FakeResponse(b'{"results": []}')

    monkeypatch.setattr(urllib.request, 'urlopen', _urlopen)

    result = uniprot._http_get_json('http://example.com', n_retries=3, _sleep=lambda s: None)

    assert result == {'results': []}
    assert attempts['n'] == 3


def test_http_get_json_gives_up_after_retries(monkeypatch):
    attempts = {'n': 0}

    def _urlopen(request, timeout=None):
        attempts['n'] = attempts['n'] + 1
        raise _http_error(500)

    monkeypatch.setattr(urllib.request, 'urlopen', _urlopen)

    with pytest.raises(APIException):
        uniprot._http_get_json('http://example.com', n_retries=2, _sleep=lambda s: None)

    # the original attempt plus two retries
    assert attempts['n'] == 3


def test_http_get_json_respects_retry_after(monkeypatch):
    slept = []
    attempts = {'n': 0}

    def _urlopen(request, timeout=None):
        attempts['n'] = attempts['n'] + 1
        if attempts['n'] == 1:
            raise _http_error(429, headers={'Retry-After': '7'})
        return FakeResponse(b'{}')

    monkeypatch.setattr(urllib.request, 'urlopen', _urlopen)
    uniprot._http_get_json('http://example.com', _sleep=slept.append)

    assert slept == [7.0]


def test_http_get_json_backs_off_exponentially(monkeypatch):
    slept = []
    attempts = {'n': 0}

    def _urlopen(request, timeout=None):
        attempts['n'] = attempts['n'] + 1
        if attempts['n'] < 4:
            raise _http_error(502)
        return FakeResponse(b'{}')

    monkeypatch.setattr(urllib.request, 'urlopen', _urlopen)
    uniprot._http_get_json('http://example.com', n_retries=5, _sleep=slept.append)

    assert slept == [1.0, 2.0, 4.0]


def test_http_get_json_does_not_retry_client_error(monkeypatch):
    attempts = {'n': 0}
    body = json.dumps({'messages': ["Accession 'NOTANACC' has invalid format."]}).encode()

    def _urlopen(request, timeout=None):
        attempts['n'] = attempts['n'] + 1
        raise _http_error(400, body=body)

    monkeypatch.setattr(urllib.request, 'urlopen', _urlopen)

    with pytest.raises(APIException) as excinfo:
        uniprot._http_get_json('http://example.com', n_retries=3, _sleep=lambda s: None)

    # a 400 is our fault, not a blip - so no retrying, and UniProt's own
    # message must be surfaced
    assert attempts['n'] == 1
    assert 'invalid format' in str(excinfo.value)
    assert '400' in str(excinfo.value)


def test_http_get_json_handles_network_failure(monkeypatch):
    def _urlopen(request, timeout=None):
        raise urllib.error.URLError('name resolution failed')

    monkeypatch.setattr(urllib.request, 'urlopen', _urlopen)

    with pytest.raises(APIException) as excinfo:
        uniprot._http_get_json('http://example.com', n_retries=1, _sleep=lambda s: None)

    assert 'online' in str(excinfo.value)


def test_http_get_json_handles_bad_json(monkeypatch):
    monkeypatch.setattr(urllib.request, 'urlopen',
                        lambda request, timeout=None: FakeResponse(b'<html>not json</html>'))

    with pytest.raises(APIException) as excinfo:
        uniprot._http_get_json('http://example.com')

    assert 'JSON' in str(excinfo.value)


def test_http_get_json_rejects_negative_retries():
    with pytest.raises(APIException):
        uniprot._http_get_json('http://example.com', n_retries=-1)


def test_extract_error_message_falls_back_to_raw_body():
    error = _http_error(500, body=b'internal explosion')
    assert 'internal explosion' in uniprot._extract_error_message(error)
