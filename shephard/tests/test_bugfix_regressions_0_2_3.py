"""
SHEPHARD regression suite : explicit tests for every bug fixed for the 0.2.3
release.

Each test is named after the bug it guards and fails if the underlying bug is
ever reintroduced. Note that the bugs fixed on the original `streamline` branch
are covered in test_bugfix_regressions.py.

Bugs covered
------------
 1.  exceptions.py  : every SHEPHARD exception subclassed bare Exception, so
                      `except ShephardException` caught almost nothing
 2.  domain.py      : Domain.get_sites_by_type() ignored return_list
 3.  domain_tools.py: the extend_ends block used a float index and compared a
                      list slice against an int, so it crashed or silently
                      did nothing
 4.  experimental.py: enrichment was inverted for a depleted domain, and sites
                      were counted as positions-with-sites rather than sites
 5.  si_tracks.py   : track lines were written with a trailing delimiter, which
                      broke round-tripping with any non-whitespace delimiter
 6.  si_tracks.py   : an invalid mode was raised inside the per-line try, so
                      skip_bad swallowed it
 7.  si_tracks.py   : write_tracks_from_list([]) raised IndexError
 8.  si_protein_attributes.py : write_protein_attributes_from_dictionary could
                      not consume the dictionary the parser produces
 9.  si_*.py        : skip_bad did not cover malformed key:value attributes
 10. si_*.py        : attribute keys were never sanitized on write, and
                      full_clean_string hardcoded '\\t' regardless of delimiter
 11. si_sites.py    : add_sites_from_dictionary never validated the proteome
 12. fasta.py       : auto-generated unique_IDs were ints, not strs
 13. fasta.py       : non-string attributes raised AttributeError on write
 14. fasta.py       : protein names gained a trailing space on round-trip
 15. metapredict_api.py : disorder_threshold was accepted, documented, and
                      never passed to metapredict
 16. albatross_api.py / metapredict_api.py : a missing optional dependency
                      only printed at import, then failed later with NameError
 17. attribute_tools.py : skip_failed did not cover TypeError

Holehouse Lab - Washington University in St. Louis
"""

import os

import pytest

import shephard
from shephard.proteome import Proteome
from shephard.exceptions import (ShephardException, SiteException,
                                 TrackException, DomainException,
                                 ProteinException, ProteomeException,
                                 UtilitiesException, InterfaceException,
                                 APIException)
from shephard.interfaces import (interface_tools, si_tracks, si_sites,
                                 si_domains, si_protein_attributes)
from shephard.tools import domain_tools, track_tools, experimental, attribute_tools
from shephard.apis import fasta, metapredict_api, albatross_api

DATA = shephard.get_data('test_data')
SEQ = 'ACDEFGHIKLMNPQRSTVWY'


def _proteome(seq=SEQ, uid='P1'):
    return Proteome([{'sequence': seq, 'name': 'demo',
                      'unique_ID': uid, 'attributes': {}}])


# --- bug 1 -----------------------------------------------------------------
def test_bug01_exception_hierarchy():
    # every SHEPHARD exception must be catchable as a ShephardException
    for exception_class in [SiteException, TrackException, DomainException,
                            ProteinException, ProteomeException,
                            UtilitiesException, InterfaceException,
                            APIException]:
        assert issubclass(exception_class, ShephardException)

    # and this must hold in practice, not just in theory
    p = _proteome()
    with pytest.raises(ShephardException):
        p.protein('P1').add_site(500, 'phos')


# --- bug 2 -----------------------------------------------------------------
def test_bug02_domain_get_sites_by_type_return_list():
    p = _proteome().protein('P1')
    p.add_domain(2, 10, 'test_domain')
    p.add_site(4, 'phos')
    p.add_site(4, 'phos')          # two sites at ONE position
    p.add_site(6, 'phos')

    d = p.domain('test_domain_2_10')

    as_dict = d.get_sites_by_type('phos')
    as_list = d.get_sites_by_type('phos', return_list=True)

    assert isinstance(as_dict, dict)
    assert isinstance(as_list, list)

    # two occupied positions, but three sites
    assert len(as_dict) == 2
    assert len(as_list) == 3


# --- bug 3 -----------------------------------------------------------------
def test_bug03_build_domains_extend_ends():
    # a track that is disordered everywhere except the first and last residues
    p = _proteome(seq='A'*100)
    values = [0.0]*5 + [1.0]*90 + [0.0]*5
    p.protein('P1').add_track('t', values=values)

    def binerize(v):
        return track_tools.binerize(v, 0.5)

    # without extend_ends we should not reach the termini...
    d = domain_tools.build_domains_from_track_values(p, 't', binerize, 'idr',
                                                     minimum_region_size=10,
                                                     verbose=False)
    assert d['P1'][0]['start'] == 6
    assert d['P1'][0]['end'] == 95

    # ... and with extend_ends set the domain runs to both termini. Note this
    # used to raise a TypeError (float index)
    d = domain_tools.build_domains_from_track_values(p, 't', binerize, 'idr',
                                                     minimum_region_size=10,
                                                     extend_ends=5,
                                                     verbose=False)
    assert d['P1'][0]['start'] == 1
    assert d['P1'][0]['end'] == 100

    # an extend_ends value larger than the sequence must not crash either
    domain_tools.build_domains_from_track_values(p, 't', binerize, 'idr',
                                                 minimum_region_size=10,
                                                 extend_ends=500,
                                                 verbose=False)


# --- bug 4 -----------------------------------------------------------------
def test_bug04_experimental_enrichment_direction():
    # sites everywhere EXCEPT in the domain, so the domain is depleted
    p = _proteome(seq='A'*100)
    protein = p.protein('P1')
    protein.add_domain(1, 20, 'test_domain')

    for i in range(30, 100):
        protein.add_site(i, 'phos')

    d = protein.domain('test_domain_1_20')

    # a domain with no sites but sites elsewhere is maximally DEPLETED - this
    # used to return min_enrichment for the opposite (enriched) case and
    # max_enrichment here
    assert experimental.get_site_density_in_domain_normalized_by_protein(
        d, 'phos', sample_size=50) == 0.01

    # the percentile of a site-free domain must be at the bottom of the
    # distribution, and must always return a value (it used to fall off the
    # end of the function and return None when no exact float match was found)
    pct = experimental.get_site_density_percentile_normalized_by_protein(d, 'phos')
    assert pct is not None
    assert 0 < pct < 1


def test_bug04_experimental_enriched_domain():
    # sites ONLY in the domain, so the domain is enriched
    p = _proteome(seq='A'*100)
    protein = p.protein('P1')
    protein.add_domain(1, 20, 'test_domain')

    for i in range(1, 21):
        protein.add_site(i, 'phos')

    d = protein.domain('test_domain_1_20')

    enrichment = experimental.get_site_density_in_domain_normalized_by_protein(
        d, 'phos', sample_size=50)

    # every residue in the domain has a site and most of the protein has none,
    # so this must come out enriched
    assert enrichment > 1

    # a fully site-covered domain must sit at the top of the distribution
    assert experimental.get_site_density_percentile_normalized_by_protein(
        d, 'phos') == 1


# --- bug 5 -----------------------------------------------------------------
def test_bug05_track_write_no_trailing_delimiter(tmp_path):
    p = _proteome()
    p.protein('P1').add_track('vals', values=[float(i) for i in range(20)])

    outfile = str(tmp_path / 'tracks.csv')
    si_tracks.write_tracks(p, outfile, 'vals', delimiter=',')

    line = open(outfile).readlines()[0].strip()
    assert not line.endswith(',')

    # and the round-trip must work with a non-whitespace delimiter, which it
    # could not when an empty trailing field was written
    p2 = _proteome()
    si_tracks.add_tracks_from_file(p2, outfile, mode='values', delimiter=',')
    assert p2.protein('P1').track('vals').values == [float(i) for i in range(20)]


def test_bug05_symbols_track_write_no_trailing_delimiter(tmp_path):
    p = _proteome()
    p.protein('P1').add_track('syms', symbols=list('-'*20))

    outfile = str(tmp_path / 'tracks_symbols.csv')
    si_tracks.write_tracks(p, outfile, 'syms', delimiter=',')

    p2 = _proteome()
    si_tracks.add_tracks_from_file(p2, outfile, mode='symbols', delimiter=',')
    assert p2.protein('P1').track('syms').symbols == list('-'*20)


# --- bug 6 -----------------------------------------------------------------
def test_bug06_bad_track_mode_always_raises():
    track_file = f'{DATA}/TS1_tracks_pscore.tsv'

    # skip_bad must NOT swallow an invalid mode
    with pytest.raises(InterfaceException):
        si_tracks._TracksInterface(track_file, mode='not_a_mode', skip_bad=True)

    p = _proteome()
    with pytest.raises(Exception):
        si_tracks.add_tracks_from_file(p, track_file, mode='not_a_mode')


# --- bug 7 -----------------------------------------------------------------
def test_bug07_write_tracks_from_empty_list(tmp_path):
    outfile = str(tmp_path / 'empty_tracks.tsv')

    # used to raise IndexError on track_list[0]
    si_tracks.write_tracks_from_list([], outfile)

    assert os.path.exists(outfile)
    assert os.path.getsize(outfile) == 0


# --- bug 8 -----------------------------------------------------------------
def test_bug08_write_protein_attributes_from_parsed_dictionary(tmp_path):
    p = _proteome()
    p.protein('P1').add_attribute('gene', 'ABC1')

    infile = str(tmp_path / 'attributes.tsv')
    si_protein_attributes.write_protein_attributes(p, infile)

    # the parser hands back {UID : [dict, ...]}, and the writer must be able to
    # consume exactly that (it used to raise AttributeError on the list)
    parsed = si_protein_attributes.add_protein_attributes_from_file(
        p, infile, return_dictionary=True)

    outfile = str(tmp_path / 'attributes_out.tsv')
    si_protein_attributes.write_protein_attributes_from_dictionary(parsed, outfile)

    p2 = _proteome()
    si_protein_attributes.add_protein_attributes_from_file(p2, outfile)
    assert p2.protein('P1').attribute('gene') == 'ABC1'

    # a plain {UID : dict} must keep working too
    outfile2 = str(tmp_path / 'attributes_out2.tsv')
    si_protein_attributes.write_protein_attributes_from_dictionary(
        {'P1': {'gene': 'ABC1'}}, outfile2)

    p3 = _proteome()
    si_protein_attributes.add_protein_attributes_from_file(p3, outfile2)
    assert p3.protein('P1').attribute('gene') == 'ABC1'


# --- bug 9 -----------------------------------------------------------------
def test_bug09_skip_bad_covers_malformed_attributes(tmp_path):
    # the final field has no ':' so it cannot be parsed as a key-value pair
    good = 'P1\t1\t10\tidr\tkey:value\n'
    bad = 'P1\t1\t10\tidr\tthis_is_not_a_key_value_pair\n'

    infile = str(tmp_path / 'domains.tsv')
    with open(infile, 'w') as fh:
        fh.write(good)
        fh.write(bad)

    # with skip_bad the malformed line is skipped rather than raising
    parsed = si_domains._DomainsInterface(infile, skip_bad=True)
    assert len(parsed.data['P1']) == 1

    # ... and without skip_bad it still raises
    with pytest.raises(InterfaceException):
        si_domains._DomainsInterface(infile, skip_bad=False)


# --- bug 10 ----------------------------------------------------------------
def test_bug10_attribute_keys_are_cleaned_on_write(tmp_path):
    p = _proteome()

    # a key containing the delimiter and a colon used to produce a file that
    # could not be read back in
    p.protein('P1').add_attribute('bad\tkey:name', 'value')

    outfile = str(tmp_path / 'attributes.tsv')
    si_protein_attributes.write_protein_attributes(p, outfile)

    p2 = _proteome()
    si_protein_attributes.add_protein_attributes_from_file(p2, outfile)

    assert len(p2.protein('P1').attributes) == 1
    assert p2.protein('P1').attribute('bad key-name') == 'value'


def test_bug10_full_clean_string_respects_delimiter():
    # the delimiter passed must be the one that gets cleaned
    assert interface_tools.full_clean_string('a,b', ',') == 'a b'

    # and the default remains the tab character
    assert interface_tools.full_clean_string('a\tb') == 'a b'

    # colons are always cleaned, since they delimit key:value pairs
    assert interface_tools.full_clean_string('a:b', ',') == 'a-b'


# --- bug 11 ----------------------------------------------------------------
def test_bug11_add_sites_from_dictionary_validates_proteome():
    # passing something that is not a Proteome must raise, as it does for
    # every other *_from_dictionary function
    with pytest.raises(InterfaceException):
        si_sites.add_sites_from_dictionary('not_a_proteome', {})


# --- bug 12 ----------------------------------------------------------------
def test_bug12_fasta_autogenerated_unique_ids_are_strings():
    P = fasta.fasta_to_proteome(f'{DATA}/testset_1.fasta')

    for uid in P.proteins:
        assert isinstance(uid, str)

    # and the proteins must be retrievable using those IDs
    for uid in P.proteins:
        assert P.protein(uid) is not None


# --- bug 13 / 14 -----------------------------------------------------------
def test_bug13_fasta_roundtrip_with_nonstring_attributes(tmp_path):
    p = _proteome()
    p.protein('P1').add_attribute('copy_number', 12345)      # an int, not a str
    p.protein('P1').add_attribute('missing', None)

    outfile = str(tmp_path / 'out.fasta')

    # used to raise AttributeError because int has no .replace()
    fasta.proteome_to_fasta(outfile, p, include_attributes_in_header=True)

    p2 = fasta.shephard_fasta_to_proteome(outfile)
    assert p2.protein('P1').attribute('copy_number') == '12345'


def test_bug14_fasta_roundtrip_preserves_name(tmp_path):
    p = _proteome()
    p.protein('P1').add_attribute('gene', 'ABC1')

    outfile = str(tmp_path / 'out.fasta')
    fasta.proteome_to_fasta(outfile, p, include_attributes_in_header=True)

    p2 = fasta.shephard_fasta_to_proteome(outfile)

    # the name used to pick up a trailing space on every round-trip
    assert p2.protein('P1').name == 'demo'


# --- bug 15 ----------------------------------------------------------------
def test_bug15_disorder_threshold_is_passed_to_metapredict(monkeypatch):
    captured = {}

    class _FakeDisorder:
        disordered_domain_boundaries = []
        folded_domain_boundaries = []
        disorder = [0.0]*len(SEQ)

    class _FakeMeta:
        @staticmethod
        def predict_disorder(uid2seq, **kwargs):
            captured.update(kwargs)
            return {k: _FakeDisorder() for k in uid2seq}

    monkeypatch.setattr(metapredict_api, 'meta', _FakeMeta)
    monkeypatch.setattr(metapredict_api, '_check_metapredict', lambda: None)

    p = _proteome()
    metapredict_api.annotate_proteome_with_disordered_domains(
        p, disorder_threshold=0.75, show_progress_bar=False)

    assert captured['disorder_threshold'] == 0.75

    captured.clear()
    metapredict_api.annotate_proteome_with_disorder_tracks_and_disordered_domains(
        p, disorder_threshold=0.3, show_progress_bar=False)

    assert captured['disorder_threshold'] == 0.3


# --- bug 16 ----------------------------------------------------------------
def test_bug16_missing_optional_dependency_raises_api_exception(monkeypatch):
    p = _proteome()

    # a missing metapredict must give an informative APIException rather than a
    # NameError raised somewhere deep in the call
    monkeypatch.setattr(metapredict_api, 'meta', None)
    with pytest.raises(APIException):
        metapredict_api.annotate_proteome_with_disorder_track(p)

    monkeypatch.setattr(albatross_api, 'batch_predict', None)
    with pytest.raises(APIException):
        albatross_api.annotate_proteome_with_dimensions(p)


# --- bug 17 ----------------------------------------------------------------
def test_bug17_cast_attributes_skips_type_errors():
    p = _proteome()
    protein = p.protein('P1')
    protein.add_attribute('good', '1.5')
    protein.add_attribute('bad', None)          # float(None) is a TypeError

    # skip_failed must cover TypeError as well as ValueError
    attribute_tools.cast_attributes(protein, cast_type=float, skip_failed=True)

    assert protein.attribute('good') == 1.5
    assert protein.attribute('bad') is None

    # and without skip_failed the error still propagates
    protein.add_attribute('good', '1.5', safe=False)
    with pytest.raises(TypeError):
        attribute_tools.cast_attributes(protein, cast_type=float, skip_failed=False)
