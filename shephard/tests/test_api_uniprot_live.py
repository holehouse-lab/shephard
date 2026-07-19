"""
SHEPHARD: live tests against the real UniProt REST API.

These tests are SKIPPED by default, because they need network access and
because they assert on data that UniProt can legitimately change between
releases. Run them with:

    SHEPHARD_TEST_NETWORK=1 pytest shephard/tests/test_api_uniprot_live.py

They are worth running before a release, since they are the only check that
our understanding of the UniProt REST API still matches reality - the offline
tests all use canned responses.

Holehouse Lab - Washington University in St. Louis
"""

import json
import os
import urllib.request

import pytest

from shephard.proteome import Proteome
from shephard.apis import uniprot
from shephard.exceptions import APIException

pytestmark = pytest.mark.skipif(
    os.environ.get('SHEPHARD_TEST_NETWORK') != '1',
    reason='live UniProt tests only run when SHEPHARD_TEST_NETWORK=1')

# human p53 - heavily annotated, and stable enough to assert on
TEST_ACCESSION = 'P04637'


def _uniprot_sequence(accession):
    """
    Pulls just the sequence for an accession, so the test Proteome matches
    what UniProt holds.
    """
    url = f'{uniprot.UNIPROT_REST_URL}/{accession}.json?fields=sequence'
    raw = urllib.request.urlopen(url, timeout=30).read()
    return json.loads(raw)['sequence']['value']


@pytest.fixture(scope='module')
def live_sequence():
    return _uniprot_sequence(TEST_ACCESSION)


def _proteome(accession, sequence):
    return Proteome([{'sequence': sequence, 'name': accession,
                      'unique_ID': accession, 'attributes': {}}])


def test_live_single_protein_annotation(live_sequence):
    P = _proteome(TEST_ACCESSION, live_sequence)
    protein = P.protein(TEST_ACCESSION)

    summary = uniprot.annotate_protein_with_uniprot(protein,
                                                    experimental_structures=True,
                                                    secondary_structure_track=True,
                                                    structure_coverage_track=True)

    assert summary['domains'] > 0
    assert summary['sites'] > 0
    assert summary['tracks'] == 2
    assert summary['skipped'] == 0

    # p53 has well-known regions and modified residues
    assert len(protein.get_domains_by_type('uniprot_region')) > 0
    assert len(protein.get_sites_by_type('uniprot_modified_residue', return_list=True)) > 0

    # every annotation must carry its provenance
    for d in protein.domains:
        assert d.attribute('source') == 'uniprot'
        assert d.attribute('uniprot_accession') == TEST_ACCESSION

    # structures come with usable metadata
    structures = protein.get_domains_by_type('uniprot_experimental_structure')
    assert len(structures) > 0
    assert structures[0].attribute('pdb_id') != ''
    assert structures[0].attribute('method') != ''

    # tracks span the protein exactly
    assert len(protein.track('uniprot_secondary_structure').symbols) == len(protein)
    assert len(protein.track('uniprot_structure_coverage').values) == len(protein)

    # p53 is a well-studied protein, so some of it must be structurally covered
    assert sum(protein.track('uniprot_structure_coverage').values) > 0


def test_live_batch_annotation():
    accessions = ['P04637', 'P38398', 'P0DTC2']

    entries = []
    for accession in accessions:
        entries.append({'sequence': _uniprot_sequence(accession),
                        'name': accession, 'unique_ID': accession,
                        'attributes': {}})

    P = Proteome(entries)

    summary = uniprot.annotate_proteome_with_uniprot(P, chunk_size=2)

    assert len(summary) == 3
    for accession in accessions:
        assert summary[accession]['domains'] + summary[accession]['sites'] > 0
        assert 'error' not in summary[accession]


def test_live_missing_accession_raises():
    # correctly formatted, but not a real entry
    P = _proteome('A0A000A000', 'ACDEFGHIKL')

    with pytest.raises(APIException):
        uniprot.annotate_protein_with_uniprot(P.protein('A0A000A000'))


def test_live_sequence_mismatch_is_caught():
    # right accession, wrong sequence - annotations would be positionally
    # meaningless, so this must be refused
    P = _proteome(TEST_ACCESSION, 'ACDEFGHIKLMNPQRSTVWY')

    with pytest.raises(APIException) as excinfo:
        uniprot.annotate_protein_with_uniprot(P.protein(TEST_ACCESSION))

    assert 'does not match' in str(excinfo.value)
