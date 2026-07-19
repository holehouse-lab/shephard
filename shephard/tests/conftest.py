"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

import shephard
from shephard.interfaces import si_sites, si_domains, si_tracks, si_protein_attributes
from shephard.apis import uniprot  

import os

import pytest
import sys

# Several legacy tests write outputs to a relative "output_test/" directory,
# which only resolves when pytest is invoked from inside the tests directory.
# This session-scoped, autouse fixture pins the working directory to the
# tests directory (and guarantees output_test/ exists) so the suite passes
# regardless of where pytest is run from. Newer tests use tmp_path and
# absolute paths and are unaffected.
_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture(scope='session', autouse=True)
def _stable_working_directory():
    original = os.getcwd()
    os.chdir(_TESTS_DIR)
    os.makedirs(os.path.join(_TESTS_DIR, 'output_test'), exist_ok=True)
    try:
        yield
    finally:
        os.chdir(original)

TS1_FILE = ['testset_1.fasta', 
            'TS1_domains_idr.tsv', 
            'TS1_domains_pscore.tsv', 
            'TS1_sites.tsv', 
            'TS1_tracks_pscore.tsv', 
            'TS1_protein_attributes.tsv',
            'testset_1_ptms.tsv']

test_data_dir = shephard.get_data('test_data')


def build_proteome(fn):
    return uniprot.uniprot_fasta_to_proteome(f'{test_data_dir}/{fn}')


@pytest.fixture
def TS1(request):    
    TS1_proteome = build_proteome(TS1_FILE[0])
    return TS1_proteome


@pytest.fixture
def TS1_domains(request):    
    TS1_proteome = build_proteome(TS1_FILE[0])
    si_domains.add_domains_from_file(TS1_proteome, f'{test_data_dir}/{TS1_FILE[1]}')
    
    return TS1_proteome


@pytest.fixture
def TS1_domains2(request):    
    TS1_proteome = build_proteome(TS1_FILE[0])
    si_domains.add_domains_from_file(TS1_proteome, f'{test_data_dir}/{TS1_FILE[1]}')
    si_domains.add_domains_from_file(TS1_proteome, f'{test_data_dir}/{TS1_FILE[2]}')
    
    return TS1_proteome

@pytest.fixture
def TS1_domains2_sites(request):    
    TS1_proteome = build_proteome(TS1_FILE[0])
    si_domains.add_domains_from_file(TS1_proteome, f'{test_data_dir}/{TS1_FILE[1]}')
    si_domains.add_domains_from_file(TS1_proteome, f'{test_data_dir}/{TS1_FILE[2]}')
    si_sites.add_sites_from_file(TS1_proteome, f'{test_data_dir}/{TS1_FILE[3]}')
    
    return TS1_proteome


@pytest.fixture
def TS1_domains2_sites_tracks(request):    
    TS1_proteome = build_proteome(TS1_FILE[0])

    si_domains.add_domains_from_file(TS1_proteome, f'{test_data_dir}/{TS1_FILE[1]}')
    si_domains.add_domains_from_file(TS1_proteome, f'{test_data_dir}/{TS1_FILE[2]}')
    si_sites.add_sites_from_file(TS1_proteome, f'{test_data_dir}/{TS1_FILE[3]}')
    si_tracks.add_tracks_from_file(TS1_proteome, f'{test_data_dir}/{TS1_FILE[4]}', 'values')
    
    return TS1_proteome


@pytest.fixture
def TS1_domains2_sites_tracks_protein_attributes(request):    
    TS1_proteome = build_proteome(TS1_FILE[0])

    si_domains.add_domains_from_file(TS1_proteome, f'{test_data_dir}/{TS1_FILE[1]}')
    si_domains.add_domains_from_file(TS1_proteome, f'{test_data_dir}/{TS1_FILE[2]}')
    si_sites.add_sites_from_file(TS1_proteome, f'{test_data_dir}/{TS1_FILE[3]}')
    si_sites.add_sites_from_file(TS1_proteome, f'{test_data_dir}/{TS1_FILE[6]}')
    si_tracks.add_tracks_from_file(TS1_proteome, f'{test_data_dir}/{TS1_FILE[4]}', 'values')
    si_protein_attributes.add_protein_attributes_from_file(TS1_proteome, f'{test_data_dir}/{TS1_FILE[5]}')
    
    return TS1_proteome

