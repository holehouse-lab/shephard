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

import pytest
import sys

TS1_FILE = ['testset_1.fasta', 
            'TS1_domains_idr.tsv', 
            'TS1_domains_pscore.tsv', 
            'TS1_sites.tsv', 
            'TS1_tracks_pscore.tsv', 
            'TS1_protein_attributes.tsv',
            'testset_1_ptms.tsv']

test_data_dir = shephard.get_data('test_data')


def build_proteome(fn):
    return uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir,fn))


@pytest.fixture
def TS1(request):    
    TS1_proteome = build_proteome(TS1_FILE[0])
    return TS1_proteome


@pytest.fixture
def TS1_domains(request):    
    TS1_proteome = build_proteome(TS1_FILE[0])
    si_domains.add_domains_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[1]))
    
    return TS1_proteome


@pytest.fixture
def TS1_domains2(request):    
    TS1_proteome = build_proteome(TS1_FILE[0])
    si_domains.add_domains_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[1]))
    si_domains.add_domains_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[2]))
    
    return TS1_proteome

@pytest.fixture(scope='session', autouse=True)
def TS1_domains2_sites(request):    
    TS1_proteome = build_proteome(TS1_FILE[0])
    si_domains.add_domains_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[1]))
    si_domains.add_domains_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[2]))
    si_sites.add_sites_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[3]))
    
    return TS1_proteome


@pytest.fixture
def TS1_domains2_sites_tracks(request):    
    TS1_proteome = build_proteome(TS1_FILE[0])

    si_domains.add_domains_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[1]))
    si_domains.add_domains_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[2]))
    si_sites.add_sites_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[3]))
    si_tracks.add_tracks_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[4]), 'values')
    
    return TS1_proteome


@pytest.fixture
def TS1_domains2_sites_tracks_protein_attributes(request):    
    TS1_proteome = build_proteome(TS1_FILE[0])

    si_domains.add_domains_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[1]))
    si_domains.add_domains_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[2]))
    si_sites.add_sites_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[3]))
    si_sites.add_sites_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[6]))
    si_tracks.add_tracks_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[4]), 'values')
    si_protein_attributes.add_protein_attributes_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[5]))
    
    return TS1_proteome


#si_sites.add_sites(proteome, 'proteomescout_ptm_sites.csv')
"""
@pytest.fixture(scope='session', autouse=True)
def GS6_CP(request):
    GS6_CO = cttrajectory.CTTrajectory("%s/%s"%(test_data_dir, GS6_FILES[1]), "%s/%s"%(test_data_dir, GS6_FILES[0])) 
    GS6_CP = GS6_CO.proteinTrajectoryList[0]    
    return GS6_CP


@pytest.fixture(scope='session', autouse=True)
def NTL9_CO(request):    
    NTL9_CO = cttrajectory.CTTrajectory("%s/%s"%(test_data_dir, NTL9_FILES[1]), "%s/%s"%(test_data_dir, NTL9_FILES[0])) 
    return NTL9_CO

@pytest.fixture(scope='session', autouse=True)
def NTL9_CP(request):    
    NTL9_CO = cttrajectory.CTTrajectory("%s/%s"%(test_data_dir, NTL9_FILES[1]), "%s/%s"%(test_data_dir, NTL9_FILES[0])) 
    NTL9_CP = NTL9_CO.proteinTrajectoryList[0]    
    return NTL9_CP
"""
