"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

import pytest

names_list = ['O00401', 'O00470', 'O00472', 'O00499', 'O00629', 'O00712', 'O00716', 'O14786', 'Q9UJX3']

def test_proteome_len(TS1):    
    assert len(TS1) == 9
    assert len(TS1) == len(TS1.proteins)

    
def test_proteins_retrieve(TS1):
    # checks that all the unique IDs are there as expected
    for protein in TS1.proteins:
        assert protein in names_list


def test_record_and_unique_ID_match(TS1):
    # checks when we iterate we get a protein
    for p in TS1.proteins:
        assert p == TS1.protein(p).unique_ID



def test_iterator(TS1):
    for p in TS1:
        assert str(type(p)) == "<class 'shephard.protein.Protein'>"



def test_unique_domain_types(TS1, TS1_domains, TS1_domains2, TS1_domains2_sites, TS1_domains2_sites_tracks):
    # check we corretly count the number of unique domains

    assert len(TS1.unique_domain_types) == 0
    assert len(TS1_domains.unique_domain_types) == 1
    assert len(TS1_domains2.unique_domain_types) == 2
    assert len(TS1_domains2_sites.unique_domain_types) == 2
    assert len(TS1_domains2_sites_tracks.unique_domain_types) == 2



def test_unique_site_types(TS1, TS1_domains, TS1_domains2, TS1_domains2_sites, TS1_domains2_sites_tracks):
    # check we corretly count the number of unique sites

    assert len(TS1.unique_site_types) == 0
    assert len(TS1_domains.unique_site_types) == 0
    assert len(TS1_domains2.unique_site_types) == 0
    assert len(TS1_domains2_sites.unique_site_types) == 15
    assert len(TS1_domains2_sites_tracks.unique_site_types) == 15



def test_domain_retriev(TS1_domains):
    # check we corretly can get all domains
    assert len(TS1_domains.domains) == 21


def test_domain_retriev(TS1_domains2_sites):
    # check we corretly can get all sites
    assert len(TS1_domains2_sites.sites) == 209

def test_protein(TS1_domains2_sites):
    # check we corretly access a valid protein

    # get sequence
    len(TS1_domains2_sites.protein('O00629')) == 521
    len(TS1_domains2_sites.protein('O00629').domains) == 2
    len(TS1_domains2_sites.protein('O00629').sites) == 16

            
            
