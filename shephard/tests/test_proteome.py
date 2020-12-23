"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

import pytest
import shephard
from shephard import proteome
from shephard.apis import uniprot, fasta
from shephard.exceptions import ProteomeException

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

def test_add_protein():

    # creating proteome and adding protein
    test_data_dir = shephard.get_data('test_data')
    proteome = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'))
    assert len(proteome.protein('O00401')) == 505
    assert len(proteome.protein('O00470')) == 390
    assert len(proteome.protein('O00472')) == 640
    assert len(proteome.protein('O00499')) == 593
    assert len(proteome.protein('O00629')) == 521
    assert len(proteome.protein('O00712')) == 420
    assert len(proteome.protein('O00716')) == 465
    assert len(proteome.protein('O14786')) == 923
    assert len(proteome.protein('Q9UJX3')) == 599

    # creating a proteome from a FASTA file (using defaul unique key)
    proteome = fasta.fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'))
    assert len(proteome.protein('1')) == 390
    assert len(proteome.protein('2')) == 640
    assert len(proteome.protein('3')) == 593
    assert len(proteome.protein('4')) == 521
    assert len(proteome.protein('5')) == 420
    assert len(proteome.protein('6')) == 465
    assert len(proteome.protein('7')) == 923
    assert len(proteome.protein('8')) == 599

    # create a proteome where FASTA header is used as uniqueID 
    proteome = fasta.fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'), use_header_as_unique_ID=True)
    assert len(proteome.protein('sp|O00470|MEIS1_HUMAN Homeobox protein Meis1 OS=Homo sapiens OX=9606 GN=MEIS1 PE=1 SV=1')) == 390
    assert len(proteome.protein('sp|O00472|ELL2_HUMAN RNA polymerase II elongation factor ELL2 OS=Homo sapiens OX=9606 GN=ELL2 PE=1 SV=2')) == 640
    assert len(proteome.protein('sp|O00499|BIN1_HUMAN Myc box-dependent-interacting protein 1 OS=Homo sapiens OX=9606 GN=BIN1 PE=1 SV=1')) == 593
    assert len(proteome.protein('sp|O00629|IMA3_HUMAN Importin subunit alpha-3 OS=Homo sapiens OX=9606 GN=KPNA4 PE=1 SV=1')) == 521
    assert len(proteome.protein('sp|O00712|NFIB_HUMAN Nuclear factor 1 B-type OS=Homo sapiens OX=9606 GN=NFIB PE=1 SV=2')) == 420
    assert len(proteome.protein('sp|O00716|E2F3_HUMAN Transcription factor E2F3 OS=Homo sapiens OX=9606 GN=E2F3 PE=1 SV=1')) == 465
    assert len(proteome.protein('sp|O14786|NRP1_HUMAN Neuropilin-1 OS=Homo sapiens OX=9606 GN=NRP1 PE=1 SV=3')) == 923
    assert len(proteome.protein('sp|Q9UJX3|APC7_HUMAN Anaphase-promoting complex subunit 7 OS=Homo sapiens OX=9606 GN=ANAPC7 PE=1 SV=4')) == 599

    #for protein in proteome:
    #    print("assert len(proteome.protein('%s')) == %i" % (protein.unique_ID, len(protein)))
    
    # check manually adding proteomes    
    local_seq = 'PPPPP'
    proteome.add_protein(local_seq, '5pp', 'U5P')
    assert proteome.protein('U5P').sequence == local_seq
    assert proteome.protein('U5P').name == '5pp'

    # should trigger exception
    with pytest.raises(ProteomeException):
        proteome.add_protein(local_seq, '5pp', 'U5P')
    proteome.add_protein('ASDF', '5pp', 'U5P', force_overwrite=True)
    assert proteome.protein('U5P').sequence == 'ASDF'
    
    

    
            


#test_add_protein()



