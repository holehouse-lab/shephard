"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

import pytest
from shephard.exceptions import ProteinException, DomainException
import numpy as np
import shephard
from shephard.apis import uniprot
from shephard import interfaces
from shephard.tools import domain_tools

names_list = ['O00401', 'O00470', 'O00472', 'O00499', 'O00629', 'O00712', 'O00716', 'O14786', 'Q9UJX3']

def test_domain_overlap():    
    """
    Note if any of these succeed the test suite will fail because it incrases
    """

    # read in a proteome
    test_data_dir = shephard.get_data('test_data')
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    P = uniprot.uniprot_fasta_to_proteome(fasta_file)

    # add a synthetic protein and extract it
    P.add_protein('A'*100,'test','test')    
    prot = P.protein('test')

    
    prot.add_domain(1,  20, 'test_domain_0') # 0
    prot.add_domain(5,  15, 'test_domain_1') # 1

    d0 = prot.get_domains_by_type('test_domain_0')[0]
    d1 = prot.get_domains_by_type('test_domain_1')[0]
    assert domain_tools.domain_overlap(d0, d1) is True
    assert domain_tools.domain_overlap(d1, d0) is True


    prot.add_domain(15, 25, 'test_domain_2') # 2
    d2 = prot.get_domains_by_type('test_domain_2')[0]
    assert domain_tools.domain_overlap(d0, d2) is True
    assert domain_tools.domain_overlap(d1, d2) is True
    assert domain_tools.domain_overlap(d2, d0) is True
    assert domain_tools.domain_overlap(d2, d1) is True

    prot.add_domain(20, 30, 'test_domain_3') # 3
    d3 = prot.get_domains_by_type('test_domain_3')[0]
    assert domain_tools.domain_overlap(d0, d3) is True
    assert domain_tools.domain_overlap(d2, d3) is True
    assert domain_tools.domain_overlap(d1, d3) is False

    assert domain_tools.domain_overlap(d3, d0) is True
    assert domain_tools.domain_overlap(d3, d2) is True
    assert domain_tools.domain_overlap(d3, d1) is False

    prot.add_domain(1, 4, 'test_domain_4') # 4
    d4 = prot.get_domains_by_type('test_domain_4')[0]
    assert domain_tools.domain_overlap(d0, d4) is True
    assert domain_tools.domain_overlap(d1, d4) is False

    assert domain_tools.domain_overlap(d4, d0) is True
    assert domain_tools.domain_overlap(d4, d1) is False

    print(prot.domains)
    prot.add_domain(1, 5, 'test_domain_5') # 5
    d5 = prot.get_domains_by_type('test_domain_5')[0]
    assert domain_tools.domain_overlap(d0, d5) is True
    assert domain_tools.domain_overlap(d1, d5) is True

    assert domain_tools.domain_overlap(d5, d0) is True
    assert domain_tools.domain_overlap(d5, d1) is True

def test_domain_overlap_fraction():    
    """

    """

    
    test_data_dir = shephard.get_data('test_data')
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    

    # TODO write a series of tests that test the domain_overlap_fraction
    # functionality
    
