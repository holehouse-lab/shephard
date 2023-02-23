import pytest

import shephard
from shephard.apis import uniprot  
from shephard.interfaces import si_domains
import numpy as np
from shephard.exceptions import ProteinException, DomainException
import random


test_data_dir = shephard.get_data('test_data')

# ....................................................................................................
#
#
def test_add_domains_file():


    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    domain_file = '%s/%s' % (test_data_dir, 'TS1_domains_idr.tsv')

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    si_domains.add_domains_from_file(P, domain_file)
    
    # this should fail because already added
    with pytest.raises(ProteinException):
        si_domains.add_domains_from_file(P, domain_file)

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    si_domains.add_domains_from_file(P, domain_file, autoname=True)

    print('')
    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    si_domains.add_domains_from_file(P, domain_file, autoname=False)

    # autoname allows 2 apparetly identical domain files to be added
    si_domains.add_domains_from_file(P, domain_file, autoname=True)

    # autoname allows 2 apparetly identical domain files to be added
    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    si_domains.add_domains_from_file(P, domain_file, autoname=False, skip_bad=True)

def test_add_domain_attribute():


    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    domain_file = '%s/%s' % (test_data_dir, 'TS1_domains_idr.tsv')

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    si_domains.add_domains_from_file(P, domain_file)
    prot = P.protein('O00401')
    domain = prot.domains[0]
    domain.add_attribute('test_attribute', 1)

    assert domain.attribute('test_attribute') == 1

    # this should fail
    with pytest.raises(DomainException):
        domain.add_attribute('test_attribute', 20)

    # because the operation above should have failed, this too should
    # have failed
    assert domain.attribute('test_attribute') == 1

    domain.add_attribute('test_attribute', 20, safe=False)
    assert domain.attribute('test_attribute') == 20

    assert len(domain.attributes) == 1
    domain.add_attribute('another_test_attribute', 'testval')
    assert len(domain.attributes) == 2


    with pytest.raises(DomainException):
        assert domain.attribute('does not exist') == 20

    # check this returns none
    assert domain.attribute('does not exist', safe=False) is None


    
def test_write_domain_with_attributes():

    # this setup was also tested in test_add_domain_attribute
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    domain_file = '%s/%s' % (test_data_dir, 'TS1_domains_idr.tsv')

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    si_domains.add_domains_from_file(P, domain_file)
    prot = P.protein('O00401')
    domain = prot.domains[0]
    domain.add_attribute('test_attribute_1', 1)
    domain.add_attribute('test_attribute_cat', 'cat')


def test_write_domains():

    
    TS1 = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir,'testset_1.fasta'))

    n_domains = 3

    uid2domain_info = {}
    for p in TS1:

        
        uid2domain_info[p.unique_ID] = []
        for idx in range(n_domains):
            s = random.randint(1,len(p))
            e = random.randint(1,len(p))
            if s > e:
                start = e
                end = s
            else:
                start = s
                end = e
                
            if start  == 0:
                start = start +1

            if end > len(p):
                end = len(p)-1
            

            uid2domain_info[p.unique_ID].append([start,end,idx])

            p.add_domain(start, end, f'test_domain_{idx}', attributes={'att1':'test'})

    si_domains.write_domains(TS1, 'output_test/test_domains.tsv')


    TS2 = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir,'testset_1.fasta'))
    si_domains.add_domains_from_file(TS2, 'output_test/test_domains.tsv')

    # check domain positions can be read in correctly
    for p in TS2:

        assert p.get_domains_by_type('test_domain_0')[0].start == uid2domain_info[p.unique_ID][0][0]
        assert p.get_domains_by_type('test_domain_1')[0].start == uid2domain_info[p.unique_ID][1][0]
        assert p.get_domains_by_type('test_domain_2')[0].start == uid2domain_info[p.unique_ID][2][0]


        assert p.get_domains_by_type('test_domain_0')[0].end == uid2domain_info[p.unique_ID][0][1]
        assert p.get_domains_by_type('test_domain_1')[0].end == uid2domain_info[p.unique_ID][1][1]
        assert p.get_domains_by_type('test_domain_2')[0].end == uid2domain_info[p.unique_ID][2][1]


def test_write_domains_with_domain_types():

    
    TS1 = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir,'testset_1.fasta'))

    n_domains = 2

    n_test_domains = 0
    uid2domain_info = {}
    for p in TS1:

        
        uid2domain_info[p.unique_ID] = []
        for idx in range(n_domains):
            s = random.randint(1,len(p))
            e = random.randint(1,len(p))
            if s > e:
                start = e
                end = s
            else:
                start = s
                end = e
                
            if start  == 0:
                start = start +1

            if end > len(p):
                end = len(p)-1
            

            uid2domain_info[p.unique_ID].append([start,end])

            p.add_domain(start, end, f'test_domain', attributes={'att1':'test'})
            n_test_domains = n_test_domains + 1

            p.add_domain(start, end, f'test_domain2', attributes={'att1':'test'})

    si_domains.write_domains(TS1, 'output_test/test_domains_types.tsv', domain_types=['test_domain'])

    TS2 = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir,'testset_1.fasta'))
    si_domains.add_domains_from_file(TS2, 'output_test/test_domains_types.tsv')

    assert len(TS2.domains) == n_test_domains

        
