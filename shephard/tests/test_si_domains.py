import pytest

import shephard
from shephard.apis import uniprot  
from shephard import interfaces
import numpy as np
from shephard.exceptions import ProteinException, DomainException




# ....................................................................................................
#
#
def test_add_domains_file():

    test_data_dir = shephard.get_data('test_data')
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    domain_file = '%s/%s' % (test_data_dir, 'TS1_domains_idr.tsv')

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_domains.add_domains_from_file(P, domain_file)
    
    # this should fail because already added
    with pytest.raises(ProteinException):
        interfaces.si_domains.add_domains_from_file(P, domain_file)

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_domains.add_domains_from_file(P, domain_file, autoname=True)

    print('')
    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_domains.add_domains_from_file(P, domain_file, autoname=False)

    # autoname allows 2 apparetly identical domain files to be added
    interfaces.si_domains.add_domains_from_file(P, domain_file, autoname=True)

    # autoname allows 2 apparetly identical domain files to be added
    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_domains.add_domains_from_file(P, domain_file, autoname=False, skip_bad=True)

def test_add_domain_attribute():

    test_data_dir = shephard.get_data('test_data')
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    domain_file = '%s/%s' % (test_data_dir, 'TS1_domains_idr.tsv')

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_domains.add_domains_from_file(P, domain_file)
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
    test_data_dir = shephard.get_data('test_data')
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    domain_file = '%s/%s' % (test_data_dir, 'TS1_domains_idr.tsv')

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_domains.add_domains_from_file(P, domain_file)
    prot = P.protein('O00401')
    domain = prot.domains[0]
    domain.add_attribute('test_attribute_1', 1)
    domain.add_attribute('test_attribute_cat', 'cat')
