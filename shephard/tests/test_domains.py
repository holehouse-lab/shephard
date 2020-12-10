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

names_list = ['O00401', 'O00470', 'O00472', 'O00499', 'O00629', 'O00712', 'O00716', 'O14786', 'Q9UJX3']

def test_bad_domain_add(TS1_domains2_sites_tracks):    
    """
    Note if any of these succeed the test suite will fail because it incrases
    """

    X = TS1_domains2_sites_tracks

    # add a domain where end is outside
    with pytest.raises(ProteinException):
        assert X.protein('O00401').add_domain(1,1000, 'fail')


    # add a domain where end is outside
    with pytest.raises(ProteinException):
        assert X.protein('O00401').add_domain(-1,20, 'fail')


def test_domain_add(TS1_domains2_sites_tracks):    
    """
    This function tests our ability to correctly add and then
    make reference to a domain
    
    """


    local_protein = TS1_domains2_sites_tracks.protein('O00629')
    
    domain_type = 'test_domain'
    local_protein.add_domain(1,10, domain_type)
    local_protein.add_domain(512,521, domain_type)

    
    # first check this works, if not, we all got done some real
    # bad codes in these parts
    assert local_protein.sequence[0:10] == 'MADNEKLDNQ'


    # ......................................................................
    # First try the N-terminal domain    
    #
    # this is just a sanity check because if THIS fails the
    # rest of the test definietly should not work

    local_domain = local_protein.get_domains_by_type(domain_type)[0]
    assert local_domain.sequence == 'MADNEKLDNQ'
    assert len(local_domain) == 10
    assert local_protein.get_sequence_region(1,10) == local_domain.sequence
    assert local_domain.start == 1
    assert local_domain.end == 10

    # these are more code coverage
    assert local_domain.protein.unique_ID == local_protein.unique_ID
    assert local_domain.domain_type == domain_type

    # assumes our domain name construction doesn't change. 
    assert local_domain.domain_name == "%s_%i_%i" % (domain_type, local_domain.start, local_domain.end)


    # note this should (A) work and (B) discard the N-terminal context
    assert local_protein.get_sequence_context(1,9) == local_domain.sequence

    # ......................................................................
    # Next check out C-terminal dmain
    #
    # next check out the
    local_domain = local_protein.get_domains_by_type(domain_type)[1]
    assert local_domain.sequence == 'ANVPTEGFQF'
    assert len(local_domain) == 10
    assert local_protein.get_sequence_region(512, 521) == local_domain.sequence
    assert local_domain.start == 512
    assert local_domain.end == 521

    # assumes our domain name construction doesn't change. 
    assert local_domain.domain_name == "%s_%i_%i" % (domain_type, local_domain.start, local_domain.end)

        
    # note this should (A) work and (B) discard the C-terminal context
    assert local_protein.get_sequence_context(521,9) == local_domain.sequence


def test_domain_attributes(TS1_domains2_sites_tracks):    
    
    # set up
    local_protein = TS1_domains2_sites_tracks.protein('O00629')
    domain_type = 'test_domain'
    local_protein.add_domain(20,200, domain_type, attributes={'test_1':1,'test_2':'yes' })
    local_domain = local_protein.get_domains_by_type(domain_type)[0]
    local_domain.add_attribute('added_attribute', np.array([1,2,3,4]))

    assert local_domain.attribute('test_1') == 1
    assert local_domain.attribute('test_2') == 'yes'
    assert local_domain.attribute('added_attribute').all() == np.array([1,2,3,4]).all()
    
    
    # no attribute called test_missing so should gracefully fail
    with pytest.raises(DomainException):
        assert local_domain.attribute('test_missing') == 2


    # no attribute called test_missing so should gracefully fail
    with pytest.raises(DomainException):
        local_domain.add_attribute('test_1', 2)

    # but we should be able to force override an attribute if safe=False
    local_domain.add_attribute('test_1', 2, safe=False)
    assert local_domain.attribute('test_1') == 2
        

    
