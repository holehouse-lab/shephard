"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

import pytest
from shephard.exceptions import ProteinException

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
    Note if any of these succeed the test suite will fail because it incrases
    """

    print(TS1_domains2_sites_tracks.protein('O00629').sequence)

    local_protein = TS1_domains2_sites_tracks.protein('O00629')
    
    domain_type = 'test_domain'
    local_protein.add_domain(1,10, domain_type)

    # this is just a sanity check because if THIS fails the
    # rest of the test definietly should not work
    assert local_protein.sequence[0:10] == 'MADNEKLDNQ'
    assert local_protein.get_domains_by_type(domain_type)[0].sequence == 'MADNEKLDNQ'
    
