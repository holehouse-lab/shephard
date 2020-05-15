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


