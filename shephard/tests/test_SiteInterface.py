import pytest
from shephard.exceptions import InterfaceException



import shephard
from shephard.interfaces.si_sites import _SitesInterface as SitesInterface

test_data_dir = shephard.get_data('test_data')


def test_read_correctly():

    test_file_name = f"{test_data_dir}/TS1_sites.tsv"

    # check we can read everything using default
    SI = SitesInterface(test_file_name)
    assert len(set(SI.data.keys())) == 9

    # check we can reduce parsing
    SI = SitesInterface(test_file_name, preauthorized_uids=['O00401'])
    assert len(set(SI.data.keys())) == 1

    # check we can read an empty list without throwing an error
    SI = SitesInterface(test_file_name, preauthorized_uids=[])
    assert len(set(SI.data.keys())) == 0



def test_bad_file():

    test_file_name = f"{test_data_dir}/TS1_small_proteins.tsv"

    with pytest.raises(InterfaceException):
        SI = SitesInterface(test_file_name, skip_bad=False)
    
