import pytest
from shephard.exceptions import InterfaceException

import shephard
from shephard.interfaces.si_proteins import _ProteinsInterface as ProteinsInterface

test_data_dir = shephard.get_data('test_data')


def test_read_correctly():

    test_file_name = f"{test_data_dir}/TS1_small_proteins.tsv"

    # check we can read everything using default
    SI = ProteinsInterface(test_file_name)
    assert len(set(SI.data.keys())) == 3

    # check we can reduce parsing
    SI = ProteinsInterface(test_file_name, preauthorized_uids=['O00401'])
    assert len(set(SI.data.keys())) == 1

    # check we can read an empty list without throwing an error
    SI = ProteinsInterface(test_file_name, preauthorized_uids=[])
    assert len(set(SI.data.keys())) == 0



def test_bad_file():

    test_file_name = f"{test_data_dir}/TS1_sites.tsv"

    with pytest.raises(InterfaceException):
        SI = ProteinsInterface(test_file_name, skip_bad=False)

