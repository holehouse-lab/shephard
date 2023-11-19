
import pytest
from shephard.exceptions import InterfaceException

import shephard
from shephard.interfaces.si_tracks import _TracksInterface as TracksInterface

test_data_dir = shephard.get_data('test_data')


def test_read_correctly():

    test_file_name = f"{test_data_dir}/TS1_tracks_pscore.tsv"

    # check we can read everything using default
    SI = TracksInterface(test_file_name, mode='values')
    assert len(set(SI.data.keys())) == 9

    # check we can reduce parsing
    SI = TracksInterface(test_file_name, mode='values', preauthorized_uids=['O00401'])
    assert len(set(SI.data.keys())) == 1

    # check we can read an empty list without throwing an error
    SI = TracksInterface(test_file_name, mode='values', preauthorized_uids=[])
    assert len(set(SI.data.keys())) == 0



    

def test_mode_keyword():

    test_file_name = f"{test_data_dir}/TS1_tracks_pscore.tsv"

    # check we can read everything using default
    SI = TracksInterface(test_file_name, mode='values')
    SI.data['O00470'][0]['track_data'][0] == 1.504

    SI = TracksInterface(test_file_name, mode='symbols')
    SI.data['O00470'][0]['track_data'][0] == "1.504"
    

    with pytest.raises(InterfaceException):
        SI = TracksInterface(test_file_name, mode='other', skip_bad=False)
        

