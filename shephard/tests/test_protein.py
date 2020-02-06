import pytest

names_list = ['O00401', 'O00470', 'O00472', 'O00499', 'O00629', 'O00712', 'O00716', 'O14786', 'Q9UJX3']

def test_sequence_lookup(TS1_domains2_sites_tracks):    
    

    natural_sequence = 'MSSVQQQPPPPRRVTNVGSLLLTPQENESLFTFLGKKCVTMSSAVVQLYAADRNCMWSKKCSGVACLVKDNPQRSYFLRIFDIKDGKLLWEQELYNNFVYNSPRGYFHTFAGDTCQVALNFANEEEAKKFRKAVTDLLGRRQRKSEKRRDPPNGPNLPMATVDIKNPEITTNRFYGPQVNNISHTKEKKKGKAKKKRLTKADIGTPSNFQHIGHVGWDPNTGFDLNNLDPELKNLFDMCGISEAQLKDRETSKVIYDFIEKTGGVEAVKNELRRQAPPPPPPSRGGPPPPPPPPHNSGPPPPPARGRGAPPPPPSRAPTAAPPPPPPSRPSVAVPPPPPNRMYPPPPPALPSSAPSGPPPPPPSVLGVGPVAPPPPPPPPPPPGPPPPPGLPSDGDHQVPTTAGNKAALLDQIREGAQLKKVEQNSRPVSCSGRDALLDQIRQGIQLKSVADGQESTPPTPAPTSGIVGALMEVMQKRSKAIHSSDEDEDEDDEEDFEDDDEWED'


    X = TS1_domains2_sites_tracks

    # test correct indexing works 
    assert X.protein('O00401').get_sequence_region(1,10) == 'MSSVQQQPPP'

    assert X.protein('O00401').sequence == natural_sequence
    
    """
    with pytest.raises(Exception) as e:
        assert X.protein('O00401').sequence_sliceable[0]
    assert str(e.value) == "For your protection complex slicing is disabled"

    with pytest.raises(Exception) as e:
        assert X.protein('O00401').sequence_sliceable[0:3:10]
    assert str(e.value) == "For your protection complex slicing is disabled"

    with pytest.raises(Exception) as e:
        assert X.protein('O00401').sequence_sliceable[-1]
    assert str(e.value) == "For your protection complex slicing is disabled"
    """
