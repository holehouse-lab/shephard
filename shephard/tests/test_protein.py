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

def test_sequence_lookup(TS1_domains2_sites_tracks):    
    

    natural_sequence = 'MSSVQQQPPPPRRVTNVGSLLLTPQENESLFTFLGKKCVTMSSAVVQLYAADRNCMWSKKCSGVACLVKDNPQRSYFLRIFDIKDGKLLWEQELYNNFVYNSPRGYFHTFAGDTCQVALNFANEEEAKKFRKAVTDLLGRRQRKSEKRRDPPNGPNLPMATVDIKNPEITTNRFYGPQVNNISHTKEKKKGKAKKKRLTKADIGTPSNFQHIGHVGWDPNTGFDLNNLDPELKNLFDMCGISEAQLKDRETSKVIYDFIEKTGGVEAVKNELRRQAPPPPPPSRGGPPPPPPPPHNSGPPPPPARGRGAPPPPPSRAPTAAPPPPPPSRPSVAVPPPPPNRMYPPPPPALPSSAPSGPPPPPPSVLGVGPVAPPPPPPPPPPPGPPPPPGLPSDGDHQVPTTAGNKAALLDQIREGAQLKKVEQNSRPVSCSGRDALLDQIRQGIQLKSVADGQESTPPTPAPTSGIVGALMEVMQKRSKAIHSSDEDEDEDDEEDFEDDDEWED'


    X = TS1_domains2_sites_tracks

    # test correct indexing works 
    assert X.protein('O00401').get_sequence_region(1,10) == 'MSSVQQQPPP'

    assert X.protein('O00401').sequence == natural_sequence
    

def test_build_track_values_from_sequence(TS1_domains2_sites_tracks):  

    def convert_sequence_values(sequence, targetresidues):
        """
        Function converts sequence to 1s and 0s based on targetresidues 
        """
        valuetrack=[]
    
        #iterate sequence and build track 
        for i in sequence:
            if str(i) in targetresidues['target_residues']:
                valuetrack.append(1)
            else:
                valuetrack.append(0)      

        return valuetrack

    local_residues = {'target_residues':['L','Y','M','W','F','I','V'] }

    protein = TS1_domains2_sites_tracks.protein('O00401')
    
    protein.build_track_values_from_sequence('hydrophobes', convert_sequence_values, input_dictionary=local_residues)

    expected = "[1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]"


    assert str(protein.track('hydrophobes').values) == expected
    
    # assert that sans Safe=False this fails
    with pytest.raises(ProteinException):
        assert protein.build_track_values_from_sequence('hydrophobes', convert_sequence_values, local_residues)

    # now assert that the same function with safe=False will run ok!
    assert protein.build_track_values_from_sequence('hydrophobes', convert_sequence_values, local_residues, safe=False) == None



def test_build_track_symbols_from_sequence(TS1_domains2_sites_tracks):  

    def convert_sequence_symbols(sequence, targetresidues):
        """
        Function converts sequence to 1s and 0s based on targetresidues 
        """
        valuetrack=[]
    
        #iterate sequence and build track 
        for i in sequence:
            if str(i) in targetresidues['target_residues']:
                valuetrack.append("A")
            else:
                valuetrack.append("O")      

        return valuetrack


    def convert_sequence_symbols_no_params(sequence):
        """
        Function converts sequence to 1s and 0s based on targetresidues 
        """
        valuetrack=[]
    
        #iterate sequence and build track 
        for i in sequence:
            if str(i) in ['L','Y','M','W','F','I','V']:
                valuetrack.append("A")
            else:
                valuetrack.append("O")      

        return valuetrack
    
    # --------------------------------------------------------------------------------


    local_residues = {'target_residues':['L','Y','M','W','F','I','V'] }
    
    protein = TS1_domains2_sites_tracks.protein('O00401')
    
    protein.build_track_symbols_from_sequence('hydrophobes2', convert_sequence_symbols, local_residues)
    protein.build_track_symbols_from_sequence('hydrophobes3', convert_sequence_symbols_no_params)

    print(protein.sequence)
    expected = 'AOOAOOOOOOOOOAOOAOOAAAOOOOOOOAAOAAOOOOAOAOOOAAOAAOOOOOOAAOOOOOOAOOAAOOOOOOOAAAOAAOAOOOOAAAOOOAAOOAAAOOOOOAAOOAOOOOOOAOAOAOOOOOOOOAOOOAOOAAOOOOOOOOOOOOOOOOOOAOAOOAOAOOOOAOOOOAAOOOAOOAOOOOOOOOOOOOOOOAOOOOAOOOOOAOOAOOAOAOOOOOAOAOOAOOOAOOAAOAOOAOOOOAOOOOOOOAAAOAAOOOOOAOOAOOOAOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOAOAOOOOOOOAAOOOOOOAOOOOOOOOOOOOOOAAOAOOAOOOOOOOOOOOOOOOOOOOAOOOOOOOAOOOOOOOOOAAOOAOOOOOAOOAOOOOOOAOOOOOOOAAOOAOOOAOAOOAOOOOOOOOOOOOOOOOAAOOAAOAAOOOOOOAOOOOOOOOOOOOOOAOOOOOAOO'

    assert "".join(protein.track('hydrophobes2').symbols) == expected
    assert "".join(protein.track('hydrophobes3').symbols) == expected
    
    # assert that sans Safe=False this fails
    with pytest.raises(ProteinException):
        assert protein.build_track_symbols_from_sequence('hydrophobes2', convert_sequence_symbols, local_residues)

    # now assert that the same function with safe=False will run ok!
    assert protein.build_track_symbols_from_sequence('hydrophobes2', convert_sequence_symbols, local_residues, safe=False) == None

