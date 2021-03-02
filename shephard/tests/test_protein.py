"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

import pytest
import shephard
from shephard.exceptions import ProteinException
from shephard.apis import uniprot  

names_list = ['O00401', 'O00470', 'O00472', 'O00499', 'O00629', 'O00712', 'O00716', 'O14786', 'Q9UJX3']

test_data_dir = shephard.get_data('test_data')

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

def test_get_track_values(TS1_domains2_sites_tracks):  

    p = TS1_domains2_sites_tracks.protein('O00401')

    assert len(p.get_track_values('pscore')) == len(TS1_domains2_sites_tracks.protein('O00401'))
    assert len(p.get_track_values('pscore', start=1)) == len(TS1_domains2_sites_tracks.protein('O00401'))
    assert len(p.get_track_values('pscore', end=len(p))) == len(TS1_domains2_sites_tracks.protein('O00401'))
    
    # check values work
    tvs = p.get_track_values('pscore', start=1, end=2)
    assert tvs[0] + 0.522 < 0.0001

    # get 26th position
    tvs = p.get_track_values('pscore', start=26, end=26)
    assert tvs[0] + 0.466 < 0.0001

    # should fail as start < 1
    with pytest.raises(ProteinException):
        p.get_track_values('pscore', start=0, end=26)

    # should fail as end passed end
    with pytest.raises(ProteinException):
        p.get_track_values('pscore', start=1, end=len(p)+1)
    


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

    expected = 'AOOAOOOOOOOOOAOOAOOAAAOOOOOOOAAOAAOOOOAOAOOOAAOAAOOOOOOAAOOOOOOAOOAAOOOOOOOAAAOAAOAOOOOAAAOOOAAOOAAAOOOOOAAOOAOOOOOOAOAOAOOOOOOOOAOOOAOOAAOOOOOOOOOOOOOOOOOOAOAOOAOAOOOOAOOOOAAOOOAOOAOOOOOOOOOOOOOOOAOOOOAOOOOOAOOAOOAOAOOOOOAOAOOAOOOAOOAAOAOOAOOOOAOOOOOOOAAAOAAOOOOOAOOAOOOAOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOAOAOOOOOOOAAOOOOOOAOOOOOOOOOOOOOOAAOAOOAOOOOOOOOOOOOOOOOOOOAOOOOOOOAOOOOOOOOOAAOOAOOOOOAOOAOOOOOOAOOOOOOOAAOOAOOOAOAOOAOOOOOOOOOOOOOOOOAAOOAAOAAOOOOOOAOOOOOOOOOOOOOOAOOOOOAOO'

    assert "".join(protein.track('hydrophobes2').symbols) == expected
    assert "".join(protein.track('hydrophobes3').symbols) == expected
    
    # assert that sans Safe=False this fails
    with pytest.raises(ProteinException):
        assert protein.build_track_symbols_from_sequence('hydrophobes2', convert_sequence_symbols, local_residues)

    # now assert that the same function with safe=False will run ok!
    assert protein.build_track_symbols_from_sequence('hydrophobes2', convert_sequence_symbols, local_residues, safe=False) == None


def test_getter_properies():
    TS1 = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir,'testset_1.fasta'))

    prot = TS1.protein('O00470')

    assert len(TS1.protein('O00470')) == 390
    
    assert prot.name == "sp|O00470|MEIS1_HUMAN Homeobox protein Meis1 OS=Homo sapiens OX=9606 GN=MEIS1 PE=1 SV=1"
    assert str(prot.proteome) == "[Proteome]: Sequence dataset with 9 protein records"
    assert prot.unique_ID == 'O00470'
    assert prot.get_residue(1) == 'M'
    assert prot.get_residue(2) == 'A'
    assert prot.get_residue(390) == 'M'
    
    with pytest.raises(ProteinException):
        assert prot.get_residue(0) == 'M'
        
    #assert protein.get_sequence = '
    assert prot.sequence == 'MAQRYDDLPHYGGMDGVGIPSTMYGDPHAARSMQPVHHLNHGPPLHSHQYPHTAHTNAMAPSMGSSVNDALKRDKDAIYGHPLFPLLALIFEKCELATCTPREPGVAGGDVCSSESFNEDIAVFAKQIRAEKPLFSSNPELDNLMIQAIQVLRFHLLELEKVHELCDNFCHRYISCLKGKMPIDLVIDDREGGSKSDSEDITRSANLTDQPSWNRDHDDTASTRSGGTPGPSSGGHTSHSGDNSSEQGDGLDNSVASPSTGDDDDPDKDKKRHKKRGIFPKVATNIMRAWLFQHLTHPYPSEEQKKQLAQDTGLTILQVNNWFINARRRIVQPMIDQSNRAVSQGTPYNPDGQPMGGFVMDGQQHMGIRAPGPMSGMGMNMGMEGQWHYM'

    assert prot.get_sequence_region(1,5) == 'MAQRY'
    assert prot.get_sequence_region(390,390) == 'M'

    with pytest.raises(ProteinException):    
        assert prot.get_sequence_region(390,391) == 'M'

    assert prot.get_sequence_region(380,390) == 'NMGMEGQWHYM'

    assert prot.get_sequence_context(1,10) == 'MAQRYDDLPHY'
    assert prot.get_sequence_context(10,3) == 'DLPHYGG'
    assert prot.get_sequence_context(390,5) == 'GQWHYM'
    assert prot.get_sequence_context(50,450) == 'MAQRYDDLPHYGGMDGVGIPSTMYGDPHAARSMQPVHHLNHGPPLHSHQYPHTAHTNAMAPSMGSSVNDALKRDKDAIYGHPLFPLLALIFEKCELATCTPREPGVAGGDVCSSESFNEDIAVFAKQIRAEKPLFSSNPELDNLMIQAIQVLRFHLLELEKVHELCDNFCHRYISCLKGKMPIDLVIDDREGGSKSDSEDITRSANLTDQPSWNRDHDDTASTRSGGTPGPSSGGHTSHSGDNSSEQGDGLDNSVASPSTGDDDDPDKDKKRHKKRGIFPKVATNIMRAWLFQHLTHPYPSEEQKKQLAQDTGLTILQVNNWFINARRRIVQPMIDQSNRAVSQGTPYNPDGQPMGGFVMDGQQHMGIRAPGPMSGMGMNMGMEGQWHYM'

    with pytest.raises(ProteinException):    
        prot.get_sequence_context(400, 3)

    assert prot.check_sequence_is_valid() is True

    assert len(prot.attributes) == 0

    assert prot.attribute('TEST', safe=False) is None

    assert len(prot.tracks) == 0
    assert len(prot.track_names) == 0
    assert prot.track('TEST', safe=False) is None
    #assert prot.get_track_values('TEST', safe=False) is None
    #assert prot.get_track_symbols('TEST', safe=False) is None

    


def test_add_site():

    # build new proteome
    TS1 = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir,'testset_1.fasta'))
