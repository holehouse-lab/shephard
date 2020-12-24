import pytest
import shephard
from shephard import proteome
from shephard.apis import fasta
from shephard.exceptions import ProteomeException, APIException


def test_fasta_to_proteome_part_1():

    test_data_dir = shephard.get_data('test_data')
    print(test_data_dir)

    P = fasta.fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'))
    assert len(P.protein('1')) == 390
    assert len(P) == 9

    test_UID = 0
    for i in P.proteins:
        assert i == str(test_UID)
        test_UID = test_UID + 1

    P = fasta.fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'), proteome=P)
    assert len(P) == 18

    test_UID = 0
    for i in P.proteins:
        assert i == str(test_UID)
        test_UID = test_UID + 1


    ##
    ## This block checks that removing a protein from the integer-indexed added proteins
    ## really removes it and that adding new proteins in does correctly start counting 
    ## in the right place
    P.remove_protein(10)
    P = fasta.fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'), proteome=P)
        
    test_UID = 0
    for i in P.proteins:
        if test_UID == 10:
            with pytest.raises(ProteomeException):
                assert P.protein(test_UID)
            test_UID = test_UID + 1

        assert i == str(test_UID)
        test_UID = test_UID + 1

    
def test_fasta_to_proteome_part_2():
    def header_parser(s):
        return s.split('|')[1]

    test_data_dir = shephard.get_data('test_data')
    print(test_data_dir)

    # check we've read in the proteome OK
    P = fasta.fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'), build_unique_ID=header_parser)
    assert len(P) == 9

    # this should trigger an exception because we're adding in duplicates
    with pytest.raises(ProteomeException):
        P = fasta.fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'), build_unique_ID=header_parser, proteome=P)

    # this should not...
    P = fasta.fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'), proteome=P)
    assert len(P) == 18

    # this also should not but should NOT add new sequences in
    P = fasta.fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'), build_unique_ID=header_parser, proteome=P, force_overwrite=True)
    assert len(P) == 18

    P = fasta.fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'), proteome=P)
    assert len(P) == 27

    P = fasta.fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'), proteome=P, use_header_as_unique_ID=True)
    assert len(P) == 36

    # this SHOULD trigger an exception because we shouldn't
    with pytest.raises(APIException):
        P = fasta.fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'), proteome=P, use_header_as_unique_ID=True, build_unique_ID=header_parser)

    expected_protein_uids = ['O00401', 'O00470', 'O00472', 'O00499', 'O00629', 'O00712', 'O00716', 'O14786', 'Q9UJX3', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', 'sp|O00401|WASL_HUMAN Neural Wiskott-Aldrich syndrome protein OS=Homo sapiens OX=9606 GN=WASL PE=1 SV=2', 'sp|O00470|MEIS1_HUMAN Homeobox protein Meis1 OS=Homo sapiens OX=9606 GN=MEIS1 PE=1 SV=1', 'sp|O00472|ELL2_HUMAN RNA polymerase II elongation factor ELL2 OS=Homo sapiens OX=9606 GN=ELL2 PE=1 SV=2', 'sp|O00499|BIN1_HUMAN Myc box-dependent-interacting protein 1 OS=Homo sapiens OX=9606 GN=BIN1 PE=1 SV=1', 'sp|O00629|IMA3_HUMAN Importin subunit alpha-3 OS=Homo sapiens OX=9606 GN=KPNA4 PE=1 SV=1', 'sp|O00712|NFIB_HUMAN Nuclear factor 1 B-type OS=Homo sapiens OX=9606 GN=NFIB PE=1 SV=2', 'sp|O00716|E2F3_HUMAN Transcription factor E2F3 OS=Homo sapiens OX=9606 GN=E2F3 PE=1 SV=1', 'sp|O14786|NRP1_HUMAN Neuropilin-1 OS=Homo sapiens OX=9606 GN=NRP1 PE=1 SV=3', 'sp|Q9UJX3|APC7_HUMAN Anaphase-promoting complex subunit 7 OS=Homo sapiens OX=9606 GN=ANAPC7 PE=1 SV=4']

    # CHECK all the unique IDs expected can be read in
    for i in expected_protein_uids:
        P.protein(i) 

    

    


    
    
    
