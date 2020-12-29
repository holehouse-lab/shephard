import pytest

import shephard
from shephard.apis import uniprot  
from shephard import interfaces
import numpy as np
from shephard.exceptions import ProteinException 




# ....................................................................................................
#
#
def test_add_track_from_file():

    test_data_dir = shephard.get_data('test_data')
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    track_file = '%s/%s' % (test_data_dir, 'TS1_tracks_pscore.tsv')
    
    
    def testfunct(P):
        assert P.protein('O00401').track('pscore').values_region(50,50)[0] == -0.971
        assert P.protein('O00470').track('pscore').values_region(50,50)[0] == 0.324
        assert P.protein('O00472').track('pscore').values_region(50,50)[0] == -0.661
        assert P.protein('O00499').track('pscore').values_region(50,50)[0] == 0.205

    
        assert np.sum(P.protein('O00499').track('pscore').values_region(500,593)) == 43.647999999999996

        with pytest.raises(ProteinException):
            assert 1 == P.protein('O00499').track('pscore').values_region(594, 594)[0]

        
    # test various permutations of reading in (all should work)

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P, track_file, mode='values', safe=False, verbose=True, skip_bad=False)
    testfunct(P)

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P, track_file, mode='values', verbose=True, skip_bad=False)
    testfunct(P)

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P, track_file, mode='values', safe=False, skip_bad=False)
    testfunct(P)

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P, track_file, mode='values', safe=False, skip_bad=True)
    testfunct(P)

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P, track_file, mode='values', safe=True, skip_bad=True)
    testfunct(P)

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P, track_file, mode='values', safe=False, verbose=True)
    testfunct(P)

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P, track_file, mode='values', verbose=True)
    testfunct(P)

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P, track_file, mode='values', skip_bad=True)
    testfunct(P)

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P, track_file, mode='values', safe=False)
    testfunct(P)


# ....................................................................................................
#
#
def test_add_track_from_dictionary():

    test_data_dir = shephard.get_data('test_data')
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    
    P = uniprot.uniprot_fasta_to_proteome(fasta_file)

    
    random_vals = np.random.rand(len(P.protein('O00401')))
    track_dict = {'O00401': [['random_vals', random_vals]]}
    interfaces.si_tracks.add_tracks_from_dictionary(P, track_dict, mode='values')
    assert len(P.protein('O00401')) == len(P.protein('O00401').track('random_vals').values)

    newstring = ''
    for i in P.protein('O00401').sequence:

        if i in ['S','T','Y']:
            newstring = newstring + "P"
        else:
            newstring = newstring + "-"

    track_dict = {'O00401': [['phosres', newstring]]}
    interfaces.si_tracks.add_tracks_from_dictionary(P, track_dict, mode='symbols')


    with pytest.raises(TypeError):
        assert len(P.protein('O00401')) == len(P.protein('O00401').track('phosres').values)

    assert "".join(P.protein('O00401').track('phosres').symbols_region(1,10)) == "-PP-------"
    
    with pytest.raises(ProteinException):
        assert "".join(P.protein('O00401').track('phosres').symbols_region(0,10)) == "-PP-------"

    with pytest.raises(ProteinException):
        local_p = P.protein('O00401')
        print(local_p.track('phosres').symbols_region(len(local_p), len(local_p)+10))



# ....................................................................................................
#
#
def test_write_track_part1():

    test_data_dir = shephard.get_data('test_data')
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    track_file = '%s/%s' % (test_data_dir, 'TS1_tracks_pscore.tsv')
    
    # write
    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P, track_file, mode='values')
    interfaces.si_tracks.write_track(P, 'output_test/out_tracks.tsv', 'pscore')

    P2 = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P2, 'output_test/out_tracks.tsv', mode='values')

    for idx in P.proteins:
        protein = P.protein(idx)
        protein2 = P2.protein(idx)
        
        assert np.sum(protein.track('pscore').values) == np.sum(protein2.track('pscore').values)

# ....................................................................................................
#
#
def test_write_track_part2():

    test_data_dir = shephard.get_data('test_data')
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    pscore_track_file = '%s/%s' % (test_data_dir, 'TS1_tracks_pscore.tsv')


    # build proteome
    P = uniprot.uniprot_fasta_to_proteome(fasta_file)

    # create an additional values track
    special_ID = 'O00401'
    random_vals = np.random.rand(len(P.protein(special_ID)))
    track_name = 'random_vals'
    outname = 'output_test/track_%s.tsv'%(track_name)

    track_dict = {special_ID: [[track_name, random_vals]]}
    interfaces.si_tracks.add_tracks_from_dictionary(P, track_dict, mode='values')
    interfaces.si_tracks.write_track(P, outname,track_name)

    P2 = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P2, pscore_track_file, mode='values')
    interfaces.si_tracks.add_tracks_from_file(P2, outname, mode='values')

    for protein in P2:
        if protein.unique_ID == special_ID:
            assert len(protein.tracks) == 2                
        else:
            assert len(protein.tracks) == 1       


def test_write_track_part3():

    test_data_dir = shephard.get_data('test_data')
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    pscore_track_file = '%s/%s' % (test_data_dir, 'TS1_tracks_pscore.tsv')


    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P, pscore_track_file, mode='values')

    # ..................................................................
    ## First create and write a random vals track
    # create an additional values track
    special_ID = 'O00401'
    random_vals = np.random.rand(len(P.protein(special_ID)))
    track_name_rv = 'random_vals'
    outname_rv = 'output_test/track_%s.tsv'%(track_name_rv)

    track_dict = {special_ID: [[track_name_rv, random_vals]]}
    interfaces.si_tracks.add_tracks_from_dictionary(P, track_dict, mode='values')
    interfaces.si_tracks.write_track(P, outname_rv, track_name_rv)


    track_name = 'phospho_track'
    outname = 'output_test/track_%s.tsv'%(track_name)


    # create an additional symbols tracks and write track
    newstring = ''
    for i in P.protein(special_ID).sequence:

        if i in ['S','T','Y']:
            newstring = newstring + "P"
        else:
            newstring = newstring + "-"

    track_dict = {special_ID: [[track_name, newstring]]}
    interfaces.si_tracks.add_tracks_from_dictionary(P, track_dict, mode='symbols')    
    interfaces.si_tracks.write_track(P, outname, track_name)

    P3 = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P3, outname, mode='symbols')
    interfaces.si_tracks.add_tracks_from_file(P3, outname_rv, mode='values')
    
    for protein in P3:
        if protein.unique_ID == special_ID:
            assert len(protein.tracks) == 2               
        else:
            assert len(protein.tracks) == 0      

    


        
    



