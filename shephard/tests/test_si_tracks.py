import pytest

import shephard
from shephard.apis import uniprot  
from shephard import interfaces
import numpy as np
from shephard.exceptions import ProteinException, InterfaceException
import os
import random



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
    track_dict = {'O00401': [{'track_name':'random_vals', 'track_data':random_vals}]}
    interfaces.si_tracks.add_tracks_from_dictionary(P, track_dict, mode='values')
    assert len(P.protein('O00401')) == len(P.protein('O00401').track('random_vals').values)

    newstring = ''
    for i in P.protein('O00401').sequence:

        if i in ['S','T','Y']:
            newstring = newstring + "P"
        else:
            newstring = newstring + "-"

    track_dict = {'O00401': [{'track_name':'phosres', 'track_data':newstring}]}
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

    fn = 'output_test/out_tracks.tsv'
    interfaces.si_tracks.write_tracks(P, fn, 'pscore')

    b = os.path.getsize(fn)
    assert b == 32531 # file size in bytes if correctly written

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

    track_dict = {special_ID: [{'track_name':track_name, 'track_data':random_vals}]}
    interfaces.si_tracks.add_tracks_from_dictionary(P, track_dict, mode='values')
    interfaces.si_tracks.write_tracks(P, outname,track_name)

    b = os.path.getsize(outname)
    assert b == 3050 # file size in bytes if correctly written

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

    track_dict = {special_ID: [{'track_name':track_name_rv, 'track_data':random_vals}]}
    interfaces.si_tracks.add_tracks_from_dictionary(P, track_dict, mode='values')
    interfaces.si_tracks.write_tracks(P, outname_rv, track_name_rv)


    track_name = 'phospho_track'
    outname = 'output_test/track_%s.tsv'%(track_name)


    # create an additional symbols tracks and write track
    newstring = ''
    for i in P.protein(special_ID).sequence:

        if i in ['S','T','Y']:
            newstring = newstring + "P"
        else:
            newstring = newstring + "-"

    track_dict = {special_ID: [{'track_name':track_name, 'track_data':newstring}]}
    interfaces.si_tracks.add_tracks_from_dictionary(P, track_dict, mode='symbols')    
    interfaces.si_tracks.write_tracks(P, outname, track_name)

    b = os.path.getsize(outname)
    assert b == 1032 # file size in bytes if correctly written

    P3 = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P3, outname, mode='symbols')
    interfaces.si_tracks.add_tracks_from_file(P3, outname_rv, mode='values')
    
    for protein in P3:
        if protein.unique_ID == special_ID:
            assert len(protein.tracks) == 2               
        else:
            assert len(protein.tracks) == 0      


def test_write_tracks_separate_files():

    ##
    ## This code block just creates a proteome where we have 2 different
    ## typs of tracks
    test_data_dir = shephard.get_data('test_data')
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    pscore_track_file = '%s/%s' % (test_data_dir, 'TS1_tracks_pscore.tsv')

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P, pscore_track_file, mode='values')

    newstring = ''
    for i in P.protein('O00401').sequence:

        if i in ['S','T','Y']:
            newstring = newstring + "P"
        else:
            newstring = newstring + "-"
    track_dict = {'O00401': [{'track_name':'phosres', 'track_data':newstring}]}
    interfaces.si_tracks.add_tracks_from_dictionary(P, track_dict, mode='symbols')

    interfaces.si_tracks.write_all_tracks_separate_files(P, 'output_test')



def test_write_tracks_separate_files():

    ##
    ## This code block just creates a proteome where we have 2 different
    ## typs of tracks
    test_data_dir = shephard.get_data('test_data')
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    pscore_track_file = '%s/%s' % (test_data_dir, 'TS1_tracks_pscore.tsv')

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P, pscore_track_file, mode='values')

    newstring = ''
    for i in P.protein('O00401').sequence:

        if i in ['S','T','Y']:
            newstring = newstring + "P"
        else:
            newstring = newstring + "-"
    track_dict = {'O00401': [{'track_name':'phosres', 'track_data':newstring}]}
    interfaces.si_tracks.add_tracks_from_dictionary(P, track_dict, mode='symbols')

    interfaces.si_tracks.write_all_tracks_separate_files(P, 'output_test')

    P2 = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P2, 'output_test/shephard_track_pscore.tsv', mode='values')

    for idx in P2.proteins:
        assert P2.protein(idx).track('pscore').values == P.protein(idx).track('pscore').values


    P2 = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P2, 'output_test/shephard_track_phosres.tsv', mode='symbols')

    for idx in ['O00401']:
        for position in range(1,len(P2.protein(idx))+1):
            POS1 = P.protein(idx).track('phosres').symbols_region(position,position)
            POS2 = P2.protein(idx).track('phosres').symbols_region(position,position)
            assert POS1 == POS2

        
def test_write_tracks_single_file_symbols():

    test_data_dir = shephard.get_data('test_data')
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    
    P = uniprot.uniprot_fasta_to_proteome(fasta_file)

    # build phosphostring
    newstring = ''
    for i in P.protein('O00401').sequence:

        if i in ['S','T','Y']:
            newstring = newstring + "P"
        else:
            newstring = newstring + "-"
    track_dict = {'O00401': [{'track_name':'phosres', 'track_data':newstring}]}
    interfaces.si_tracks.add_tracks_from_dictionary(P, track_dict, mode='symbols')

    # build phosphostring
    newstring = ''
    for i in P.protein('O00472').sequence:

        if i in ['Y','F','W']:
            newstring = newstring + "A"
        else:
            newstring = newstring + "-"
    track_dict = {'O00472': [{'track_name':'aro', 'track_data':newstring}]}
    
    interfaces.si_tracks.add_tracks_from_dictionary(P, track_dict, mode='symbols')
    interfaces.si_tracks.write_all_symbols_tracks_single_file(P, 'output_test/all_symbols_test.tsv')
    
    P2 = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P2, 'output_test/all_symbols_test.tsv', mode='symbols')

    idx = 'O00401'
    for position in range(1,len(P2.protein(idx))+1):
        POS1 = P.protein(idx).track('phosres').symbols_region(position,position)
        POS2 = P2.protein(idx).track('phosres').symbols_region(position,position)
        assert POS1 == POS2


    idx = 'O00472'
    for position in range(1,len(P2.protein(idx))+1):
        POS1 = P.protein(idx).track('aro').symbols_region(position,position)
        POS2 = P2.protein(idx).track('aro').symbols_region(position,position)
        assert POS1 == POS2


def test_write_tracks_single_file_values():

    test_data_dir = shephard.get_data('test_data')
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    
    P = uniprot.uniprot_fasta_to_proteome(fasta_file)

    # build phosphostring

    r1 = random.randint(1,20)
    r2 = random.random()
    for protein in P:

        new_vals = []
        for i in protein.sequence:

            if i in ['S','T','Y']:
                new_vals.append(r1)
            else:
                new_vals.append(r2)

        protein.add_track('test_track', values=new_vals)

    # write all the tracks out
    interfaces.si_tracks.write_all_values_tracks_single_file(P, 'output_test/all_values_test.tsv')
    
    # read the tracks into a new proteome
    P2 = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P2, 'output_test/all_values_test.tsv', mode='values')

    # check everything worked ok
    for protein in P2:
        print(protein)
        for i in range(1,len(protein)+1):

            if protein.residue(i) in ['S','T','Y']:
                assert protein.track('test_track').value(i) == r1
            else:
                assert protein.track('test_track').value(i) == np.round(r2,3)



def test_write_values_tracks_from_list():

    test_data_dir = shephard.get_data('test_data')
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    
    P = uniprot.uniprot_fasta_to_proteome(fasta_file)

    # build phosphostring

    r1 = random.randint(1,20)
    r2 = random.random()
    for protein in P:

        new_vals = []
        for i in protein.sequence:

            if i in ['S','T','Y']:
                new_vals.append(r1)
            else:
                new_vals.append(r2)

        protein.add_track('test_track', values=new_vals)

    track_list = []
    for protein in P:
        track_list.append(protein.track('test_track'))

    # write all the tracks out
    interfaces.si_tracks.write_tracks_from_list(track_list, 'output_test/all_values_test_from_list.tsv')

    # read the tracks into a new proteome
    P2 = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P2, 'output_test/all_values_test_from_list.tsv', mode='values')

    # check everything worked ok
    for protein in P2:
        print(protein)
        for i in range(1,len(protein)+1):

            if protein.residue(i) in ['S','T','Y']:
                assert protein.track('test_track').value(i) == r1
            else:
                assert protein.track('test_track').value(i) == np.round(r2,3)



def test_write_symbols_tracks_from_list():

    test_data_dir = shephard.get_data('test_data')
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    
    P = uniprot.uniprot_fasta_to_proteome(fasta_file)

    # build phosphostring

    r1 = random.randint(1,20)
    r2 = random.random()
    for protein in P:

        new_symbols =''
        for i in protein.sequence:

            if i in ['S','T','Y']:
                new_symbols = new_symbols + "P"
            else:
                new_symbols = new_symbols + "-"

        protein.add_track('test_track', symbols=new_symbols)

    track_list = []
    for protein in P:
        track_list.append(protein.track('test_track'))

    # write all the tracks out
    interfaces.si_tracks.write_tracks_from_list(track_list, 'output_test/all_values_test_from_list.tsv')

    # read the tracks into a new proteome
    P2 = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_tracks.add_tracks_from_file(P2, 'output_test/all_values_test_from_list.tsv', mode='symbols')

    # check everything worked ok
    for protein in P2:
        print(protein)
        for i in range(1,len(protein)+1):

            if protein.residue(i) in ['S','T','Y']:
                assert protein.track('test_track').symbol(i) == 'P'
            else:
                assert protein.track('test_track').symbol(i) == '-'



def test_write_symbols_tracks_from_list():

    test_data_dir = shephard.get_data('test_data')
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    
    P = uniprot.uniprot_fasta_to_proteome(fasta_file)

    # build phosphostring

    r1 = random.randint(1,20)
    r2 = random.random()
    for protein in P:



        # add symbol
        if random.random() < 0.5:            
            new_symbols =''
            for i in protein.sequence:

                if i in ['S','T','Y']:
                    new_symbols = new_symbols + "P"
                else:
                    new_symbols = new_symbols + "-"

            protein.add_track('test_track_1', symbols=new_symbols)

        # add value
        else:
            new_vals = []
            for i in protein.sequence:

                if i in ['S','T','Y']:
                    new_vals.append(r1)
                else:
                    new_vals.append(r2)

            protein.add_track('test_track_2', values=new_vals)
    
    track_list = []
    for protein in P:
        track_list.append(protein.track('test_track_1', safe=False))
        track_list.append(protein.track('test_track_2', safe=False))

    with pytest.raises(InterfaceException):
        interfaces.si_tracks.write_tracks_from_list(track_list, 'output_test/should_never_be_written.tsv')


def test_write_tracks_sanity_checks():

    with pytest.raises(InterfaceException):
        interfaces.si_tracks.write_tracks_from_list([1,None,2], 'output_test/should_never_be_written.tsv')

    
