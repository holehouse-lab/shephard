import pytest

import shephard
from shephard.apis import uniprot  
from shephard import interfaces
import numpy as np
from shephard.exceptions import ProteinException 





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


    

    


