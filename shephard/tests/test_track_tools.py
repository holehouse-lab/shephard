"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

import pytest
from shephard.exceptions import ProteinException, DomainException
import numpy as np
import shephard
from shephard.apis import uniprot
from shephard import interfaces
from shephard.tools import track_tools

names_list = ['O00401', 'O00470', 'O00472', 'O00499', 'O00629', 'O00712', 'O00716', 'O14786', 'Q9UJX3']


def test_domain_to_track_1(TS1_domains2_sites_tracks):    
    """
    Note if any of these succeed the test suite will fail because it incrases
    """



    # build tracks from domains
    track_dict = track_tools.build_track_from_domains(TS1_domains2_sites_tracks)

    # check they're the right length
    for k in track_dict:
        assert len(track_dict[k]) == len(TS1_domains2_sites_tracks.protein(k))

    # check positions inside the domain are in fact set to 1
    for k in track_dict:
        for domain in TS1_domains2_sites_tracks.protein(k).domains:
            assert track_dict[k][domain.start-1] == "1"
            assert track_dict[k][domain.end-1] == "1"
            assert track_dict[k][domain.end-len(domain)] == "1"



    


def test_domain_to_track_2(TS1_domains2_sites_tracks):    
    """

    """



    # build tracks from domains
    track_dict = track_tools.build_track_from_domains(TS1_domains2_sites_tracks, domain_type='pscore_domain')


    # check they're the right length
    for k in track_dict:
        assert len(track_dict[k]) == len(TS1_domains2_sites_tracks.protein(k))

    # check positions inside the domain are in fact set to 1
    for protein in TS1_domains2_sites_tracks:
        k = protein.unique_ID

        for domain in protein.domains:

            if domain.domain_type == 'pscore_domain':
                print(domain.protein.unique_ID)
                assert track_dict[k][domain.start-1] == "1"
                assert track_dict[k][domain.end-1] == "1"
                assert track_dict[k][domain.end-len(domain)] == "1"



    # has no pscore_domains but does have an IDR domain
    k = 'O00629'
    for domain in TS1_domains2_sites_tracks.protein(k).domains:
        if domain.domain_type == 'IDR':
            
            # should NOT be a 1 because these domains are not pscore domains
            assert track_dict[k][domain.start-1] == "0"
            assert track_dict[k][domain.end-1] == "0"
            assert track_dict[k][domain.end-len(domain)] == "0"




