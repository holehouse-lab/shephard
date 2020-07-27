"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

import pytest

names_list = ['O00401', 'O00470', 'O00472', 'O00499', 'O00629', 'O00712', 'O00716', 'O14786', 'Q9UJX3']

def test_load_protein_attribute(TS1_domains2_sites_tracks_protein_attributes):    

    assert TS1_domains2_sites_tracks_protein_attributes.protein('O00401').attribute('GO3') == 'tags from a separate'

    assert len(TS1_domains2_sites_tracks_protein_attributes.protein('O00401').attributes) == 5




