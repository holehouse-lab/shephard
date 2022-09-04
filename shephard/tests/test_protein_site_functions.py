"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

import pytest
import shephard
from shephard.interfaces import si_sites
from shephard.exceptions import ProteinException
from shephard.apis import uniprot  

test_data_dir = shephard.get_data('test_data')


def test_get_site_by_position():
    TS1 = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir,'testset_1.fasta'))

    si_sites.add_sites_from_file(TS1, '%s/%s' % (test_data_dir, 'ts1_bonus_sites.tsv'))

    pos = TS1.protein('O00470').get_sites_by_position(72)    
    assert len(pos[72]) == 1
    assert pos[72][0].site_type == 'Sumoylation'


    pos = TS1.protein('O00470').get_sites_by_position(72, wiggle=5)    
    assert len(pos[72]) == 1
    assert pos[72][0].site_type == 'Sumoylation'


    pos = TS1.protein('O00470').get_sites_by_position(72, wiggle=5, return_list=True)    
    assert len(pos) == 1
    assert pos[0].site_type == 'Sumoylation'



    pos = TS1.protein('O00470').get_sites_by_position(196)    
    assert len(pos[196]) == 1
    assert pos[196][0].site_type == 'Phosphoserine'


    pos = TS1.protein('O00470').get_sites_by_position(197, wiggle=5)    
    assert len(pos) == 2
    assert pos[196][0].site_type == 'Phosphoserine'


    pos = TS1.protein('O00470').get_sites_by_position(197, wiggle=5, return_list=True)    
    assert len(pos) == 2
    assert pos[0].site_type == 'Phosphoserine'
    assert pos[1].site_type == 'Phosphoserine'


    pos = TS1.protein('O00470').get_sites_by_position(100, wiggle=5, return_list=True)    
    assert len(pos) == 0



def test_get_site_by_range():
    TS1 = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir,'testset_1.fasta'))

    si_sites.add_sites_from_file(TS1, '%s/%s' % (test_data_dir, 'ts1_bonus_sites.tsv'))

    pos = TS1.protein('O00470').get_sites_by_range(55,79)    
    assert len(pos[72]) == 1
    assert pos[72][0].site_type == 'Sumoylation'


    pos = TS1.protein('O00470').get_sites_by_range(70,80, wiggle=5)    
    assert len(pos[72]) == 1
    assert pos[72][0].site_type == 'Sumoylation'


    pos = TS1.protein('O00470').get_sites_by_range(55,79, wiggle=5, return_list=True)    
    assert len(pos) == 1
    assert pos[0].site_type == 'Sumoylation'



    pos = TS1.protein('O00470').get_sites_by_range(150,250)    
    assert len(pos[196]) == 1
    assert pos[196][0].site_type == 'Phosphoserine'


    pos = TS1.protein('O00470').get_sites_by_range(150,250, wiggle=5)    
    assert len(pos) == 2
    assert pos[196][0].site_type == 'Phosphoserine'


    pos = TS1.protein('O00470').get_sites_by_range(150,250, wiggle=5, return_list=True)    
    assert len(pos) == 2
    assert pos[0].site_type == 'Phosphoserine'
    assert pos[1].site_type == 'Phosphoserine'

    pos = TS1.protein('O00470').get_sites_by_range(30,60, wiggle=5, return_list=True)    
    assert len(pos) == 0


def test_get_site_by_type():
    TS1 = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir,'testset_1.fasta'))

    si_sites.add_sites_from_file(TS1, '%s/%s' % (test_data_dir, 'ts1_bonus_sites.tsv'))

    
    pos = TS1.protein('O00470').get_sites_by_type('Phosphoserine')
    assert len(pos[196]) == 1
    assert pos[196][0].site_type == 'Phosphoserine'
    assert pos[198][0].site_type == 'Phosphoserine'

    pos = TS1.protein('O00470').get_sites_by_type('Phosphoserine', return_list=True)
    assert len(pos) == 2


    pos = TS1.protein('O00470').get_sites_by_type('PhosphoX')
    assert len(pos) == 0

    pos = TS1.protein('O00470').get_sites_by_type('PhosphoX', return_list=True)
    assert len(pos) == 0

    
    # add a site, check we can find it, remove it and check it's gone!
    TS1.protein('O00470').add_site(1, 'PhosphoX')
    
    pos = TS1.protein('O00470').get_sites_by_type('PhosphoX', return_list=True)
    assert len(pos) == 1

    TS1.protein('O00470').remove_site(pos[0])
    pos = TS1.protein('O00470').get_sites_by_type('PhosphoX', return_list=True)
    assert len(pos) == 0


def test_get_site_by_type_and_range():
    TS1 = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir,'testset_1.fasta'))

    si_sites.add_sites_from_file(TS1, '%s/%s' % (test_data_dir, 'ts1_bonus_sites.tsv'))

    
    pos = TS1.protein('O00470').get_sites_by_type_and_range('Phosphoserine', 190, 200)
    assert len(pos[196]) == 1
    assert pos[196][0].site_type == 'Phosphoserine'
    assert pos[198][0].site_type == 'Phosphoserine'

    pos = TS1.protein('O00470').get_sites_by_type_and_range('Phosphoserine', 190, 200, return_list=True)
    assert len(pos) == 2


    pos = TS1.protein('O00470').get_sites_by_type_and_range('PhosphoX', 1, 50)
    assert len(pos) == 0

    pos = TS1.protein('O00470').get_sites_by_type_and_range('PhosphoX', 1, 50, return_list=True)
    assert len(pos) == 0

    
    # add a site, check we can find it, remove it and check it's gone!
    TS1.protein('O00470').add_site(1, 'PhosphoX')
    
    pos = TS1.protein('O00470').get_sites_by_type_and_range('PhosphoX', 1, 50, return_list=True)
    assert len(pos) == 1

    TS1.protein('O00470').remove_site(pos[0])
    pos = TS1.protein('O00470').get_sites_by_type_and_range('PhosphoX', 1, 50, return_list=True)
    assert len(pos) == 0




    

    
    
