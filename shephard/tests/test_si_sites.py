import pytest

from shephard.exceptions import ProteinException
from shephard.exceptions import InterfaceException
import shephard
from shephard.interfaces import si_sites
from shephard.apis import uniprot  


test_data_dir = shephard.get_data('test_data')


def test_si_site_add_from_file():

    TS1 = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir,'testset_1.fasta'))

    si_sites.add_sites_from_file(TS1, '%s/%s' % (test_data_dir, 'ts1_bonus_sites.tsv'))

    O00470 = TS1.protein('O00470')

    # check length
    assert len(O00470) == 390
    assert len(O00470.site(1)) == 1
    assert len(O00470.site(390)) == 1

    # get site a position 2 that doesn't exist
    with pytest.raises(ProteinException) as e_info:
        O00470.site(2)

    # get site at position 2 that doesn't exist
    assert O00470.site(2, safe=False) == []


def test_si_site_add_from_file_test_robustness():

    TS1 = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir,'testset_1.fasta'))

    # this should fail - skip_bad influences skipping
    with pytest.raises(ProteinException) as e_info:
        si_sites.add_sites_from_file(TS1, '%s/%s' % (test_data_dir, 'ts1_bonus_sites_bad.tsv'), skip_bad=True)

    with pytest.raises(ProteinException) as e_info:
        si_sites.add_sites_from_file(TS1, '%s/%s' % (test_data_dir, 'ts1_bonus_sites_bad.tsv'), skip_bad=False)

    # but if we load in with safe=False will skip over the bad site
    si_sites.add_sites_from_file(TS1, '%s/%s' % (test_data_dir, 'ts1_bonus_sites_bad.tsv'), safe=False)
    O00470 = TS1.protein('O00470')

    assert len(O00470) == 390



        
def test_si_site_add_file_read_robustness():

    TS1 = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir,'testset_1.fasta'))

    # thsi should fail - malformatted line
    with pytest.raises(InterfaceException) as e_info:    
        si_sites.add_sites_from_file(TS1, '%s/%s' % (test_data_dir, 'ts1_bonus_sites_bad2.tsv'), skip_bad=False)


    si_sites.add_sites_from_file(TS1, '%s/%s' % (test_data_dir, 'ts1_bonus_sites_bad2.tsv'), skip_bad=True)
    assert len(TS1.sites) == 4


def test_write_sites():

    
    TS1 = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir,'testset_1.fasta'))

    for p in TS1:
        
        for s in range(1,len(p)+1):
            if p.residue(s) in ['S','P','E']:
                p.add_site(s, 'test_site')

    si_sites.write_sites(TS1, 'output_test/test_sites.tsv')


    TS2 = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir,'testset_1.fasta'))
    si_sites.add_sites_from_file(TS2, 'output_test/test_sites.tsv')

    for p in TS2:
        
        for s in range(1,len(p)+1):
            if p.residue(s) in ['S','P','E']:
                assert len(p.site(s)) == 1



def test_write_sites_with_attributes():

    
    TS1 = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir,'testset_1.fasta'))

    for p in TS1:
        
        for s in range(1,len(p)+1):
            if p.residue(s) in ['S','P','E']:
                p.add_site(s, 'test_site', symbol='A', value=0.2, attributes={'test_at_1':0, 'test_at_2':'yeah', 'test_at_3':1/33})

    si_sites.write_sites(TS1, 'output_test/test_sites.tsv')


    TS2 = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir,'testset_1.fasta'))
    si_sites.add_sites_from_file(TS2, 'output_test/test_sites.tsv')

    for p in TS2:
        
        for s in range(1,len(p)+1):
            if p.residue(s) in ['S','P','E']:
                assert len(p.site(s)) == 1
                assert p.site(s)[0].attribute('test_at_1') == '0'
                assert p.site(s)[0].attribute('test_at_2') == 'yeah'
                assert p.site(s)[0].attribute('test_at_3') == str(1/33)
                assert p.site(s)[0].symbol == 'A'
                assert p.site(s)[0].value == 0.2


def test_write_sites_from_list():

    
    TS1 = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir,'testset_1.fasta'))

    for p in TS1:
        
        for s in range(1,len(p)+1):
            if p.residue(s) in ['M','A','K']:
                p.add_site(s, 'test_site', symbol='A', value=0.2, attributes={'test_at_1':0, 'test_at_2':'yeah', 'test_at_3':1/33})

    all_sites = TS1.sites
    
    si_sites.write_sites_from_list(all_sites, 'output_test/test_sites.tsv')


    TS2 = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir,'testset_1.fasta'))
    si_sites.add_sites_from_file(TS2, 'output_test/test_sites.tsv')

    for p in TS2:
        
        for s in range(1,len(p)+1):
            if p.residue(s) in ['M','A','K']:
                assert len(p.site(s)) == 1
                assert p.site(s)[0].attribute('test_at_1') == '0'
                assert p.site(s)[0].attribute('test_at_2') == 'yeah'
                assert p.site(s)[0].attribute('test_at_3') == str(1/33)
                assert p.site(s)[0].symbol == 'A'
                assert p.site(s)[0].value == 0.2


            
            
