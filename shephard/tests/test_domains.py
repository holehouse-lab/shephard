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

names_list = ['O00401', 'O00470', 'O00472', 'O00499', 'O00629', 'O00712', 'O00716', 'O14786', 'Q9UJX3']

def test_bad_domain_add(TS1_domains2_sites_tracks):    
    """
    Note if any of these succeed the test suite will fail because it incrases
    """

    X = TS1_domains2_sites_tracks

    # add a domain where end is outside
    with pytest.raises(ProteinException):
        assert X.protein('O00401').add_domain(1,1000, 'fail')


    # add a domain where end is outside
    with pytest.raises(ProteinException):
        assert X.protein('O00401').add_domain(-1,20, 'fail')


def test_domain_add(TS1_domains2_sites_tracks):    
    """
    This function tests our ability to correctly add and then
    make reference to a domain
    
    """


    local_protein = TS1_domains2_sites_tracks.protein('O00629')
    
    domain_type = 'test_domain'
    local_protein.add_domain(1,10, domain_type)
    local_protein.add_domain(512,521, domain_type)

    
    # first check this works, if not, we all got done some real
    # bad codes in these parts
    assert local_protein.sequence[0:10] == 'MADNEKLDNQ'


    # ......................................................................
    # First try the N-terminal domain    
    #
    # this is just a sanity check because if THIS fails the
    # rest of the test definietly should not work

    local_domain = local_protein.get_domains_by_type(domain_type)[0]
    assert local_domain.sequence == 'MADNEKLDNQ'
    assert len(local_domain) == 10
    assert local_protein.get_sequence_region(1,10) == local_domain.sequence
    assert local_domain.start == 1
    assert local_domain.end == 10

    # these are more code coverage
    assert local_domain.protein.unique_ID == local_protein.unique_ID
    assert local_domain.domain_type == domain_type

    # assumes our domain name construction doesn't change. 
    assert local_domain.domain_name == "%s_%i_%i" % (domain_type, local_domain.start, local_domain.end)


    # note this should (A) work and (B) discard the N-terminal context
    assert local_protein.get_sequence_context(1,9) == local_domain.sequence

    # ......................................................................
    # Next check out C-terminal dmain
    #
    # next check out the
    local_domain = local_protein.get_domains_by_type(domain_type)[1]
    assert local_domain.sequence == 'ANVPTEGFQF'
    assert len(local_domain) == 10
    assert local_protein.get_sequence_region(512, 521) == local_domain.sequence
    assert local_domain.start == 512
    assert local_domain.end == 521

    # assumes our domain name construction doesn't change. 
    assert local_domain.domain_name == "%s_%i_%i" % (domain_type, local_domain.start, local_domain.end)

        
    # note this should (A) work and (B) discard the C-terminal context
    assert local_protein.get_sequence_context(521,9) == local_domain.sequence


def test_domain_attributes(TS1_domains2_sites_tracks):    
    
    # set up
    local_protein = TS1_domains2_sites_tracks.protein('O00629')
    domain_type = 'test_domain'
    local_protein.add_domain(20,200, domain_type, attributes={'test_1':1,'test_2':'yes' })
    local_domain = local_protein.get_domains_by_type(domain_type)[0]
    local_domain.add_attribute('added_attribute', np.array([1,2,3,4]))

    assert local_domain.attribute('test_1') == 1
    assert local_domain.attribute('test_2') == 'yes'
    assert local_domain.attribute('added_attribute').all() == np.array([1,2,3,4]).all()
    
    
    # no attribute called test_missing so should gracefully fail
    with pytest.raises(DomainException):
        assert local_domain.attribute('test_missing') == 2


    # no attribute called test_missing so should gracefully fail
    with pytest.raises(DomainException):
        local_domain.add_attribute('test_1', 2)

    # but we should be able to force override an attribute if safe=False
    local_domain.add_attribute('test_1', 2, safe=False)
    assert local_domain.attribute('test_1') == 2
        

    
def test_inside_domain(TS1_domains2_sites_tracks):
    """
    This function tests our ability to correctly identify if a 
    position of site is in a domain 
    
    """
    
    local_protein = TS1_domains2_sites_tracks.protein('O00629')
    
    # get first domain 1-63 
    local_domain = local_protein.domains[0]
    
    # get first site 2
    local_site_position = local_protein.sites[0].position
    
    # this should pass
    assert local_domain.inside_domain(local_site_position) == True
    
    # this should also pass and it is False (!=) 
    assert local_domain.inside_domain(100) == False
    
    ## CURRENTLY NO CHECKS WRITTEN FOR THIS TEST 
    # passed value is outside the length of the protien
    # with pytest.raises(ProteinException):
    #    assert local_domain.inside_domain(4000) == False
    
    
    
def test_domain_overlap(TS1_domains2):
    
    """
    This function tests our ability to correctly identify if a 
    two domains overlap in a single protein 
    """

    test_data_dir = shephard.get_data('test_data')
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    domain_file = '%s/%s' % (test_data_dir, 'TS1_domains_idr.tsv')

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_domains.add_domains_from_file(P, domain_file)
    
    local_protein = P.protein('O00629')
        
    # get first domain 1-10
    local_protein.add_domain(1, 10, 'overlapping_domain')
    print(local_protein.domains)

    # sub domain 10
    local_domain = local_protein.domains[1]
    
    # domain 1 - 63
    local_domain_true = local_protein.domains[0]
    
    # domain 489-521
    local_domain_false = local_protein.domains[2]
    
    
    assert local_domain.domain_overlap(local_domain_true) == True
    assert local_domain.domain_overlap(local_domain_false) == False
    
    # check for passage of domains from two different protiens 
    with pytest.raises(DomainException):
        assert local_domain.domain_overlap(TS1_domains2.protein('O00716').domains[0])
        
        
def test_domain_track_functions(TS1_domains2_sites_tracks):
    
    """
    This function tests our ability to properly interface with 
    protien associated tracks from the domain object
    
    missing here is the check of a proper get track_symbols 
    """
    
    local_domain = TS1_domains2_sites_tracks.protein('O00629').domains[0]
    local_track = TS1_domains2_sites_tracks.protein('O00629').tracks[0]
    local_trackname = local_track.name
    
    nonassociated_track = TS1_domains2_sites_tracks.protein('O00716').tracks[0]
    nonassociated_trackname = 'FAILTRACK'
   
    # check that extracted track is correct
    assert local_domain.get_track_values(local_trackname) == local_track.values[0:63] 
    
    # check that error is raise when wrong track type is passed or track has no values 
    with pytest.raises(DomainException):
        assert local_domain.get_track_symbols(local_trackname)
    
    # check for missed passed track name 
    with pytest.raises(ProteinException):
        assert local_domain.get_track_values(nonassociated_trackname)
        
    # check for missed passed track name 
    with pytest.raises(ProteinException):
        assert local_domain.get_track_values(nonassociated_trackname)
        
        
def test_domain_site_functions(TS1_domains2_sites):
    
    """
    This function tests our ability to properly interface with 
    protien associated sites from the domain object
    """
    
    local_protein = TS1_domains2_sites.protein('O00629')
    local_domain = local_protein.domains[0]
    true_site_type = 'Phosphoserine' 
    false_site_type = 'Phosphotyrosine' # this site type is outside the domain 
    
    sites_in_region = [s for s in local_protein.sites if local_domain.start <= s.position <= local_domain.end]

    ## basic self contianed functions 
    assert local_domain.sites == [s for s in local_protein.sites if local_domain.start <= s.position <= local_domain.end]
    assert local_domain.site_positions == [s.position for s in sites_in_region]
    
    ## fuctions that depend on user input 
    # remember local_domain.site() returns a list cuz there can be multible sites at one posision 
    assert local_domain.site(2)[0] == local_protein.sites[0]
    
    # remember local_domain.get_sites_by_type() returns a dictionary cuz there can be multible sites locations
    sites_type_in_region = [s for s in sites_in_region if s.site_type == true_site_type]
    assert local_domain.get_sites_by_type(true_site_type) == {s.position:[s] for s in sites_type_in_region}
    
    sites_type_in_region = [s for s in sites_in_region if s.site_type == false_site_type]
    assert local_domain.get_sites_by_type(false_site_type) == {}
