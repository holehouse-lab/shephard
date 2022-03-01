import pytest

import shephard
from shephard.apis import uniprot  
from shephard import interfaces
import numpy as np
from shephard.exceptions import ProteinException 
import os



def test_remove_track():

    ##
    ## This code block just creates a proteome where we have 2 different
    ## typs of tracks
    test_data_dir = shephard.get_data('test_data')
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    pscore_track_file = '%s/%s' % (test_data_dir, 'TS1_tracks_pscore.tsv')

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)


    assert len(P.unique_track_names) == 0
    interfaces.si_tracks.add_tracks_from_file(P, pscore_track_file, mode='values')

    assert len(P.unique_track_names) == 1
    
    for protein in P:
        t = protein.track('pscore')
        protein.remove_track(t)
        break

    assert len(P.unique_track_names) == 1

    for protein in P:
        t = protein.track('pscore', safe=False)
        if t is not None:
            protein.remove_track(t, safe=False)
        
    assert len(P.unique_track_names) == 0


def test_remove_site():
    test_data_dir = shephard.get_data('test_data')
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    pscore_track_file = '%s/%s' % (test_data_dir, 'TS1_tracks_pscore.tsv')

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)

    assert len(P.unique_site_types) == 0

    interfaces.si_sites.add_sites_from_file(P, '%s/%s' % (test_data_dir, 'ts1_bonus_sites.tsv'))

    assert len(P.unique_site_types) == 4

    for protein in P:
        for site in protein.sites:
            protein.remove_site(site)

    assert len(P.unique_site_types) == 0
    interfaces.si_sites.add_sites_from_file(P, '%s/%s' % (test_data_dir, 'ts1_bonus_sites.tsv'))
    assert len(P.unique_site_types) == 4


    for protein in P:
        for site in protein.sites:
            if site.site_type == 'Sumoylation':
                protein.remove_site(site)

                # this should raise an exception
                with pytest.raises(ProteinException):
                    protein.remove_site(site)
                protein.remove_site(site, safe=False)


                # this should raise an exception
                with pytest.raises(ProteinException):
                    protein.remove_site('THIS IS NOT VALID')
                protein.remove_site('STILL NOT VALID', safe=False)
                


    assert len(P.unique_site_types) == 3            

    for protein in P:
        for site in protein.sites:
            if site.site_type == 'Sumoylation':
                protein.remove_site(site)


def test_remove_domains():
    test_data_dir = shephard.get_data('test_data')
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    domain_file = '%s/%s' % (test_data_dir, 'TS1_domains_idr.tsv')

    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    assert len(P.unique_domain_types) == 0

    interfaces.si_domains.add_domains_from_file(P, domain_file)
    assert len(P.unique_domain_types) == 1


    domain_file = '%s/%s' % (test_data_dir, 'TS1_domains_idr_and_others.tsv')
    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_domains.add_domains_from_file(P, domain_file)
    assert len(P.unique_domain_types) == 2

    
    protein = P.protein('O00401')
    protein.add_domain(1,40,'TEST')
    assert len(P.unique_domain_types) == 3
    domain_obj = protein.domains[1]
    protein.remove_domain(domain_obj)
    assert len(P.unique_domain_types) == 2
