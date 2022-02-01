"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""


def test_site_track_functions(TS1_domains2_sites_tracks_protein_attributes):
    # check we corretly count the number of unique domains

    # deep cut that relies on a bunch of things working
    assert TS1_domains2_sites_tracks_protein_attributes.protein('O00401').site(242)[1].attribute('ENZYME',safe=False) == 'TNK21'

    assert TS1_domains2_sites_tracks_protein_attributes.protein('O00401').site(307)[1].symbol == 'R'

    print(TS1_domains2_sites_tracks_protein_attributes.protein('O00401').sequence[1])
    assert TS1_domains2_sites_tracks_protein_attributes.protein('O00401').site(307)[1].symbol == TS1_domains2_sites_tracks_protein_attributes.protein('O00401').site(307)[1].residue

    


def test_site_domain_functions(TS1_domains2_sites_tracks_protein_attributes):

    c= 0
    for site in TS1_domains2_sites_tracks_protein_attributes.sites:

        d = site.get_domains()
        
        for domain in d:
            assert domain.start <= site.position
            assert domain.end >= site.position
            c= c+1


    # 126 is expected number for this file
    assert c == 126
                
        
