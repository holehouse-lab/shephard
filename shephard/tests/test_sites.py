def test_unique_domain_types(TS1_domains2_sites_tracks_protein_attributes):
    # check we corretly count the number of unique domains

    # deep cut that relies on a bunch of things working
    assert TS1_domains2_sites_tracks_protein_attributes.protein('O00401').site(242)[1].attribute('ENZYME',safe=False) == 'TNK21'

    assert TS1_domains2_sites_tracks_protein_attributes.protein('O00401').site(307)[1].symbol == 'R'

    print(TS1_domains2_sites_tracks_protein_attributes.protein('O00401').sequence[1])
    assert TS1_domains2_sites_tracks_protein_attributes.protein('O00401').site(307)[1].symbol == TS1_domains2_sites_tracks_protein_attributes.protein('O00401').site(307)[1].residue

    

