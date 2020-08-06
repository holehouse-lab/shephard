"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

def get_site_density_in_domain_normalized_by_protein(domain, site_types, sample_size=20, max_enrichment=100, min_enrichment=0.01):
    """

    """
    
    # first calculate the density of sites in the domain
    density_of_sites_in_domain = len(domain.get_sites_by_type(site_types))/len(domain)
    
    subsample=[]

    # next, for sample_size iteractions we're going to...
    count=0
    for i in range(0, sample_size):

        # next we're going to randomly select another region of the protein that
        # is the same size as the domain...
        maxlen = len(domain.protein)
        domain_size = len(domain)    
        d_start = random.randint(0,maxlen-domain_size)        
        d_end   = d_start + domain_size

        # and then count how many hits in the site
        count = count + len(domain.protein.get_sites_by_type_and_range(site_types, d_start, d_end))


    # we do this to avoid dividing by zero - i.e. compute on average how many sites one finds
    # and then use that number to compute the density across the domain)
    mean_expected = (count/sample_size)/len(domain)

    # if we didn't find sites anywhere
    if density_of_sites_in_domain == 0:
        if count == 0:
            return 1
        else:
            return max_enrichment
    if count == 0:
        return min_enrichment

    return density_of_sites_in_domain/mean_expected


def get_site_density_percentile_normalized_by_protein(domain, site_types):

    domain_length = len(domain)

    # get density of sites in the domain
    density_of_sites_in_domain = len(domain.get_sites_by_type(site_types))/domain_length

    # 
    density_vector  = site_tools.build_site_density_vector(domain.protein, site_types=site_types, window_size=domain_length, append_leading_lagging=False)
    density_vector.sort()

    c=0
    for i in density_vector:
        c=c+1
        if i == density_of_sites_in_domain:
            return c/len(density_vector)




    
                
            
