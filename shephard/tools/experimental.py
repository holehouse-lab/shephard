"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

import random

from shephard.tools import site_tools
from shephard.exceptions import UtilitiesException


def get_site_density_in_domain_normalized_by_protein(domain, site_types, sample_size=20, max_enrichment=100, min_enrichment=0.01):
    """
    Function that computes how enriched a domain is for one or more site
    types, relative to randomly-chosen regions of the same length taken
    from the same protein.

    The enrichment is the density of sites in the domain divided by the
    mean density of sites across ``sample_size`` random regions. Note
    that because the background is sampled randomly this function is
    stochastic - increase ``sample_size`` for a more stable estimate.

    Parameters
    ---------------
    domain : shephard.domain.Domain
        Domain object in which site density is calculated.

    site_types : str or list of str
        One or more site types to be counted.

    sample_size : int (default = 20)
        Number of random same-length regions used to build the background
        expectation.

    max_enrichment : float (default = 100)
        Value returned if the domain contains sites but the background
        sampling found none (i.e. enrichment would otherwise be infinite).

    min_enrichment : float (default = 0.01)
        Value returned if the domain contains no sites but the background
        sampling did find some (i.e. maximal depletion).

    Returns
    ---------------
    float
        Returns the fold enrichment of sites in the domain relative to the
        randomly-sampled background. A value of 1 means the domain matches
        the background expectation.

    """

    domain_size = len(domain)
    maxlen = len(domain.protein)

    # first calculate the density of sites in the domain. Note return_list=True is
    # important here - without it we would count the number of positions that have
    # one or more sites, rather than the number of sites
    density_of_sites_in_domain = len(domain.get_sites_by_type(site_types, return_list=True))/domain_size

    # next, for sample_size iterations we're going to...
    count = 0
    for _ in range(0, sample_size):

        # ... randomly select another region of the protein that is the same size as
        # the domain. Note we're working in real-world (1-indexed, inclusive) positions
        # here, hence the +1s
        d_start = random.randint(1, maxlen - domain_size + 1)
        d_end   = d_start + domain_size - 1

        # and then count how many sites fall in that region
        count = count + len(domain.protein.get_sites_by_type_and_range(site_types, d_start, d_end, return_list=True))

    # compute on average how many sites one finds and then use that number to compute
    # the expected density across a domain-sized region
    mean_expected = (count/sample_size)/domain_size

    # if the domain has no sites then we're either uninformative (no sites anywhere)
    # or maximally depleted (sites in the background but not here)
    if density_of_sites_in_domain == 0:
        if count == 0:
            return 1
        else:
            return min_enrichment

    # domain has sites but the background found none, so enrichment is unbounded
    if count == 0:
        return max_enrichment

    return density_of_sites_in_domain/mean_expected


def get_site_density_percentile_normalized_by_protein(domain, site_types):
    """
    Function that reports where a domain's site density sits in the
    distribution of site densities calculated across all same-length
    windows in the parent protein.

    Parameters
    ---------------
    domain : shephard.domain.Domain
        Domain object in which site density is calculated.

    site_types : str or list of str
        One or more site types to be counted.

    Returns
    ---------------
    float
        Returns a value between 0 and 1 which reports the fraction of
        same-length windows in the protein whose site density is less
        than or equal to that of the domain. A value of 1 means the
        domain is the most site-dense region in the protein.

    """

    domain_length = len(domain)

    # get density of sites in the domain. Note we deliberately count POSITIONS that
    # have one or more sites here (i.e. the length of the returned dictionary) and
    # not the number of sites, because build_site_density_vector() below treats site
    # presence as a binary per-residue phenomenon. The two numbers must be counted
    # the same way for the comparison to mean anything
    density_of_sites_in_domain = len(domain.get_sites_by_type(site_types))/domain_length

    # build the distribution of site densities over every same-length window
    density_vector  = site_tools.build_site_density_vector(domain.protein, site_types=site_types, window_size=domain_length, append_leading_lagging=False)

    # this should not happen for a valid domain, but avoids a divide-by-zero if it does
    if len(density_vector) == 0:
        raise UtilitiesException(f'Could not build a site density vector for the protein associated with domain {domain}')

    density_vector.sort()

    # note we count with <= rather than testing for equality; the domain density is a
    # float built from a different calculation, so exact equality is not guaranteed
    c = 0
    for i in density_vector:
        if i <= density_of_sites_in_domain:
            c = c + 1
        else:
            break

    return c/len(density_vector)




    
                
            
