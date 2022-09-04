"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

from shephard import general_utilities
import numpy as np

def build_site_density_vector(protein, site_types=None, window_size=30, append_leading_lagging=True):
    """
    Function that constructs a sliding-window density vector of sites along
    a protein.

    site_types is a list of one or site types.

    This tool is stateless - i.e. it does not alter the passed protein
    but instead only generates a numerical list which could be added as 
    a track.

    Parameters
    ----------------

    protein : shephard.protein.Protein object
        Protein object over which sites are identified

    site_type : str or list of strings
        One or more possible site_types that may be found in the protein. 
        Either a single string or a list of strings can be passed, allowing 
        for one or more sites to be grouped together

    window_size : int
        Size of sliding window over which site density is calculated

    append_leading_lagging : Bool
        Flag that if true will mean the function returns a numerical vector
        equal in length of the protein. If false, will return a shorter 
        vector and not add leading/lagging values.

    Returns
    -------------
    list

        Returns a list of values equal to the length of the protein, where
        the value at each position reports on the local denisty of sites 
        averaged over the window_size.

    """

    # give ability to take a single string or a list of strings if we wish
    # to compare against multype site types
    if site_types is not None:
        site_types = general_utilities.string_to_list_of_strings(site_types)

    nres = len(protein)

    # first build an empty vector of 0s equal to the length of the sequence    
    all_res = [0]*len(protein)

    # then for each site assign a '1' to positions where a site exists
    for position in protein.site_positions:
        sites = protein.site(position)

        # if we didn't specifit a site type assume all sites are fair game
        if site_types is None:
            all_res[position - 1] = 1

        # else validate against the set of available sites and IF the site
        # matches the requested
        else:

            # for every site object found at $position
            for s in sites:
                if s.site_type in site_types:
                    all_res[position - 1] = 1
                            
    # finally we're going to calculate the density of sites from this vector
    # note we're treating site presence as a binary phenomenon - ie a residue
    # has a site or does not
    density_vector = []
    for pos in range(0, (nres - window_size)+1):

        local_density = np.sum(all_res[pos:pos+window_size])/window_size        

        density_vector.append(local_density)

    # having built a density vector that is nres-window_size+1 in length, we now need to extend the N and C
    # termini such that the len(denisty_vector) = len(protein)

    # this code creates leading/lagging values to fill in the missing ones such that the actual value
    # of the density vector reports on the density half-way across the window_size and the 
    # track length = nres
    leading_values = [density_vector[0]]*int(window_size/2)
    lagging_values = [density_vector[-1]]*(nres - (len(density_vector) + len(leading_values)))
          
    # this line then combines the leading, desnity and lagging lists into a single list, and we
    # then add this numerical list as a track
    final = leading_values + density_vector + lagging_values    

    return final
     

    
