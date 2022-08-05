"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (alex.holehouse@wustl.edu, g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

import numpy as np
from shephard import general_utilities

## ------------------------------------------------------------------------
##
def binerize(input_vector, threshold, mode='above'):
    """
    Function that takes an input vector (which can be any iterable datatype) and
    will cyle over each element and assess if the value is above (or below) the 
    threshold. If above the value is set to 1. If below, the value is set to 0.
    

    Parameters
    --------------
    input_vector : list, numpy.ndarray (1D), other
        Any data type that can be iterated through and the content evaluated 
        using a greater than or less than symbol.

    threshold : float/int
        Threshold used on a per-element basis to define a 1 or a 0.

    mode : str  (default = 'above')
        Selector that determines if threshold is a min or max to define 
        0 vs. 1. Must be one of 'above' or 'below'.
        
    Returns
    ----------
    list
        Returns a list of 1s, or 0s that is the same length of the 
        input_vector.

    """
    
    return_list = []

    # check the passed mode keyword
    general_utilities.valid_keyword('mode', mode, ['above', 'below'])

    if mode == 'above':
        for i in input_vector:
            if i > threshold:
                return_list.append(1)
            else:
                return_list.append(0)
    else:
        for i in input_vector:
            if i < threshold:
                return_list.append(1)
            else:
                return_list.append(0)
                        
    return return_list



## ------------------------------------------------------------------------
##
def build_track_from_domains(proteome, domain_type=None):

    """
    Function that builds a dictionary of unique IDs to binary symbol tracks
    that map the domains in the protein 1s or 0s.

    Note this function returns a dictionary which can then be used to update 
    the underlying Proteome object, but does not itself alter the Proteome
    object.


    Parameters
    ----------------
    proteome : Proteome
        The Proteome which is going to be scanned for each track. Note that
        the underlying Proteome is not altered by this function.

    domain_type : str
        The domain type to be used. If none provided then all domains are
        used. Default = None.


    Returns
    ---------------
    dict
        Returns a dictionary where keys are unique IDs for every protein in 
        the passed proteome and the values is a list of str where positions are
        1 (if in a domain) or 0 (if not in a domain).


    """
    
    tracks = {}

    # for each protein  in the proteome
    for protein in proteome:

        # construct an empty track
        raw = np.zeros((len(protein)), dtype=int)

        # for each domain in the protein
        for d in protein.domains:


            if domain_type is None:

                # assign the domain positions to 1
                raw[d.start-1:d.end] = 1

            elif domain_type == d.domain_type:
                raw[d.start-1:d.end] = 1

        # finally convert to a list of strings to create a symbols track
        tracks[protein.unique_ID] = [str(i) for i in raw]

        
    return tracks
                             
    
