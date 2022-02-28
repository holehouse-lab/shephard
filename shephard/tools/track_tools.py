"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (alex.holehouse@wustl.edu, g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""


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

