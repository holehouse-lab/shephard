from shephard import general_utilities

def binerize(input_vector, threshold, mode='above'):
    """
    Function that takes an input vector (which can be any iterable datatype) and
    will cyle over each element and assess if the value is above (or below) the 
    threshold. If above the value is set to 1. If below, the value is set to 0.
    

    Parameters
    --------------
    input_vector : list, numpy.ndarray (1D), other
        Any data type that can be iterated through and the content evaluated using 
        a greater than or less than symbol.

    threshold : float/int
        Threshold used on a per-element basis to 

    mode : str ('above' or 'below')
        Selector that determines if threshold is a min or max to define 0 vs. 1

    Returns
    ----------
    list
        Returns a list of 1s, or 0s that is the same length of the input_vector

    """
    
    return_list = []

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

