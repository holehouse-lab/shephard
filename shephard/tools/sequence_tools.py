"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""
import re

## ------------------------------------------------------------------------
##
def build_mega_string(object_list, return_as_list=False): 
    
    """
    This takes a list of protein or domain SHEPHARD objects 
    and builds a single long str object that concatinates all
    object.sequence elements together (i.e. a "megastring").

    Allowed types for the object_list are Protein and Domain
    objects.

    This string can be used for simple statistical analysis
    of composition.

        
    Parameters
    ----------
    object_list : list 
        List of SHEPHARD objects with object.sequence variable,
        for example, a list of Domains or a list of Proteins

    return_as_list : bool
        If provided, rather than a single megastring, the 
        function returns a list of sequences from the objects
        in question.
        
    Returns
    ----------
    str or list
        Returns either a concatinated str object of the amino
        acid sequences associated with the passed object
    """


    megastring = ''
    for obj in object_list:  
        megastring = megastring + obj.sequence

    return megastring

## ------------------------------------------------------------------------
##
def find_string_positions(A,B):
    """
    Returns list of start positions where stringA is in stringB - including overlaps 
    """
    return [s.start() for s in re.finditer('(?=%s)' % (A), B)]
