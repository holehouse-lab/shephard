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
def find_string_positions(query, target, protein_indexing=True):
    """
    Returns list of start positions where 
    stringA is in stringB - including overlaps.

    Note that by default the indices use 1-indexing so that this works
    directly with protein sequence numbering. However, for manipulating
    Python strings this may be undesirable and 0 indexing may be better,
    in which case setting protein_indexing = False will address this.

    Practically, this uses the re regex expression under the hood 
    and searches left-to-right across the target, so if you want 
    to get fancier with your searching you can always pass in a 
    regular expression.

    Examples
    ----------------

    Conveninet regular expression syntax includes:

        1. ``'.'`` for wildcards (e.g. ``'L.P'`` would match an L and P around any other character

        2. ``[A|C]`` for requiring matching of a subset of residues (e.g. residue A and C).

    But the python re module has a fairly complex pattern matching
    ability

    Parameters
    --------------
    query : str
        The search query. 

    target : str
        The string that we'll search for 1 or more entries
        of the query

    protein_indexing : bool
        Flag which, if set to True, means the first residue in 
        a string indexes at '1' instead of '0' (as would be normal
        in Python. If set to False, then indexing is done from 0.
     

    Returns
    ----------
    list
         Returns a list with the start positions
    
    """
    if protein_indexing is True:
        offset = 1
    else:
        offset = 0

    return [s.start() + offset for s in re.finditer('(?=%s)' % (query), target)]
