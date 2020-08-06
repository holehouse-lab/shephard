"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

from shephard.exceptions import InterfaceException


## ------------------------------------------------------------------------
##
def check_proteome(p, function_name):
    """
    Function that takes takes in some object and tests if its a Proteome object (or not).
    If yes returns None, but if no raise an exception with an error message that includes the
    name of the parent function that is calling this validation function

    Parameters
    --------------
    p : Proteome object 
        Unknown object that will be tested to see if it's a Proteome object or not

    function_name : string
        Passed string that makes it easy for the user to debug which function has failed

    Returns
    ---------
    No returns but the function will raise an exception if the passed object is not a Proteome
    object

    """
    if "<class 'shephard.proteome.Proteome'>" == str(p.__class__):
        return None
    else:
        raise InterfaceException('First argument passed to function [%s] was not a proteome' %(function_name))



## ------------------------------------------------------------------------
##
def parse_key_value_pairs(split_line, filename, linecount, line):
    """
    Helper function for parsing input files that have attributes (for Domains and Sites).

    The function takes in a list of key-value pairs (where key-values are split by a ':'
    symbol) and returns a parsed dictionary. Note that values will always be strings.

    Parameters
    -------------
    split_line : list of strings
        Each string in the list should have the format <KEY> : <VALUE> 

    filename : string
        Name of the file the calling function is parsing. Only used when raising an exception.

    linecount : int
        Current line number the file processing is on. Only used when raising an exception.

    line : string
        Full line that the file processing is on. Again, only used when raising an exception

    """

    attributes = {}

    for idx in range(len(split_line)):
        try:
            sentry = split_line[idx].split(':')
        except Exception:
                
            # should update this to also display the actual error...
            raise InterfaceException('Failed parsing key-value pairs in file [%s] on line [%i]... line printed below:\n%s'%(filename, linecount, line))
            
        k = sentry[0].strip()
        v = sentry[1].strip()
                        
        attributes[k] = v

    return attributes
            
