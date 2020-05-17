"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

from . import exceptions
from shephard.exceptions import ShephardException


## ------------------------------------------------------------------------
##
def valid_keyword(keywordname, kw, kwlist):
    if kw not in kwlist:
        raise ShephardException("Provided keyword '%s' for argument %s but expected one of %s" % (kw, keywordname, str(kwlist)))


## ------------------------------------------------------------------------
##
def cast_or_none(value, cast_type):
    """
    Function that takes a value and a type (cast_type). If the value is None
    the fuction returns None. Otherwise the the cast_type will try and cast 
    the value. 
    
    Parameters
    ------------
    value : unknown
        Some kind of variable that can (hopefully) by 
    
    cast_type : some type of castable type
        A Python type that could be used to cast a variable, such as str, int,
        float etc.


    Returns
    -------
    
    None or casted type
    """

    if value is None:
        return None
    else:
    
        # envelope this in a try/except block to help the user out
        try:
            return cast_type(value)
        except Exception as e:
            raise ShephardException('In cast_or_none unable to cast passed variable, error below:\n%s' % str(e))
            


## ------------------------------------------------------------------------
##
def string_to_list_of_strings(inval):
    """
    Function that converts a string or list of
    strings to a list of strings. 

    """

    if type(inval) == str:
        return [inval]
    elif type(inval) == list:
        return inval

    else:
        raise exceptions.ShephardException('Expected either a single string or a list of strings')



## ------------------------------------------------------------------------
##
def numerical_average(l):
    return numerical_sum(l)/len(l)


## ------------------------------------------------------------------------
##
def numerical_sum(l):

    t = 0

    for i in l:
        t = t + i

    return t



