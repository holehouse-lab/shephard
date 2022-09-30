"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (alex.holehouse@wustl.edu, g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

from shephard.exceptions import InterfaceException


## ------------------------------------------------------------------------
##
def check_protein(p, function_name):
    """
    Function that takes takes in some object and tests if its a Protein
    object (or not). If yes returns None, but if no raise an exception with 
    an error message that includes the name of the parent function that is 
    calling this validation function

    Parameters
    --------------
    p : Protein object 
        Unknown object that will be tested to see if it's a Proteome object
        or not

    function_name : string
        Passed string that makes it easy for the user to debug which function 
        has failed

    Returns
    ---------
    None
        No returns but the function will raise an exception if the passed object 
        is not a Proteome object
    

    """
    if "<class 'shephard.protein.Protein'>" == str(p.__class__):
        return None
    else:
        raise InterfaceException(f'First argument passed to function {function_name} was not a Protein object')


## ------------------------------------------------------------------------
##
def check_proteome(p, function_name):
    """
    Function that takes takes in some object and tests if its a Proteome
    object (or not). If yes returns None, but if no raise an exception with 
    an error message that includes the name of the parent function that is 
    calling this validation function

    Parameters
    --------------
    p : Proteome object 
        Unknown object that will be tested to see if it's a Proteome object
        or not

    function_name : string
        Passed string that makes it easy for the user to debug which function 
        has failed

    Returns
    ---------
    None
        No returns but the function will raise an exception if the passed object 
        is not a Proteome object
    

    """
    if "<class 'shephard.proteome.Proteome'>" == str(p.__class__):
        return None
    else:
        raise InterfaceException(f'First argument passed to function {function_name} was not a Proteome object')

## ------------------------------------------------------------------------
##
def check_domain(d, function_name):
    """
    Function that takes takes in some object and tests if its a Domain
    object (or not). If yes returns None, but if no raise an exception with 
    an error message that includes the name of the parent function that is 
    calling this validation function

    Parameters
    --------------
    d : Domain object 
        Unknown object that will be tested to see if it's a Domain object
        or not

    function_name : string
        Passed string that makes it easy for the user to debug which function 
        has failed

    Returns
    ---------
    None
        No returns but the function will raise an exception if the passed object 
        is not a Domain object
    

    """
    if "<class 'shephard.domain.Domain'>" == str(d.__class__):
        return None
    else:
        raise InterfaceException(f'First argument passed to function {function_name} was not a Domain object')


## ------------------------------------------------------------------------
##
def check_track(t, function_name):
    """
    Function that takes takes in some object and tests if its a Track
    object (or not). If yes returns None, but if no raise an exception with 
    an error message that includes the name of the parent function that is 
    calling this validation function

    Parameters
    --------------
    t : Track object 
        Unknown object that will be tested to see if it's a Track object
        or not

    function_name : string
        Passed string that makes it easy for the user to debug which function 
        has failed

    Returns
    ---------
    None
        No returns but the function will raise an exception if the passed object 
        is not a Track object
    

    """
    if "<class 'shephard.track.Track'>" == str(t.__class__):
        return None
    else:
        raise InterfaceException(f'First argument passed to function {function_name} was not a track object')


## ------------------------------------------------------------------------
##
def check_site(s, function_name):
    """
    Function that takes takes in some object and tests if its a Site
    object (or not). If yes returns None, but if no raise an exception with 
    an error message that includes the name of the parent function that is 
    calling this validation function

    Parameters
    --------------
    s : Site object 
        Unknown object that will be tested to see if it's a Site object
        or not

    function_name : string
        Passed string that makes it easy for the user to debug which function 
        has failed

    Returns
    ---------
    None
        No returns but the function will raise an exception if the passed object 
        is not a Site object
    

    """
    if "<class 'shephard.site.Site'>" == str(s.__class__):
        return None
    else:
        raise InterfaceException(f'First argument passed to function {function_name} was not a site object')
    


## ------------------------------------------------------------------------
##
def full_clean_string(instring):
    r"""
    Wrapper function that takes in a string (or a variable that can be
    cast to a string) and ensure it contains neither tab nor colon characters,
    both of which would render a attribute or free-form text string in a 
    SHEPHARD file invalid.

    This function actually calls `clean_string()` twice with parameters to
    replace, ':' and '\\t' characters.

    Parameters
    -----------
    instring : {str, str-castable}
        String (or string castable) variable to be checked

    Returns
    ---------
    str
        Returns a string which will be essentially identical to the input 
        string but cleaned up to remove tab and colon characteracs .

    """

    s = clean_string(instring)
    s = clean_string(s, ':', '-')
    return s



## ------------------------------------------------------------------------
##
def clean_string(instring, delimiter='\t', replace_char=' '):
    r"""
    Function that takes in a variable  that is to be written to a SHEPHARD 
    complient file and ensures it has no delimiter characters in it. This 
    avoids the scenario where - for example - your protein name has a 
    tab in it which introduces a new field into the SHEPHARD data file.

    In general, we assume instring will be a string (str). However, if it
    is not this function also casts it to a string, so that does not need
    to be done prior to a variable being passed in.

    Furthermore, we'd generally recommend any string that will be written
    out as a user-inputable string to a SHEPHARD file be passed through
    this

    If such delimiter characters are found they're replaced by a replace_char, 
    which by default is just a space    
    
    Parameters
    ---------------
    instring : str
        String (or string castable) variable to be checked

    delimiter : str, default='\\t'
        Character or string that we do NOT want to find in the string

    replace_char : str, default=' '
        Character or string that, if a delimiter is found, will 
        replace the  delimiter

    Returns
    ---------
    str
        Returns a string which will be essentially identical to the input 
        string but cleaned up to remove any bad delimiters if they exist.
        

    """

    if type(instring) is not str:
        instring = str(instring)

    return instring.replace(delimiter, replace_char)



## ------------------------------------------------------------------------
##
def parse_key_value_pairs(split_line, filename, linecount, line):
    """
    Helper function for parsing input files that have attributes (for 
    Domains and Sites).

    The function takes in a list of key-value pairs (where key-values 
    are split by a ':' symbol) and returns a parsed dictionary. Note 
    that values will always be strings.
    
    Parameters
    -------------
    split_line : list of strings
        Each string in the list should have the format <KEY> : <VALUE> 

    filename : string
        Name of the file the calling function is parsing. Only used when 
        raising an exception.

    linecount : int
        Current line number the file processing is on. Only used when 
        raising an exception.

    line : string
        Full line that the file processing is on. Again, only used when 
        raising an exception

    """

    attributes = {}

    for idx in range(len(split_line)):
        try:
            sentry = split_line[idx].split(':')
        except Exception:
                
            # should update this to also display the actual error...
            raise InterfaceException(f'Failed parsing key-value pairs in file [{filename}] on line [{linecount}]... line printed below:\n{line}.\nPassed split_line:[{split_line}]')
            
        # note we enclose this in a try/catch block
        try:
            k = sentry[0].strip()
            v = sentry[1].strip()
                        
            attributes[k] = v
        except IndexError:
            raise InterfaceException(f'Failed parsing key-value pairs in file [{filename}] on line [{linecount}]... line printed below:\n{line}.\nPassed split_line:[{split_line}]')


    return attributes


## ------------------------------------------------------------------------
##
def is_comment_line(line):
    """
    Function that checks if a line should be skipped because it's a comment 
    line or not.

    Parameters
    -------------
    line : str
        A line from an input file

    Returns
    -----------
    bool
        If the line starts with a '#' character then the function returns
        True, else returns False


    """
    line = line.strip()

    if line[0] == '#':
        return True
        

            
