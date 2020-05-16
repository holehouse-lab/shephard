"""
SHEPHARD: 
Sequence-based Hierachical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

from .interface_exceptions import InterfaceException
from . import interface_tools 
from shephard.exceptions import ProteinException

class _DomainsInterface:

    """
    Class whose sole purpose is to encapsulate and then store
    parsed Domains files. This is a hidden class and is not 
    accessible outside of this file
    
    """

    def __init__(self, filename, delimiter='\t', skip_bad=True):
        """
        Expect files of the followin format:

        Unique_ID, start, stop, domain_type, key1:value1, key2:value2, ..., keyn:valuen

        Note that the first four arguments are required, while all of the key:value pairs 
        are optional. Key value must be separated by a ':', but any delimiter (other than ':') 
        is allowed

        """

        if delimiter == ':':
            raise InterfaceException('When parsing domain file cannot use ":" as a delimeter because this is used to delimit key/value pairs (if provided)')

        with open(filename,'r') as fh:
            content = fh.readlines()

        ID2domain={}

        linecount=0
        for line in content:

            linecount=linecount+1

            sline = line.strip().split(delimiter)
            
            try:
                unique_ID = sline[0].strip()
                start = int(sline[1].strip())
                end = int(sline[2].strip())
                domain_type = sline[3].strip()
                attribute_dictionary = {}
            except Exception:

                msg = 'Failed parsing file [%s] on line [%i]... line printed below:\n%s'%(filename, linecount, line)

                # if we're skipping bad things then...
                if skip_bad:
                    print('Warning: %s'%msg)
                    print('Skipping this line...')
                    continue
                else:
                    raise InterfaceException('Error: %s'%msg)
            
            # if some key/value pairs were included then parse these out one at a time
            if len(sline) > 4:
                attribute_dictionary = interface_tools.parse_key_value_pairs(sline[4:], filename, linecount, line)
                                          
            if unique_ID in ID2domain:
                ID2domain[unique_ID].append([start, end, domain_type, attribute_dictionary])
            else:
                ID2domain[unique_ID] =[[start, end, domain_type, attribute_dictionary]]

        self.data = ID2domain



##############################################
##                                          ##
##     PUBLIC FACING FUNCTIONS BELOW        ##
##                                          ##
##############################################


## ------------------------------------------------------------------------
##
def add_domains_from_file(proteome, filename, delimiter='\t', autoname=False, safe=True, skip_bad=True, verbose=True):
    """
    Function that takes a correctly formatted shephard 'domains' file and reads 
    all domains into the passed proteome.

    Expect Domain files to have the following format:

    One domain per line where:
    
    >>> Unique_ID    start    stop    domain_type    key_1:value_1    key_2:value_2, ...,     key_n:value_n

    A couple of key points here:
    - The default delimiter is tabs ('\t') but this can be changed with the delimiter argument
    - The first four arguments are required, while all of the key:value pairs are optional
    - Key value must be separated by a ':'. As a result any column delimiter (other than ':') can be 
      used, but ':' is reserved for this role
    
    Parameters
    ----------
    proteome : Proteome Object
        Proteome object to which domains will be added

    filename : str
        Name of the shephard domains file to read

    delimiter : str 
        String used as a delimiter on the input file. Default = '\t'

    autoname : boolean
        If autoname is set to true, this function ensures each domain ALWAYS has a unique
        name - i.e. the allows for multiple domains to be perfecly overlapping in position
        and type. This is generally not going to be required and/or make sense, but having
        this feature in place is useful. In general we want to avoid this as it makes it 
        easy to include duplicates which by default are prevented when autoname = False. 
        Default = False.
    
    safe : boolean 
        If set to True over-writing domains will raise an exception. If False, overwriting
        a domain will silently over-write. Default = True.

    skip_bad : boolean
        Flag that means if bad lines (lines that trigger an exception) are encountered the code 
        will just skip them. By default this is true, which adds a certain robustness to file 
        parsing, but could also hide errors. Note that if lines are skipped a warning will be 
        printed (regardless of verbose flag). Default = True
    

    verbose : boolean
        Flag that defines how 'loud' output is. Will warn about errors on adding domains.

    Returns
    -----------
    None
        No return value, but domains are added to the Proteome object passed as the first argument
            
    """        


    # check first argument is a proteome
    interface_tools.check_proteome(proteome, 'add_domains_from_file (si_domains)')

    # next read in the file
    domains_interface = _DomainsInterface(filename, delimiter, skip_bad=skip_bad)

    # finally add the domains from the dictionary generated by the DomainsInterface parser
    add_domains_from_dictionary(proteome, domains_interface.data, autoname=autoname, safe=safe, verbose=verbose)



## ------------------------------------------------------------------------
##
def add_domains_from_dictionary(proteome, domain_dictionary, autoname=False, safe=True, verbose=True):
    """
    Function that takes a correctly formatted Domains dictionary and will add those 
    domains to the proteins in the proteome.

    Domains dictionaries are key-value pairs, where the key is a unique_ID associated 
    with a given protein, and the value is a list of lists. Each sublist has four positions

    [0] = start position
    [1] = end position 
    [2] = domain type
    [3] = attribute dictionary

    The start and end positions should be locations within the sequence defined by the unique_ID, 
    and if they are out of the sequence bounds this will throw an exception. Domain type is a string
    that nams the type of domain. The attribute dictionary is an arbitrary key-value pair dictionary 
    where key-values map an arbitrary key to an aribitrary value. 

    In this way, each domain that maps to a give unique_ID will be added. Note the attribute is
    optional.

    Parameters
    ----------
    proteome : Proteome Object
        Proteome object to which domains will be added

    domain_dictionary : dict
        Dictionary that maps unique_IDs to domain lists [start, end, type, attribute_dictionary].


    autoname : boolean
        If autoname is set to true, this function ensures each domain ALWAYS has a unique
        name - i.e. the allows for multiple domains to be perfecly overlapping in position
        and type. This is generally not going to be required and/or make sense, but having
        this feature in place is useful. In general we want to avoid this as it makes it 
        easy to include duplicates which by default are prevented when autoname = False. 
        Default = False.
    
    safe : boolean 
        If set to True over-writing domains will raise an exception. If False, overwriting
        a domain will silently over-write. Default = True.

    skip_bad : boolean
        Flag that means if bad lines (lines that trigger an exception) are encountered the code 
        will just skip them. By default this is true, which adds a certain robustness to file 
        parsing, but could also hide errors. Note that if lines are skipped a warning will be 
        printed (regardless of verbose flag). Default = True.

    Returns
    -----------
    None
        No return value, but domains are added to the Proteome object passed as the first argument
    
    """
    # Note - the safe keyword is actually dealt with in this function in conjunction with the Verbose
    # keyword, so we pass safe=False to the add_domain function and then catch the exception in this
    # function.


    # check first argument is a proteome
    interface_tools.check_proteome(proteome, 'add_domains (si_domains)')
    
    for protein in proteome:
        if protein.unique_ID in domain_dictionary:
            for domain in domain_dictionary[protein.unique_ID]:

                # extract
                start       = domain[0]
                end         = domain[1]
                domain_type = domain[2]
                ad          = domain[3]
                
                # try and add the domain...
                try:
                    protein.add_domain(start, end, domain_type, attribute_dictionary=ad, safe=False, autoname=autoname)
                except ProteinException as e:

                    msg='- skippiing domain at %i-$i on %s' %(start, end, protein)
                    if safe:
                        # if safe=True fail on any errors
                        print('Error %s' %(msg)) 
                        raise e
                    else:
                        if verbose:
                            print('Warning %s' %(msg))
                        continue



                

## ------------------------------------------------------------------------
##
def write_domains(proteome, filename, delimiter='\t', skip_bad=False):
    """
    Function that writes out domains to file in a standardized format.
    
    <TO DO>

    """

    with open(filename, 'w') as fh:
        for protein in proteome:
            for d in protein.domains:
                start = protein.domain(d).start
                end = protein.domain(d).end
                domain_type = protein.domain(d).domain_type
                fh.write('%s%s%i%s%i%s%s\n' % (protein.unique_ID, 
                                               delimiter, 
                                               start,
                                               delimiter,
                                               end,
                                               delimiter,
                                               domain_type))
    
