"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (alex.holehouse@wustl.edu, g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""


from . import interface_tools 
from shephard.exceptions import InterfaceException, ProteinException, DomainException
import shephard.exceptions as shephard_exceptions

class _DomainsInterface:

    """
    Class whose sole purpose is to encapsulate and then store
    parsed Domains files. This is a hidden class and is not 
    accessible outside of this file
    
    """

    def __init__(self, filename, delimiter='\t', skip_bad=True):
        """
        Expect files of the following format:

        Unique_ID, start, stop, domain_type, key1:value1, key2:value2, ..., keyn:valuen

        Note that the first four arguments are required, while all of the key:value pairs 
        are optional. Key value must be separated by a ':', but any delimiter (other than ':') 
        is allowed.

        When created, this constructor parses the keyfile to generate a .data class object, 
        which itself maps a uniqueID to a list of domain dictionaries.

        Domain dictionaries have the following key-value pairs

        REQUIRED:
        start                : int (domain start position)
        end                  : int (domain end position)
        domain_type          : string (domain type)

        OPTIONAL:
        attributes           : dictionary of arbitrary key-value pairs 
                               that will be associated with the domain


        

        """

        if delimiter == ':':
            raise InterfaceException('When parsing domain file cannot use ":" as a delimiter because this is used to delimit key/value pairs (if provided)')

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
                attributes = {}
            except Exception as e:

                msg = 'Failed parsing file [%s] on line [%i].\n\nException raised: %s\n\nline printed below:\n%s'%(filename, linecount, str(e), line)

                # if we're skipping bad things then...
                if skip_bad:
                    shephard_exceptions.print_warning(msg + "\nSkipping this line...")
                    continue
                else:
                    raise InterfaceException(msg)
            
            # if some key/value pairs were included then parse these out one at a time
            if len(sline) > 4:
                attributes = interface_tools.parse_key_value_pairs(sline[4:], filename, linecount, line)
                                          
            if unique_ID in ID2domain:
                ID2domain[unique_ID].append({'start':start, 'end':end, 'domain_type':domain_type, 'attributes':attributes})
            else:
                ID2domain[unique_ID] =[{'start':start, 'end':end, 'domain_type':domain_type, 'attributes':attributes}]

        self.data = ID2domain



##############################################
##                                          ##
##     PUBLIC FACING FUNCTIONS BELOW        ##
##                                          ##
##############################################


## ------------------------------------------------------------------------
##
def add_domains_from_file(proteome, filename, delimiter='\t', autoname=False, return_dictionary=False, safe=True, skip_bad=True, verbose=True):
    """
    Function that takes a correctly formatted shephard 'domains' file and reads 
    all domains into the passed Proteome.

    Expect Domain files to have the following format:

    One domain per line where:
    
    >>> Unique_ID    start    stop    domain_type    key_1:value_1    key_2:value_2     ...  key_n:value_n

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

    autoname : bool
        If autoname is set to True, this function ensures each domain ALWAYS has a unique
        name - i.e. the allows for multiple domains to be perfectly overlapping in position
        and type. This is generally not going to be required and/or make sense, but having
        this feature in place is useful. In general we want to avoid this as it makes it 
        easy to include duplicates which by default are prevented when autoname = False. 
        Default = False.

    return_dictionary : bool
        If set to true, this function will return the domains dictionary and will NOT add that
        dictionary to the proteome - i.e. the function basically becomes a parser for SHEPHARD-compliant
        domains files. Default = False

    safe : bool
        If set to True then any exceptions raised during the domain-adding process (i.e. after file
        parsing) are acted on. If set to false, exceptions simply mean the domain in question is skipped. 
        Note if set to False, pre-existing domains with the same name would be silently overwritten (although 
        this is not consider an error), while overwriting will trigger an exception in safe=True.
        There are various reasons domain addition could fail (start/end position outside of the 
        protein limits etc) and so if verbose=True then the cause of an exception is also printed to 
        screen. It is highly recommend that if you choose to use safe=False you also set verbose=True. 
        Default = True.
        
    skip_bad : boolean
        Flag that means if bad lines (lines that trigger an exception) are encountered the code 
        will just skip them. By default this is true, which adds a certain robustness to file 
        parsing, but could also hide errors. Note that if lines are skipped a warning will be 
        printed (regardless of verbose flag). Default = True
    
    verbose : boolean
        Flag that defines how 'loud' output is. Will warn about errors on adding domains.

    Returns
    -----------
    None or dict
        If return_dictionary is set to False (default) then this function has no return
        value, but the domains are added to the Proteome object passed as the first argument. If
        return_dictionary is set to True the function returns the parsed domains dictionary without
        adding the newly-read domains to the proteome.
            
            
    """        

    # check first argument is a proteome
    interface_tools.check_proteome(proteome, 'add_domains_from_file (si_domains)')

    # next read in the file
    domains_interface = _DomainsInterface(filename, delimiter, skip_bad=skip_bad)

    if return_dictionary:
        return domains_interface.data

    # finally add the domains from the dictionary generated by the DomainsInterface parser
    add_domains_from_dictionary(proteome, domains_interface.data, autoname=autoname, safe=safe, verbose=verbose)



## ------------------------------------------------------------------------
##
def add_domains_from_dictionary(proteome, domain_dictionary, autoname=False, safe=True, verbose=True):
    """
    Function that takes a correctly formatted Domains dictionary and will add those 
    domains to the proteins in the Proteome.

    Domains dictionaries are key-value pairs, where the key is a unique_ID associated 
    with a given protein, and the value is a list of dictionaries. Each subdictionary has 
    four key-value pairs:

    'start' = start position (int showing start of the domain, starting at 1)
    'end' = end position (int showing end of the domain, inclusive)
    'domain_type' = domain type (string that names the domain)
    'attributes' = dictionary of arbitrary key:value pairings (optional)

    The start and end positions should be locations within the sequence defined by the unique_ID, 
    and if they are out of the sequence bounds this will throw an exception. Domain type is a string
    that names the type of domain. The attributes dictionary is an arbitrary key-value pair dictionary 
    where key-values map an arbitrary key to an arbitrary value (read in as strings).

    In this way, each domain that maps to a give unique_ID will be added. Note the attribute is
    optional.

    Parameters
    ----------
    proteome : Proteome object
        Proteome object to which domains will be added

    domain_dictionary : dict
        Dictionary that maps unique_IDs to a list of one or more domain dictionaries

    autoname : bool
        If autoname is set to true, this function ensures each domain ALWAYS has a unique
        name - i.e. the allows for multiple domains to be perfecly overlapping in position
        and type. This is generally not going to be required and/or make sense, but having
        this feature in place is useful. In general we want to avoid this as it makes it 
        easy to include duplicates which by default are prevented when autoname = False. 
        Default = False.
    
    safe : bool
        If set to True then any exceptions raised during the Domain-adding process are acted
        on. If set to False, exceptions simply mean the domain in question is skipped. 
        Note if set to False, pre-existing Domains with the same name would be silently overwritten (although 
        this is not consider an error), while overwriting will trigger an exception in safe=True
        There are various reasons Domain addition could fail (start/end position outside of the 
        protein limits etc.) and so if verbose=True then the cause of an exception is also printed to 
        screen. It is highly recommend that if you choose to use safe=False you also set verbose=True. 
        Default = True.
    
    verbose : bool
        Flag that defines how 'loud' output is. Will warn about errors on adding domains.

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

                start       = domain['start']
                end         = domain['end']
                domain_type = domain['domain_type']
                                    
                try:
                    ad          = domain['attributes']
                except:
                    ad = {}
                
                # try and add the domain...
                try:
                    protein.add_domain(start, end, domain_type, attributes=ad, safe=safe, autoname=autoname)
                except (ProteinException, DomainException) as e:

                    msg='- skipping domain at %i-%i on %s' %(start, end, protein)
                    if safe:
                        shephard_exceptions.print_and_raise_error(msg, e)
                    else:
                        if verbose:
                            shephard_exceptions.print_warning(msg)
                            continue
                

## ------------------------------------------------------------------------
##
def write_domains(proteome, filename, delimiter='\t'):
    """
    Function that writes out domains to file in a standardized format. Note that
    attributes are converted to a string, which for simple attributes is reasonable
    but is not really a viable stratergy for complex objects, although this will 
    not yeild and error.
    
    Parameters
    -----------

    proteome :  Proteome object
        Proteome object from which the domains will be extracted from

    filename : str
        Filename that will be used to write the new domains file

    delimiter : str
        Character (or characters) used to separate between fields. Default is '\t'
        Which is recommended to maintain compliance with default `add_domains_from
        file()` function

    Returns
    --------
    None
        No return type, but generates a new file with the complete set of domains
        from this proteome written to disk.

    """

    with open(filename, 'w') as fh:
        for protein in proteome:
            for d in protein.domains:

                # systematically construct each line in the file 
                line = ''
                line = line + str(protein.unique_ID) + delimiter

                start = d.start
                line = line + str(start) + delimiter

                end = d.end
                line = line + str(end) + delimiter

                domain_type = d.domain_type                
                line = line + str(domain_type)

                if d.attributes:
                    for k in d.attributes:
                        line = line + delimiter
                        line = line + str(k) + ":" + str(d.attribute(k))

                line = line + "\n"

                fh.write('%s'%(line))
