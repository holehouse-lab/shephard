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

MAX_BAD_COUNT  = 10

class _DomainsInterface:

    """
    Class whose sole purpose is to encapsulate and then store
    parsed Domains files. This is a hidden class and is not 
    accessible outside of this file
    
    """

    def __init__(self, filename, delimiter='\t', skip_bad=True, preauthorized_uids=None):
        r"""
        Expect files of the following format:

        Unique_ID, start, stop, domain_type, key1:value1, key2:value2, ..., keyn:valuen

        Note that the first four arguments are required, while all of the 
        key:value pairs are optional. Key value must be separated by a ':', 
        but any delimiter (other than ':') 
        
        is allowed.

        When created, this constructor parses the keyfile to generate a .data 
        class object, 
        which itself maps a uniqueID to a list of domain dictionaries.

        Domain dictionaries have the following key-value pairs

        REQUIRED::
            start                : int (domain start position)
            end                  : int (domain end position)
            domain_type          : string (domain type)

            OPTIONAL:
            attributes           : dictionary of arbitrary key-value pairs 
                                   that will be associated with the domain

        Parameters
        ----------------
        
        filename : str
            Name of the shephard domains file to read

        delimiter : str (default = '\\t')
            String used as a delimiter on the input file. 

        skip_bad : bool (default = True)
            Flag that means if bad lines (lines that trigger an exception) 
            are encountered the code will just skip them. By default this is 
            true, which adds a certain robustness to file parsing, but could 
            also hide errors. Note that if lines are skipped a warning will be 
            printed (regardless of verbose flag). 

        preauthorized_uids : list of str (default = None)
            List of unique_IDs that are allowed to be added to the domains
            dictionary. If None then all domains are allowed. Avoids parsing
            lines that are not needed into the interface objects


        """

        bad_count = 0

        if delimiter == ':':
            raise InterfaceException('When parsing domain file cannot use ":" as a delimiter because this is used to delimit key/value pairs (if provided)')

        with open(filename,'r') as fh:
            content = fh.readlines()


        # convert the preauthorized uids to a set for faster lookup.         
        if preauthorized_uids is not None:
            preauthorized_uids = set(preauthorized_uids)
            
        ID2domain = {}

        linecount = 0
        for line in content:

            linecount = linecount + 1

            # skip comment lines
            if interface_tools.is_comment_line(line):
                continue

            sline = line.strip().split(delimiter)
            
            try:
                unique_ID = sline[0].strip()

                # check if UID associated with this line is found in the
                # preauthorized list. If not, skip this line 
                if preauthorized_uids is not None and unique_ID not in preauthorized_uids:
                    continue
                
                start = int(sline[1].strip())
                end = int(sline[2].strip())
                domain_type = sline[3].strip()
                attributes = {}
                
            except Exception as e:

                msg = f'Failed parsing file [{filename}] on line [{linecount}].\n\nException raised: {str(e)}\n\nline printed below:\n{line}'

                # if we're skipping bad things then...
                if skip_bad and bad_count < MAX_BAD_COUNT:
                    bad_count = bad_count + 1
                    shephard_exceptions.print_warning(msg + f"\nSkipping this line (count {bad_count} of {MAX_BAD_COUNT} ...)")
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
    r"""
    Function that takes a correctly formatted shephard 'domains' file and 
    reads all domains into the passed Proteome.
    
    Expect Domain files to have the following format:

    One domain per line where with the format::
    
            1       2    3       4            5            6                   n
        Unique_ID start stop domain_type key_1:value_1 key_2:value_2 ... key_n:value_n

    A couple of key points here:

        * The default delimiter is tabs ('\\t') but this can be changed with 
          the delimiter argument. 
          
        * The first four elements in the each line are required, while all of 
          the key:value pairs are optional

        * Attribute key-value pairs must be separated by a ``:`` character. 
          As a result any column delimiter (other than ``:``) can be used, 
          but ``:`` is reserved for this role
          
    Parameters
    ----------
    proteome : shephard.proteome.Proteome
        Proteome object to which domains will be added

    filename : str
        Name of the shephard domains file to read

    delimiter : str (default = '\\t')
        String used as a delimiter on the input file.

    autoname : bool (default = False)
        If autoname is set to True, this function ensures each domain ALWAYS 
        has a unique name - i.e. the allows for multiple domains to be 
        perfectly overlapping in position and type. This is generally not 
        going to be required and/or make sense, but having this feature in 
        place is useful. In general we want to avoid this as it makes it 
        easy to include duplicates which by default are prevented when 
        autoname = False. 

    return_dictionary : bool, default=False
        If set to true, this function will return the domains dictionary 
        and will NOT add that dictionary to the proteome - i.e. the 
        function basically becomes a parser for SHEPHARD-compliant 
        domains files. 

    safe : bool (default = True)
        If set to True then any exceptions raised during the domain-adding 
        process (i.e. after file parsing) are acted on. If set to false, 
        exceptions simply mean the domain in question is skipped. Note if 
        set to False, pre-existing domains with the same name would be 
        silently overwritten (although this is not consider an error), 
        while overwriting will trigger an exception in safe=True. There 
        are various reasons domain addition could fail (start/end position 
        outside of the protein limits etc) and so if verbose=True then the
        cause of an exception is also printed to screen. It is highly 
        recommend that if you choose to use safe=False you also set 
        verbose=True.       

    skip_bad : bool (default = True)
        Flag that means if bad lines (lines that trigger an exception) are 
        encountered the code will just skip them. By default this is true, 
        which adds a certain robustness to file parsing, but could also 
        hide errors. Note that if lines are skipped a warning will be 
        printed (regardless of verbose flag). 

    verbose : bool (default  = True)
        Flag that defines how 'loud' output is. Will warn about errors on 
        adding domains.

    Returns
    -----------
    None or dict
        If return_dictionary is set to False (default) then this function 
        has no return value, but the domains are added to the Proteome object 
        passed as the first argument. If return_dictionary is set to True 
        the function returns the parsed domains dictionary without
        adding the newly-read domains to the proteome.

    """        

    # check first argument is a proteome
    interface_tools.check_proteome(proteome, 'add_domains_from_file (si_domains)')

    # next read in the file
    domains_interface = _DomainsInterface(filename,
                                          delimiter = delimiter,
                                          skip_bad=skip_bad,
                                          preauthorized_uids = proteome.proteins)

    if return_dictionary:
        return domains_interface.data

    # finally add the domains from the dictionary generated by the
    # DomainsInterface parser
    add_domains_from_dictionary(proteome, domains_interface.data, autoname=autoname, safe=safe, verbose=verbose)



## ------------------------------------------------------------------------
##
def add_domains_from_dictionary(proteome, domain_dictionary, autoname=False, safe=True, verbose=True):
    """
    Function that takes a correctly formatted Domains dictionary and will add 
    those domains to the proteins in the Proteome.
    
    Domains dictionaries are key-value pairs, where the key is a unique_ID 
    associated with a given protein, and the value is a list of dictionaries. 
    Each subdictionary has four key-value pairs::

       * 'start' = start position (int showing start of the domain, starting at 1)

       * 'end' = end position (int showing end of the domain, inclusive)

       * 'domain_type' = domain type (string that names the domain)

       * 'attributes' = dictionary of arbitrary key:value pairings (optional)

    The start and end positions should be locations within the sequence 
    defined  by the unique_ID, and if they are out of the sequence bounds 
    this will throw an exception. Domain type is a string that names the type 
    of domain. The attributes dictionary is an arbitrary key-value pair 
    dictionary where key-values map an arbitrary key to an arbitrary value 
    (read in as strings).
    
    In this way, each domain that maps to a give unique_ID will be added. 
    Note the attribute is optional.
    
    Parameters
    ----------
    proteome : Proteome object
        Proteome object to which domains will be added

    domain_dictionary : dict
        Dictionary that maps unique_IDs to a list of one or more domain 
        dictionaries

    autoname : bool (default = False)
        If autoname is set to true, this function ensures each domain 
        ALWAYS has a unique name - i.e. the allows for multiple domains 
        to be perfecly overlapping in position and type. This is generally 
        not going to be required and/or make sense, but having this feature 
        in place is useful. In general we want to avoid this as it makes it 
        easy to include duplicates which by default are prevented when 
        autoname = False. 
    
    safe : bool (default = True)
        If set to True then any exceptions raised during the Domain-adding 
        process are acted on. If set to False, exceptions simply mean the 
        domain in question is skipped. Note if set to False, pre-existing 
        Domains with the same name would be silently overwritten (although 
        this is not consider an error), while overwriting will trigger an 
        exception in safe=True There are various reasons Domain addition 
        could fail (start/end position outside of the protein limits etc.) 
        and so if verbose=True then the cause of an exception is also printed 
        to screen. It is highly recommend that if you choose to use safe=False 
        you also set verbose=True. 
    
    verbose : bool (default = True)
        Flag that defines how 'loud' output is. Will warn about errors on 
        adding domains.

    Returns
    -----------
    None
        No return value, but domains are added to the Proteome object passed 
        as the first argument.
    
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
def add_domain_attributes_from_file(proteome, filename, delimiter='\t', safe=True, add_new=True, skip_bad=True, verbose=True):
    r"""
    Function that takes a correctly formatted 'domain' files and reads 
    all domain attributes adding them to domains in the passed proteome, 
    if new domains are inclued the add_new flag determins if new domains 
    are added.
    
    The function expects domain attribute files to have the following 
    format:

    One domain defined per line (although the same protein can appear 
    multiple times)::
        
       Unique_ID,  domain_name, key1:value1, key2:value2, ..., keyn:valuen

    A couple of key points here:

    * The default delimiter is tabs ('\\t') but this can be changed with the delimiter argument

    * Key value must be separated by a ':', as a result, any delimiter (other than ':') can be used, but ':' is reserved for this role.

    Parameters
    ------------
    proteome : Proteome Object
        Proteome object to which attributes will be added

    filename : str
        Name of the shephard protein attributes file to read

    delimiter : str (default = '\t')
        String used as a delimiter on the input file. 

    add_new : boolean (default = True)
        If set to True then any new found domains are added to their 
        associated protein. If False any unfound domains are not added 
        and are skipped over. 

        If a new domain is passed that does not have an associated 
        protein in the passed proteome an exception will always be 
        raised regardless of the status of this parameter.

    safe : bool (default = True)
        If set to True then any exceptions raised during the 
        protein_attribute-adding process are acted on. If set to False, 
        exceptions simply mean the protein_attribute in question is skipped.         
        Note if set to False, pre-existing protein_attributes with the same 
        name would be silently overwritten (although this is not consider an 
        error), while overwriting will trigger an exception in safe=True.
        
        The only reason protein attribute addition could fail is if the 
        attribute already exists, so this is effectively a flag to define 
        if pre-existing attributes should be overwritten (False) or not 
        (True).

    skip_bad : bool (default = True)
        Flag that means if bad lines (lines that trigger an exception) are 
        encountered the code will just skip them. By default this is true, 
        which adds a certain robustness to file parsing, but could also hide
        errors. Note that if lines are skipped a warning will be printed 
        (regardless of verbose flag). 
    
    verbose : bool (default = True)
        Flag that defines how 'loud' output is. Will warn about errors on 
        adding attributes.

    Returns
    ----------- 
    None or dict
        If return_dictionary is set to False (default) then this function 
        has no return value, but the protein_attributes are added to the 
        Proteome object passed as the first argument. If return_dictionary 
        is set to True the function returns the parsed domains_dictionary 
        without adding the newly-read protein_attributes to the proteome.
    
    """        
    # check first argument is a proteome
    interface_tools.check_proteome(proteome, 'add_attributes_from_file (si_protein_attributes)')

    # next read in the domain file
    domains_interface = _DomainsInterface(filename, delimiter, skip_bad=skip_bad)

    # finally add the domains from the dictionary generated by the DomainsInterface parser
    add_domain_attributes_from_dictionary(proteome, 
                                           domains_interface.data, 
                                           add_new=add_new,
                                           safe=safe, 
                                           verbose=verbose)

## ------------------------------------------------------------------------
##
def add_domain_attributes_from_dictionary(proteome, domain_dictionary, add_new=True, safe=True, verbose=True):
    """
    Function that takes a correctly formatted Domains dictionary and will 
    add those associated attributes domains to the proteins in the Proteome.
    
    Domains dictionaries are key-value pairs, where the key is a unique_ID 
    associated  with a given protein, and the value is a list of 
    dictionaries. Each subdictionary has four key-value pairs:

       * 'protein'  the unique_ID of the protein for which to domain is associated with
                    
       * 'domain_name' = domain type (string that names the domain)

       * 'attributes'  = dictionary of arbitrary key:value pairings (optional)

    The start and end positions should be locations within the sequence 
    defined by the unique_ID, and if they are out of the sequence bounds this 
    will throw an exception. Domain type is a string that names the type of 
    domain. The attributes dictionary is an arbitrary key-value pair dictionary 
    where key-values map an arbitrary key to an arbitrary value (read in as 
    strings).

    In this way, each domain that maps to a give unique_ID will be added. Note 
    the attribute is optional.

    Parameters
    ----------
    proteome : Proteome object
        Proteome object to which domains will be added

    domain_dictionary : dict
        Dictionary that maps unique_IDs to a list of one or more domain 
        dictionaries.

    add_new : boolean (default = True)
        If set to True then any new found domains are added to their 
        associated protein. If False any unfound domains are not added 
        and are skipped over. If a new domain is passed that does not 
        have an associated protein in the passed proteome an exception 
        will always be raised regardless of the status of this parameter.
    
    safe : bool (default = True)
        If set to True then any exceptions raised during the Domain-adding 
        process are acted on. If set to False, exceptions simply mean the 
        domain in question is skipped. Note if set to False, pre-existing 
        Domains with the same name would be silently overwritten (although 
        this is not consider an error), while overwriting will trigger an 
        exception in safe=True There are various reasons Domain addition 
        could fail (start/end position outside of the protein limits etc.) 
        and so if verbose=True then the cause of an exception is also printed 
        to screen. It is highly recommend that if you choose to use safe=False 
        you also set verbose=True. 
    
    verbose : bool (default = True)
        Flag that defines how 'loud' output is. Will warn about errors 
        on adding domains.

    Returns
    -----------
    None
        No return value, but domains are added to the Proteome object passed 
        as the first argument.
    
    """
    # Note - the safe keyword is actually dealt with in this function in conjunction with the Verbose
    # keyword, so we pass safe=False to the add_domain function and then catch the exception in this
    # function.


    # check first argument is a proteome
    interface_tools.check_proteome(proteome, 'add_domains (si_domains)')
    
    # iterate proteins with new domains 
    for unique_ID in domain_dictionary:

        # check if protein in proteome  
        if unique_ID in proteome:

            # build dict of local domains with domain IDs as keys
            local_domain_dict = {"%s_%i_%i" % (d.domain_type, d.start, d.end): d for d in proteome.protein(unique_ID).domains}

            # iterate new domains
            for new_domain in domain_dictionary[unique_ID]:
                new_domain_ID = "%s_%i_%i" % (new_domain['domain_type'], new_domain['start'], new_domain['end'])
                
                # check if domain is in local domain by ID 
                if new_domain_ID in local_domain_dict:
                    local_domain = local_domain_dict[new_domain_ID]
                    try:
                        ad = new_domain['attributes']
                    except:
                        ad = {}

                    # merge attributes
                    for k, v in ad.items():

                        try:
                            local_domain.add_attribute(k, v, safe=safe)
                        except (ProteinException, DomainException) as e:

                            msg = f"- skipping attribute being added to {unique_ID} for domain type {new_domain['domain_type']}, with start={new_domain['start']} and end={new_domain['end']})"
                            if safe:
                                shephard_exceptions.print_and_raise_error(msg, e)
                            else:
                                if verbose:
                                    shephard_exceptions.print_warning(msg)
                                    continue

                    # move on to next new domain
                    continue

                # add new domain if new domain not found and flag is true 
                if add_new:
                    try:
                        ad = new_domain['attributes']
                    except:
                        ad = {}
                    
                    # try and add the domain...
                    try:
                        proteome.protein(unique_ID).add_domain(new_domain['start'], new_domain['end'], new_domain['domain_type'], attributes=ad, safe=safe)
                    except (ProteinException, DomainException) as e:

                        msg='- skipping domain at %i-%i on %s' %(new_domain['start'], new_domain['end'], proteome.protein(unique_ID))
                        if safe:
                            shephard_exceptions.print_and_raise_error(msg, e)
                        else:
                            if verbose:
                                shephard_exceptions.print_warning(msg)
                                continue


    
## ------------------------------------------------------------------------
##
def write_domains(proteome, filename, delimiter='\t', domain_types=None):
    r"""
    Function that writes out domains to a SHEPHARD domains file. Note that
    attributes are converted to a string, which for simple attributes is 
    reasonable but is not really a viable stratergy for complex objects, 
    although this will not yeild and error.
            
    Parameters
    -----------

    proteome :  Proteome object
        Proteome object from which the domains will be extracted from

    filename : str
        Filename that will be used to write the new domains file

    delimiter : str (default = '\\t')
        Character (or characters) used to separate between fields. 
        Default is '\t' Which is recommended to maintain compliance 
        with default `add_domains_from_file()` function.

    domain_types : list (default None)
        Lets you define a list of one or more domain types that will
        be written out. Domain types are passed as strings which should
        map to named domain types in the Proteome.
        
    Returns
    --------
    None
        No return type, but generates a new file with the complete 
        set of domains from this proteome written to disk.
        

    """

    # added so that we ensure domain_types is a list if passed
    if domain_types is not None:
        if type(domain_types) is not list:
            raise InterfaceException('When passing a domain_type this must be a list')
        

    with open(filename, 'w') as fh:
        for protein in proteome:
            for d in protein.domains:

                # if domain_types is passed check if each domain
                # is found in the list
                if domain_types is not None:
                    if d.domain_type not in domain_types:
                        continue

                line = __build_domain_line(d, delimiter)
                
                fh.write(line)

## ------------------------------------------------------------------------
##
def write_domains_from_list(domain_list, filename, delimiter='\t'):
    r"""
    Function that writes out domains to a SHEPHARD domains file from a list
    of Domain objects. 
    Note that attributes are converted to a string, which for simple 
    attributes is reasonable but is not really a viable stratergy for 
    complex objects, although this will not yeild and error.
            
    Parameters
    -----------

    domain_list : List of Domain objects
        List of domain objects which will be written

    filename : str
        Filename that will be used to write the new domains file

    delimiter : str (default = '\\t')
        Character (or characters) used to separate between fields. Default is 
        '\\t' which is recommended to maintain compliance with default 
        `add_domains_from_file()` function
       
    Returns
    --------
    None
        No return type, but generates a new file with the complete set of 
        domains from this proteome written to disk.

    """

    # first check if items in the list are Domain objects
    for d in domain_list:
        interface_tools.check_domain(d, 'write_domains_from_list')


    with open(filename, 'w') as fh:
        for d in domain_list:

            line = __build_domain_line(d, delimiter)

            fh.write(line)


## ------------------------------------------------------------------------
##
def __build_domain_line(d, delimiter):
    """
    Internal function that takes a Domain object and returns a line that can
    be written to a Domains file. This is called internally by functions that
    write Domains.

    Parameters
    ----------------------
    d : shephard.Domain
        Domain object being converted to a string

    delimiter : str (default = '\\t')
        Character (or characters) used to separate between fields. 
        Default is the tab character ('\\t'), which is recommended to 
        maintain compliance with default SHEPHARD file-reading functions.     

    Returns
    --------------
    str
        Returns a string that is ready to be written to file

    """

    # systematically construct each line in the file 
    line = ''
    line = line + str(d.protein.unique_ID) + delimiter

    start = d.start
    line = line + str(start) + delimiter

    end = d.end
    line = line + str(end) + delimiter

    domain_type = d.domain_type    
            
    # note last required element has no trailing delimiter
    line = line + str(domain_type) 

    if d.attributes:
        for k in d.attributes:

            # 
            atrbt = interface_tools.full_clean_string(d.attribute(k))
            line = line + delimiter + f"{k}:{atrbt}" 

    line = line + "\n"

    return line
