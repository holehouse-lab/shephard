"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (alex.holehouse@wustl.edu, g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

from . import interface_tools 
import shephard.exceptions as shephard_exceptions
from shephard.exceptions import InterfaceException, ProteinException, SiteException

MAX_BAD_COUNT  = 10

class _SitesInterface:

    def __init__(self, filename, delimiter='\t', skip_bad=True, preauthorized_uids=None):
        """
        Expect files of the following format:
        
        A SHEPHARD sites file is a tab (or other) delineated file where each 
        line has the following convention::
    
               1        2          3       4      5   [      6            7        ...     n         ] 
            Unique_ID position site_type symbol value [key_1:value_1 key_2:value_2 ... key_n:value_n ]
    
        Each line has six required values and then can have as many key:value pairs as may be
        desired.

        Note that the first four arguments are required, while all of the 
        key:value pairs are optional. Key value must be separated by a ':', 
        but any delimiter (other than ':') is allowed. 

        Parameters
        ----------------
        
        filename : str
            Name of the shephard domains file to read

        delimiter : str (default = \t)
            String used as a delimiter on the input file. 

        skip_bad : bool (default = True)
            Flag that means if bad lines (lines that trigger an exception) 
            are encountered the code will just skip them. By default this is 
            true, which adds a certain robustness to file parsing, but could 
            also hide errors. Note that if lines are skipped a warning will be 
            printed (regardless of verbose flag). 

        preauthorized_ids : list of str (default = None)
            List of unique_IDs that are allowed to be added to the sites
            dictionary. If None then all sites are allowed. Avoids parsing
            lines that are not needed into the interface objects

        """

        bad_count = 0

        if delimiter == ':':
            raise InterfaceException('When parsing site file cannot use ":" as a delimeter because this is used to delimit key/value pairs (if provided)')


        with open(filename,'r') as fh:
            content = fh.readlines()

        # convert the preauthorized uids to a set for faster lookup
        if preauthorized_uids is not None:
            preauthorized_uids = set(preauthorized_uids)

        ID2site = {}
        
        linecount=0
        for line in content:

            linecount = linecount + 1

            # skip comment lines
            if interface_tools.is_comment_line(line):
                continue

            sline = line.strip().split(delimiter)

            try:
                unique_ID = sline[0].strip()

                # check if UID associated with this line is found in the
                # preauthorized list. If  not then skip this line
                if preauthorized_uids is not None and unique_ID not in preauthorized_uids:
                    continue
                
                position = int(sline[1].strip())
                site_type = sline[2].strip()
                symbol = sline[3].strip()

                # this enables the value to be None if you
                # write a symbol where there's no value associated
                # with a site
                tmp = sline[4].strip()
                if tmp == 'None':
                    value = None
                else:
                    value = float(tmp)

                attributes = {}
                
            except Exception as e:
                msg = f'Failed parsing file [{filename}] on line [{linecount}].\n\nException raised: {str(e)}\n\nline printed below:\n{line}'

                # should update this to also display the actual error...
                if skip_bad and bad_count < MAX_BAD_COUNT:
                    bad_count = bad_count + 1
                    shephard_exceptions.print_warning(msg + f"\nSkipping this line (count {bad_count} of {MAX_BAD_COUNT} ...)")
                    continue
                else:
                    raise InterfaceException(msg)

            # if there's more parse attribute dictionary entries
            if len(sline) > 5:
                attributes = interface_tools.parse_key_value_pairs(sline[5:], filename, linecount, line)

            if unique_ID in ID2site:
                ID2site[unique_ID].append({'position':position, 'site_type':site_type, 'symbol':symbol, 'value':value, 'attributes':attributes})
            else:
                ID2site[unique_ID] =[{'position':position, 'site_type':site_type, 'symbol':symbol, 'value':value, 'attributes':attributes}]

        self.data = ID2site



##############################################
##                                          ##
##     PUBLIC FACING FUNCTIONS BELOW        ##
##                                          ##
##############################################


## ------------------------------------------------------------------------
##
def add_sites_from_file(proteome, filename, delimiter='\t', return_dictionary=False, safe=True, skip_bad=True, verbose=True):
    r"""
    Function that provides the user-facing interface for reading correctly 
    configured SHEPHARD sites files and adding those sites to the proteins 
    of interest.
    
    
    A SHEPHARD sites file is a tab (or other) delineated file where each 
    line has the following convention::
    
          1        2          3       4      5   [      6            7        ...     n         ] 
       Unique_ID position site_type symbol value [key_1:value_1 key_2:value_2 ... key_n:value_n ]
    
    Each line has six required values and then can have as many key:value pairs as may be
    desired.


    Parameters
    ----------
    proteome : Proteome
        Proteome object to which we're adding sites. Note that ONLY sites 
        for which a protein is found will be used. Protein-Site 
        cross-referencing is done using the protein's unique_ID which 
        should be the key used in the sites_dictionary

    filename : str
        Name of the shephard site file to be read

    delimiter : str (default = '\\t')
        String used as a delimiter on the input file. 

    return_dictionary : bool, default=False
        If set to true, this function will return the sites dictionary 
        and will NOT add that dictionary to the proteome - i.e. the 
        function basically becomes a parser for SHEPHARD-compliant        
        sites files. 

    safe : bool (default = True)
        If set to True then any exceptions raised during the site-adding 
        process (i.e. after file parsing) are acted on. If set to False, 
        exceptions simply mean the site in question is skipped. There are 
        various reasons site addition could fail (e.g. site falls outside 
        of protein position so if verbose=True then the cause of an exception 
        is also printed to screen. It is highly recommend that if you choose 
        to use safe=False you also set verbose=True. Default = True.
        
    skip_bad : bool (default = True)
        Flag that means if bad lines (lines that trigger an exception) are 
        encountered the code will just skip them. By default this is true, 
        which adds a certain robustness to file parsing, but could also hide 
        errors. Note that if lines are skipped a warning will be printed 
        (regardless of verbose flag). 

    verbose : bool (default = True)
        Flag that defines how 'loud' output is. Will warn about errors 
        on adding sites.
        
    Returns
    ---------
    None or dict
        If return_dictionary is set to False (default) then this function 
        has no return value, but the sites are added to the Proteome object 
        passed as the first argument. If return_dictionary is set to True 
        the function returns the parsed sites dictionary without adding the 
        newly-read sites to the proteome.

    """

    # check first argument is a proteome
    interface_tools.check_proteome(proteome, 'add_sites_from_file (si_sites)')

    # build the SitesInterface object
    sites_interface = _SitesInterface(filename,
                                      delimiter=delimiter,
                                      skip_bad=skip_bad,
                                      preauthorized_uids = proteome.proteins)

    if return_dictionary:
        return sites_interface.data


    # finally add the site from the dictionary generated by the
    # SitesInterface parser
    add_sites_from_dictionary(proteome, sites_interface.data, safe, verbose)



## ------------------------------------------------------------------------
##
def add_sites_from_dictionary(proteome, sites_dictionary, safe=True, verbose=False):
    """
    Function that takes a correctly formatted Sites dictionary and will add 
    those Sites to the proteins in the Proteome.
    
    Sites dictionaries are key-value pairs, where the key is a unique_ID 
    associated with a given Protein, and the value is a list of dictionaries. 
    Each subdirectionay has the following elements::
    
        'position'   = site position
        'site_type'  = site type
        'symbol'     = site symbol 
        'value'      = site value 
        'attributes' = site attribute dictionary

    In this way, each site that maps to a give unique_ID will be added to 
    the associated protein. The use of a list of dictionaries (as opposed
    to a simple unique_ID:site_dictionary pairing) means multiple sites 
    for a single protein can be added at once.

    Parameters
    -------------

    proteome : Proteome
        Proteome object to which we're adding sites. Note that ONLY sites 
        for which a protein is found will be used. Protein:Site 
        cross-referencing is done using the protein's unique_ID        
        which should be the key used in the sites_dictionary

    sites_dictionary : dict
        A sites dictionary (defined above) is dictionary that maps a 
        unique_ID back to a list of dictionaries, where each 
        subdictionay has five elements, desribed above.

        Recall the only type-specific values (position and value) are 
        cast automatically when a site is added by the Protein object, 
        so there is no need to do that in this function too.

        Extra key-value paris in each sub-dictionary are ignored

    safe : bool (default = True)
        If set to True then any exceptions raised during the site-adding 
        process are acted on. If set to false, exceptions simply mean the 
        site in question is skipped. There are various reasons site addition 
        could fail (notably position of the site is outside of the protein 
        limits) and so if verbose=True then the cause of an exception is
        also printed to screen. It is highly recommend that if you choose to
        use safe=False you also set verbose=True

    verbose : bool (default = False)
        Flag that defines how 'loud' output is. Will warn about errors on 
        adding sites.

    Returns
    ---------
    None
        No return value, but adds all of the passed sites to the protein
    
    """
    
    for protein in proteome:
        if protein.unique_ID in sites_dictionary:
            for site in sites_dictionary[protein.unique_ID]:

                try:
                    position = site['position']
                    site_type = site['site_type']
                    symbol = site['symbol']
                    value = site['value']
                    try:
                        ad = site['attributes'] 
                    except:
                        ad = {}
                except Exception:
                    raise InterfaceException('When sites dictionary for key [%s] was unable to extract five distinct parametes. Entry is:\n%s\n'% (protein.unique_ID, site))

                # assuming we can read all five params try and add the site
                try:
                    protein.add_site(position, site_type, symbol, value, attributes = ad)


                except ProteinException as e:
                    msg='- skipping site %s at %i on %s' %(site_type, position, protein)
                    if safe:
                        shephard_exceptions.print_and_raise_error(msg, e)
                    else:
                        if verbose:
                            shephard_exceptions.print_warning(msg)
                            continue
  

                  
## ------------------------------------------------------------------------
##
def write_sites(proteome, filename, delimiter='\t', site_types=None):
    r"""
    Function that writes out sites to file in a standardized format. Note 
    that attributes are converted to a string, which for simple attributes 
    is reasonable but is not really a viable stratergy for complex objects, 
    although this will not yeild and error.

    If a site_types list is provided, only site_types that match to
    strings in this list are written out.
    
    Parameters
    -----------
    proteome :  Proteome
        Proteome object from which the sites will be extracted from

    filename : str
        Filename that will be used to write the new sites file

    site_type : str (default = None)
        If provided, this is an identifier that allows you to specificy 
        a specific site type to write out.

    delimiter : str (default = '\\t')
        Character (or characters) used to separate between fields. 
        Default is the tab character ('\\t'), which is recommended to 
        maintain compliance with default SHEPHARD file-reading functions.     

    Returns
    --------
    None
        No return type, but generates a new file with the complete set of 
        sites from this proteome written to disk.
        

    """

    # added so that we ensure site_types is a list if passed
    if site_types is not None:
        if type(site_types) is not list:
            raise InterfaceException('When passing a site_type this must be a list')


    with open(filename, 'w') as fh:
        for protein in proteome:
            for s in protein.sites:

                # if we're using site_types and the current sites 
                if site_types is not None:
                    if s.site_type not in site_types:
                        continue

                # build a line 
                # if the passed parameter site_types is being
                # used
                line = __build_site_line(s, delimiter)

                fh.write(f"{line}")



## ------------------------------------------------------------------------
##
def write_sites_from_list(site_list, filename, delimiter='\t'):
    r"""
    Function that writes out sites to a SHEPHARD sites file from a list
    of Site objects. 
    Note that attributes are converted to a string, which for simple 
    attributes is reasonable but is not really a viable stratergy for 
    complex objects, although this will not yeild and error.
            
    Parameters
    -----------

    site_list : List of Site objects
        List of site objects which will be written

    filename : str
        Filename that will be used to write the new sites file

    delimiter : str (default = '\\t')
        Character (or characters) used to separate between fields. Default is 
        '\\t' which is recommended to maintain compliance with default 
        `add_sites_from_file()` function
       
    Returns
    --------
    None
        No return type, but generates a new file with the complete set of 
        sites from this proteome written to disk.

    """

    # first check if items in the list are site objects
    for s in site_list:
        interface_tools.check_site(s, 'write_sites_from_list')

    with open(filename, 'w') as fh:

        # for each site in the list
        for s in site_list:

            # build a line 
            # if the passed parameter site_types is being
            # used
            line = __build_site_line(s, delimiter)

            fh.write(f"{line}")



## ------------------------------------------------------------------------
##
def __build_site_line(s, delimiter):
    """
    Internal function that takes a Site object and returns a line that can
    be written to a Sites file. This is called internally by functions that
    write Sites.

    Parameters
    ----------------------
    s : shephard.Site
        Site object being converted to a string

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
    line = line + str(s.protein.unique_ID) + delimiter
    line = line + str(s.position) + delimiter
    line = line + str(s.site_type) + delimiter                
    line = line + str(s.symbol) + delimiter

    # note last required element has no trailing delimiter
    line = line + str(s.value) 
    
    if s.attributes:
        for k in s.attributes:
            atrbt = interface_tools.full_clean_string(s.attribute(k))
            line = line + delimiter + f"{k}:{atrbt}"

    line = line + "\n"

    return line
