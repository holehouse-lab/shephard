"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

from . import interface_tools 
import shephard.exceptions as shephard_exceptions
from shephard.exceptions import InterfaceException, ProteinException, SiteException

class _SitesInterface:

    def __init__(self, filename, delimiter='\t', skip_bad=True):
        """
        Expect files of the following format:

        Unique_ID, position, site type, symbol, value, key1:value1, key2:value2, ..., keyn:valuen

        Note that the first four arguments are required, while all of the key:value pairs 
        are optional. Key value must be separated by a ':', but any delimiter (other than ':') 
        is allowed

        """

        if delimiter == ':':
            raise InterfaceException('When parsing site file cannot use ":" as a delimeter because this is used to delimit key/value pairs (if provided)')


        with open(filename,'r') as fh:
            content = fh.readlines()

        ID2site={}
        
        linecount=0
        for line in content:

            linecount = linecount + 1
            sline = line.strip().split(delimiter)

            try:
                unique_ID = sline[0].strip()
                position = int(sline[1].strip())
                site_type = sline[2].strip()
                symbol = sline[3].strip()
                value = float(sline[4].strip())
                attributes = {}
            except Exception:

                msg = 'Failed parsing file [%s] on line [%i]... line printed below:\n%s'%(filename, linecount, line)

                # should update this to also display the actual error...
                if skip_bad:
                    print('Warning: %s'%msg)
                    print('Skipping this line...')
                    continue
                else:
                    raise InterfaceException('Error: %s'%msg)

            # if there's more parse attribute dictionary entries
            if len(sline) > 5:
                attributes = interface_tools.parse_key_value_pairs(sline[5:], filename, linecount, line)


            if unique_ID in ID2site:
                ID2site[unique_ID].append([position, site_type, symbol, value, attributes])
            else:
                ID2site[unique_ID] =[[position, site_type, symbol, value, attributes]]

        self.data = ID2site



##############################################
##                                          ##
##     PUBLIC FACING FUNCTIONS BELOW        ##
##                                          ##
##############################################


def add_sites_from_file(proteome, filename, delimiter='\t', safe=True, skip_bad=True, verbose=True):
    """
    Function that provides the user-facing interface for reading correctly configured SHEPHARD 
    sites files and adding those sites to the proteins of interest.
    
    A SHEPHARD sites file is a tab (or other) delineated file where each line has the following
    convention:

    Unique_ID, position, site type, symbol, value, [ key_1:value_1, key_2:value_2, ..., key_n:value_n ]
    
    Each line has six required values and then can have as many key:value pairs as may be
    desired.

    Those required values are

    unique_ID : the
    [0] = site position
    [1] = site type
    [2] = site symbol 
    [3] = site value 
    [4] = site attribute dictionary


    Parameters
    -------------
    proteome : Proteome
        Proteome object to which we're adding sites. Note that ONLY sites for which a protein
        is found will be used. Protein-Site cross-referencing is done using the protein's unique_ID
        which should be the key used in the sites_dictionary

    filename : str
        Name of the shephard site file to be read

    delimiter : str 
        String used as a delimiter on the input file. Default = '\t'

    safe : boolean 
        If set to True then any exceptions raised during the site-adding process (i.e. after file
        parsing) are acted on. If set to False, exceptions simply mean the site in question is skipped. 
        There are various reasons site addition could fail (e.g. site falls outside of protein position
        so if verbose=True then the cause of an exception is also printed to screen. It is highly 
        recommend that if you choose to use safe=False you also set verbose=True. Default = True.
        
    skip_bad : boolean
        Flag that means if bad lines (lines that trigger an exception) are encountered the code 
        will just skip them. By default this is true, which adds a certain robustness to file 
        parsing, but could also hide errors. Note that if lines are skipped a warning will be 
        printed (regardless of verbose flag). Default = True.

    verbose : boolean
        Flag that defines how 'loud' output is. Will warn about errors on adding sites.
        
    Returns
    ---------
    None
        No return value, but adds all of the passed sites to the protein
    
    """

    # check first argument is a proteome
    interface_tools.check_proteome(proteome, 'add_sites_from_file (si_sites)')

    sites_interface = _SitesInterface(filename, delimiter, skip_bad)

    # finally add the site from the dictionary generated by the SitesInterface parser
    add_sites_from_dictionary(proteome, sites_interface.data, safe, verbose)



def add_sites_from_dictionary(proteome, sites_dictionary, safe=True, verbose=False):
    """
    Function that takes a correctly formatted Sites dictionary and will add those 
    sites to the proteins in the Proteome.

    Sites dictionaries are key-value pairs, where the key is a unique_ID associated 
    with a given protein, and the value is a list of lists. Each sublist has five positions

    [0] = site position
    [1] = site type
    [2] = site symbol 
    [3] = site value 
    [4] = site attribute dictionary

    In this way, each site that maps to a give unique_ID will be added. 

    Parameters
    -------------

    proteome : Proteome
        Proteome object to which we're adding sites. Note that ONLY sites for which a protein
        is found will be used. Protein-Site cross-referencing is done using the protein's unique_ID
        which should be the key used in the sites_dictionary

    sites_dictionary : dict
        A sites dictionary is a defined dictionary that maps a unique_ID back to a list with five
        elements. Each of those elements is position-specific information for the site, specifically:

            [0] = site position
            [1] = site type
            [2] = site symbol 
            [3] = site value 
            [4] = site attribute dictionary

        Recall the only type-specific values (position and value) are cast automatically when a 
        site is added by the Protein object, so no need to do that in this function too.

        Extra elements in the each sites_dictionary value are ignored.

    safe : boolean 
        If set to True then any exceptions raised during the site-adding process are acted
        on. If set to false, exceptions simply mean the site in question is skipped. There 
        are various reasons site addition could fail (notably position of the site is  
        outside of the protein limits) and so if verbose=True then the cause of an exception 
        is also  printed to screen. It is highly recommend that if you choose to
        use safe=False you also set verbose=True
        Default = True

    verbose : boolean
        Flag that defines how 'loud' output is. Will warn about errors on adding sites.

    Returns
    ---------
    None
        No return value, but adds all of the passed sites to the protein
    
    """
    
    for protein in proteome:
        if protein.unique_ID in sites_dictionary:
            for site in sites_dictionary[protein.unique_ID]:

                try:
                    position = site[0]
                    site_type = site[1]
                    symbol = site[2]
                    value = site[3]
                    ad    = site[4]
                except:
                    raise InterfaceException('When sites dictionary for key [%s] was unable to extract five distinct parametes. Entry is:\n%s\n'% (protein.unique_ID, site))

                try:
                    protein.add_site(position, site_type, symbol, value, attributes = ad)                

                except ProteinException as e:

                    msg='- skipping site %s at %i on %s' %(site_type, position, protein)
                    if safe:
                        shephard_exceptions.print_and_raise_error(msg, e)
                    else:
                        if verbose:
                            shephard_exceptions.print_warning(msg, e)
                            continue
                    


