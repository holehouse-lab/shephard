from .interface_exceptions import InterfaceException
from . import interface_tools 
from shephard.exceptions import ProteinException

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
                attribute_dictionary = {}
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
                attribute_dictionary = interface_tools.parse_key_value_pairs(sline[5:], filename, linecount, line)


            if unique_ID in ID2site:
                ID2site[unique_ID].append([position, site_type, symbol, value, attribute_dictionary])
            else:
                ID2site[unique_ID] =[[position, site_type, symbol, value, attribute_dictionary]]

        self.data = ID2site



##############################################
##                                          ##
##     PUBLIC FACING FUNCTIONS BELOW        ##
##                                          ##
##############################################


def add_sites_from_dictionary(proteome, sites_dictionary):
    """
    Function that takes a correctly formatted Sites dictionary and will add those 
    sites to the proteins in the proteome.

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

        Recall the only type-sepecific values (position and value) are cast automatically when a 
        site is added by the Protein object, so no need to do that in this function too.

        Exta elements in the each sites_dictionary value are ignored.

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

                protein.add_site(position, site_type, symbol, value, attributes = ad)                



def add_sites_from_file(proteome, filename, delimiter='\t', skip_bad=True):
    """
    Function that provides the user-facing interface for reading correctly configured SHEPHARD 
    sites files and adding those sites to the proteins of interest.
    
    A SHEPHARD sites file is a tab (or other) deliniated file where each line has the following
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


    TBC
    

    """

    # check first argument is a proteome
    interface_tools.check_proteome(proteome, 'add_sites_from_file (si_sites)')

    sites_interface = _SitesInterface(filename, delimiter, skip_bad)

    # finally add the domains from the dictionary generated by the DomainsInterface parser
    add_sites_from_dictionary(proteome, sites_interface.data)


                        
