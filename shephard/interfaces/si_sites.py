from .interface_exceptions import InterfaceException
from . import interface_tools 
from shephard.exceptions import ProteinException

class SitesInterface:

    def __init__(self, filename, delimiter='\t'):
        """
        Expect files of the following format:

        Unique_ID, position, site type, symbol, value, key1:value1, key2:value2, ..., keyn:valuen

        Note that the first four arguments are required, while all of the key:value pairs 
        are optional. Key value must be separated by a ':', but any delimiter (other than ':') 
        is allowed

        """

        if delimiter == ':':
            raise InterfaceException('When parsing domain file cannot use ":" as a delimeter because this is used to delimit key/value pairs (if provided)')


        with open(filename,'r') as fh:
            content = fh.readlines()

        ID2site={}
        
        linecount=0
        for line in content:
            linecount=linecount+1
            sline = line.split(delimiter)

            try:
                unique_ID = sline[0].strip()
                position = int(sline[1].strip())
                site_type = sline[2].strip()
                symbol = sline[3].strip()
                value = float(sline[4].strip())
                attribute_dictionary = {}
            except Exception:

                # should update this to also display the actual error...
                raise InterfaceException('Failed parsing file [%s] on line [%i]... line printed below:\n%s'%(filename, linecount, line))

            if len(sline) > 5:
                attribute_dictionary = interface_tools.parse_key_value_pairs(sline[5:], filename, linecount, line)


            if unique_ID in ID2site:
                ID2site[unique_ID].append([position, site_type, symbol, value, attribute_dictionary])
            else:
                ID2site[unique_ID] =[[position, site_type, symbol, value, attribute_dictionary]]

        self.data = ID2site


def read_in_sites(proteome, filename, delimiter='\t', skip_bad=False):
    """
    Function that provides the user-facing interface for reading correctly configured shephard sites 
    files and adding those sites to the proteins of interest.

    A shephard sites file is a tab (or other) deliniated file where the 

    Unique_ID, position, site type, symbol, value, key1:value1, key2:value2, ..., keyn:valuen
    

    """

    # check first argument is a proteome
    interface_tools.check_proteome(proteome, 'read_in_sites (si_sites)')

    sites_interface = SitesInterface(filename, delimiter)

    for protein in proteome:
        if protein.unique_ID in sites_interface.data:
            for site in sites_interface.data[protein.unique_ID]:

                position = site[0]
                site_type = site[1]
                symbol = site[2]
                value = site[3]
                ad    = site[4]

                try:
                    protein.add_site(position, site_type, symbol, value, attributes = ad)
                except ProteinException as e:
                    if skip_bad:
                        print('Warning - skippiing site at %i on %s' %(position, protein))
                        continue
                    else:
                        raise e
                        
