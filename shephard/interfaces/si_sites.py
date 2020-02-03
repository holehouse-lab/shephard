from .interface_exceptions import InterfaceException
from . import interface_tools 

class SitesInterface:

    def __init__(self, filename, delimiter='\t'):
        """
        Expect files of the followin format:

        Unique_ID, position, site type, symbol, value


        """

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
            except Exception:

                # should update this to also display the actual error...
                raise InterfaceException('Failed parsing file [%s] on line [%i]... line printed below:\n%s'%(filename, linecount, line))


            if unique_ID in ID2site:
                ID2site[unique_ID].append([position, site_type, symbol, value])
            else:
                ID2site[unique_ID] =[[position, site_type, symbol, value]]

        self.data = ID2site


def read_in_sites(proteome, filename, delimiter='\t'):

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
                
                protein.add_site(position, site_type, symbol, value)
