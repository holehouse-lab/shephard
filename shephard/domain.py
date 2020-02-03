"""
Region Class File - ShanEnt Suite 

version: 2.1b

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Pappu Lab - Washington University in St. Louis
"""

from . import sequence_utilities
from . import exceptions

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Class that defines a sequence region 
#
class Domain:
    """
    Domains are defined sub-regions within a protein. 
    
    Proteins contain a list of 0 or more domains, and each domain is associated 
    with the protein it originates from via the linking protein object.
    


    """

    ## ------------------------------------------------------------------------
    ##          
    def __init__(self, start, end, protein, domain_type, attribute_dictionary=None):

        self._start = int(start)
        self._end   = int(end)

        self._start_i0 = int(start)-1
        self._end_i0   = int(end)-1

        self._protein = protein         
        self._domain_type = domain_type

        # set attribute dictionary IF a dictionary was passed
        if isinstance(attribute_dictionary, dict):
            self._attributes = attribute_dictionary

        # set dictionary to an empty dictionary if none was passed
        elif attribute_dictionary is None:
            self._attributes = {}

        else:
            raise exceptions.DomainException('[FATAL]: If provided, protein attribute must a dictionary')

        # update unique domain types
        protein.proteome.__update_domain_types(self._domain_type)


    ## ------------------------------------------------------------------------
    ##      
    @property
    def start(self):
        return self._start

    ## ------------------------------------------------------------------------
    ##      
    @property
    def end(self):
        return self._end

    ## ------------------------------------------------------------------------
    ##      
    @property
    def protein(self):
        return self._protein

    ## ------------------------------------------------------------------------
    ##      
    @property
    def sequence(self):

        # note that  domain boundaries are inclusive, but Python's slice selection
        # and range functions are exclusive, hence the +1 offset
        return self._protein.sequence[self._start_i0:self._end_i0+1]


    ## ------------------------------------------------------------------------
    ##      
    @property
    def domain_type(self):
        return self._domain_type


    
    ## ------------------------------------------------------------------------
    ##      
    def inside_domain(self, position):
        """
        Function that returns True/False depending on if the provided position
        lies inside the domain.

        """
        return sequence_utilities.inside_region(self.start, self.end, position)

    ###############################
    ##                           ##
    ##     SITE FUNCTIONS        ##
    ##                           ##
    ###############################

    @property
    def sites(self):
        """
        Get list of sites inside the domain
        """
        
        return list(self._protein.get_sites_by_range(self.start, self.end).keys())


    def site(self, position):
        """

        """
        return self._protein._sites[int(position)]

    def get_sites_by_type(self, site_type):
        """
        Get list of sites inside the domain
        """

        return self._protein.get_sites_by_type_and_range(site_type, self.start, self.end)
        


    ###############################
    ##                           ##
    ##     TRACK FUNCTIONS       ##
    ##                           ##
    ###############################

    ## ------------------------------------------------------------------------
    ##      
    def get_track_values(self, name):
        """
        Function that returns the region of a protein's track associated with
        this domain.

        Parameters
        --------------

        """

        return self._protein.track(name).values[self._start_i0:self._end_i0]


    ## ------------------------------------------------------------------------
    ##      
    def get_track_symbols(self, name):
        """
        Function that returns the region of a protein's track symbols associated 
        with this domain.
        

        Parameters
        --------------

        """

        return self._protein.track(name).symbols[self._start_i0:self._end_i0]
        




    ## ------------------------------------------------------------------------
    ##      
    def __repr__(self):
        return "|Domain: %s, %i-%i (len=%i)| in %s" % (self._domain_type, self.start, self.end, len(self), self.protein)

    ## ------------------------------------------------------------------------
    ##      
    def __len__(self):
        return len(self.sequence)

    
            



