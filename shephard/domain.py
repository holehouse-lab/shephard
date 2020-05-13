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

        self._protein = protein         
        self._domain_type = domain_type
        
        if start > end:
            raise DomainException("Trying to a domain to protein [%s] where start site is bigger than the end site (positions: %i-%i - this does not work!" %(str(protein), start, end))

        # check the domain falls within the region
        helper_string="Trying to add domain to protein [%s] at positions [%i-%i] - this falls outside the protein's dimensions [%i-%i]" %(protein, start, end, 1, protein._len)
        protein._check_position_is_valid(start, helper_string)
        protein._check_position_is_valid(end, helper_string) 

        # set attribute dictionary IF a dictionary was passed. Otherwise we just ignore
        # anything bassed to attribute_dictionary
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
    def attributes(self):
        """
        Provides a list of the keys associated with every attribute associated
        with this domain.

        Returns
        -------
        list
            returns a list of the attribute keys associated with the domain. 


        """
        return list(self._attributes.keys())
    
    
    ## ------------------------------------------------------------------------
    ##
    def attribute(self, name, safe=True):

        """
        Function that returns a specific attribute as defined by the name. 

        Recall that attributes are name : value pairs, where the 'value' can be 
        anything and is user defined. This function will return the value associated 
        with a given name.

        Parameters
        ----------------
        name : str
             The attribute name. A list of valid names can be found by calling the
             ``<domain>.attributes()`` (which returns a list of the valid names)

        safe : bool (default = True)
            Flag which if true with throw an exception if an attribute with the same
            name  already exists
            
        Returns
        ---------
        Unknown 
            Will either return whatever was associated with that attribute (which could be anything)
            or None if that attribute is missing.
        
        """

        # if name is in the _atributes dictionary the  return
        if name in self._attributes:
            return self._attributes[name]
        else:

            # else if safe was passed raise an exception if that attribute was missing
            if safe:
                raise exceptions.DomainException('Requesting attribute [%s] from domain [%s] but this attribute has not been assigned' % (name, str(self))) 

            # if safe not passed just return None
            else:
                return None
                


    ## ------------------------------------------------------------------------
    ##
    def add_attribute(self, name, val, safe=True):
        """
        Function that adds an attribute. Note that if safe is true, this function will
        raise an exception if the attribute is already present. If safe=False, then
        an exisiting value will be overwritten.

        Parameters
        ----------------

        name : str
            The parameter name that will be used to identfy it

        val : <anything>
            An object or primitive we wish to associate with this attribute

        safe : bool (default = True)
            Flag which if True with throw an exception if an attribute with the same
            name already exists, otherwise the newly introduced attribute will overwrite
            the previous one.

        Returns
        ---------
            None - but adds an attribute to the calling object

        """

        if safe:
            if name in self._attributes:
                raise exceptions.DomainException("Trying to add attribute [%s=%s] to domain [%s] but this attribute is already set.\nPossible options are: %s" %(name,val, str(self), str(self._attributes.keys())))
                
        self._attributes[name] = val 

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

        # note that  domain boundaries are inclusive, and the protein.sequence
        # variable is a Protseq object which corrects the indexing such that
        # slicing natively works inclusively and with the approriate indexing
        # for real-world position
        return self._protein.sequence[self._start:self._end]


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


    ######################################
    ##                                  ##
    ##     DOMAIN SITE FUNCTIONS        ##
    ##                                  ##
    ######################################

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

        Parameters
        ------------

        site_type : string
            The site type identifier for which the function will search for matching sites
        
        Returns
        --------

        list
            Returns a dictionary, where each key-value pair is:

                key - site position (integer)
                value - list of one or more site object
        
        """

        return self._protein.get_sites_by_type_and_range(site_type, self.start, self.end)
        


    #####################################
    ##                                 ##
    ##     DOMAI TRACK FUNCTIONS       ##
    ##                                 ##
    #####################################

    ## ------------------------------------------------------------------------
    ##      
    def get_track_values(self, name, safe=True):
        """
        Function that returns the region of a protein's track associated with
        this domain.

        Parameters
        --------------

        """
        

        return self._protein.track(name)._values[self._start:self._end]


    ## ------------------------------------------------------------------------
    ##      
    def get_track_symbols(self, name):
        """
        Function that returns the region of a protein's track symbols associated 
        with this domain.
        

        Parameters
        --------------

        """

        return self._protein.track(name).symbols_slice[self._start:self._end]
        




    ## ------------------------------------------------------------------------
    ##      
    def __repr__(self):
        return "|Domain: %s, %i-%i (len=%i)| in %s" % (self._domain_type, self.start, self.end, len(self), self.protein)

    ## ------------------------------------------------------------------------
    ##      
    def __len__(self):
        return len(self.sequence)

    
            



