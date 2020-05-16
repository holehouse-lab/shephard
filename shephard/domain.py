"""
SHEPHARD: 
Sequence-based Hierachical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

from . import sequence_utilities
from . import exceptions
from shephard.tools import domain_tools

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Class that defines a sequence region 
#
class Domain:
    """
    Domains are defined sub-regions within a protein. 
    
    Proteins contain a list of 0 or more domains, and each domain is associated 
    with the protein it originates from via the linking protein object.

    Domains are indexed using the same indexing as the overall protein sequence (i.e.
    a protein does not automatically start from 1), and as such the 'native' sequence
    indexing should be used for working with Proteins. This is a long-winded way of saying
    position X refers to the same residue regardless of if it's taken from the Protein or Domain
    or Track object.

    Parameters
    -------------

    start : int
        Start position in sequence (recall we index from 1)

    end : int
        End position in sequence (recall we index from 1)

    protein : Protein
        Protein object for which this Domain is part of

    domain_type :  str
        Name of the domain type - can be any free-form string

    attribute_dictionary : dict
        Dictionary where key/value pairs allow a Domain to have
        arbitrary metadata associated with it.
        
    """

    ## ------------------------------------------------------------------------
    ##          
    def __init__(self, start, end, protein, domain_type, attribute_dictionary=None):
        """ 
        """
        
        if start > end:
            raise DomainException("Trying to a domain to protein [%s] where start site is bigger than the end site (positions: %i-%i - this does not work!" %(str(protein), start, end))

        # check the domain falls within the region
        helper_string="Trying to add domain to protein [%s] at positions [%i-%i] - this falls outside the protein's dimensions [%i-%i]" %(protein, start, end, 1, protein._len)
        protein._check_position_is_valid(start, helper_string)
        protein._check_position_is_valid(end, helper_string) 

        # assign if all OK
        self._start = int(start)
        self._end   = int(end)

        self._protein = protein         
        self._domain_type = domain_type

        # set attribute dictionary IF a dictionary was passed. Otherwise we just ignore
        # anything passed to attribute_dictionary
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
            Name that will be used to identify the attribute

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
        """
        **[Property]**: Returns the start position that defines this domain
        """

        return self._start



    ## ------------------------------------------------------------------------
    ##      
    @property
    def end(self):
        """
        **[Property]**: Returns the end position that defines this domain
        """
        return self._end



    ## ------------------------------------------------------------------------
    ##      
    @property
    def protein(self):
        """
        **[Property]**: Returns the Protein that this Domain is associated with
        """
        return self._protein



    ## ------------------------------------------------------------------------
    ##      
    @property
    def sequence(self):
        """
        **[Property]**: Returns the amino acid sequence associated with this domain
        """
        return self._protein.get_sequence_region(self._start, self._end)


    ## ------------------------------------------------------------------------
    ##      
    @property
    def domain_type(self):
        """
        Returns the domain type as a string
        """
        return self._domain_type




    ######################################
    ##                                  ##
    ##     DOMAIN  FUNCTIONS             #
    ##                                  ##
    ######################################

    
    ## ------------------------------------------------------------------------
    ##      
    def inside_domain(self, position):
        """
        Function that returns True/False depending on if the provided position
        lies inside the domain.

        Parameters
        ------------
        position : int
            Position in the sequence

        Returns
        -----------
        bool
            Returns True if position is inside the domain region, else False
        

        """
        return sequence_utilities.inside_region(self.start, self.end, position)


    ## ------------------------------------------------------------------------
    ##      
    def domain_overlap(self, domain2):
        """
        Function that takes in a second domain and calculates if those two domains
        overlap at all. This is a binary check and does not compute the extent of 
        overlap.


        Parameters
        ------------
        domain2 : Domain
            The Domain object of interest

        Returns
        -----------
        bool
            Returns True if the domains overlap, False if not. Note this will throw
            an exception if the domains are from different proteins.
        

        """
        return domain_tools.domain_overlap(self, domain2)


    ######################################
    ##                                  ##
    ##     DOMAIN SITE FUNCTIONS        ##
    ##                                  ##
    ######################################

    ## ------------------------------------------------------------------------
    ##
    @property
    def sites(self):
        """
        Get list of all sites inside the domain
        
        Returns
        --------
        list
            Returns a list of positions for which sites exist in this domain
        """
        
        return list(self._protein.get_sites_by_range(self.start, self.end).keys())

    ## ------------------------------------------------------------------------
    ##
    def site(self, position):
        """
        Returns the list of sites that are found at a given position. Note that - in general
        site() should be used to retrieve sites you know exist while get_sites_by_position()
        offers a way to more safely get sites at a position. Site will throw an exception 
        if the position passed does not exist (while get_sites_by_position() will not).

        Parameters
        -------------
        position : int
            Defines the position in the sequence we want to interrogate 

        Returns
        ---------
        list
            Returns a list with between 1 and n sites. Will raise an exception if 
            the passed position cannot be found in the codebase.

        """

        ipos = int(position)
        if sequence_utilities.inside_region(domain.start, domain.end, ipos):
            return self._protein._sites[int(position)]
        else:
            raise DomainException('Passed position [%i] is outside of the domain boundaries [%i-%i]' %(ipos, domain.start, domain.end))


    ## ------------------------------------------------------------------------
    ##
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
        


    #######################################
    ##                                   ##
    ##      DOMAIN TRACK FUNCTIONS       ##
    ##                                   ##
    #######################################

    ## ------------------------------------------------------------------------
    ##      
    def get_track_values(self, name, safe=True):
        """
        Function that returns the region of a protein's values- track associated with
        this domain.
        
        If the track name is not found in this protein and safe is True, this will throw 
        an exception, otherwise (if safe=False) then if the track is missing the function        
        will return None.

        Parameters
        --------------
        
        name : str
            Track name

        safe : boolean
            If set to True, missing tracks trigger an exception, else they 
            just return None

        Returns
        ----------
        list
            Returns a list of floats that corresponds to the set of residues associated
            with the domain of interest, or None if the track does not exist and safe=False.

        """
        
        t = self._protein.track(name, safe)

        if t is not None:        
            return t.values_region(self._start, self._end)
        else:
            return None


    ## ------------------------------------------------------------------------
    ##      
    def get_track_symbols(self, name, safe=True):
        """
        Function that returns the region of a protein's symbols track associated with
        this domain.
        
        If the track name is missing and safe is True, this will throw an exception,
        otherwise (if safe=False) then if the track is missing the function
        returns None

        Parameters
        --------------
        
        name : str
            Track name

        safe : boolean
            If set to True, missing tracks trigger an exception, else they 
            just return None

        Returns
        ----------
        list
            Returns a list of strs that corresponds to the set of residues associated
            with the domain of interest.

        """

        t = self._protein.track(name, safe)
        
        if t is not None:            
            return t.symbols_region(self._start, self._end)
        else:
            return None
        

    ## ------------------------------------------------------------------------
    ##      
    def __repr__(self):
        return "|Domain: %s, %i-%i (len=%i)| in %s" % (self._domain_type, self.start, self.end, len(self), self.protein)

    ## ------------------------------------------------------------------------
    ##      
    def __len__(self):
        return len(self.sequence)

    
            



