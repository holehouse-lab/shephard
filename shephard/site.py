"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

from . import general_utilities
from . import sequence_utilities
from .exceptions import ProteinException, SiteException


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Class that defines a Site in sequence
#
class Site:
    
    def __init__(self, position, site_type, protein, symbol=None, value=None, attributes=None):

        """
        Sites are defined sub-positions within a protein that map to a 
        specific residue.
    
        Proteins contain a list of 0 or more sites, and each site is 
        associated with the protein it originates from via the linking 
        protein object.


        Parameters
        -------------

        position : int
            Position in sequence associated with this site

        site_type :  str
            Identifier for the site type

        protein : Protein
            Protein object for which this site is part of
        
        symbol : str (default = None)
            Symbol associated with the site - a string-based representation 
            of something specific to the site. For a mutation that could be 
            the residue the native residue mutates to, for example. 

        value : float (default = None)
            Value associated with the site - a numerical value (cast to a 
            float). 
        
        attributes : dict (default = None)
            The attributes dictionary provides a key-value pairing for 
            arbitrary information. This could include different types of 
            identifies, track generator functions, a set of Site partners, 
            or anything else one might wish to associated with the track as 
            a whole. 
        
        """


        # absolute position in protein associated with the site
        self._position  = int(position)
        
        # reference to the protein object from which the site was taken
        self._protein   = protein

        # the site type identifier
        self._site_type = str(site_type)

        # a numerical value associated with the site
        self._value     = general_utilities.cast_or_none(value, float)

        # a symbol associated with the site
        self._symbol    = general_utilities.cast_or_none(symbol, str)

        # verify that the attributes dictionary is a dictionary
        general_utilities.variable_is_dictionary(attributes, SiteException, 'attributes argument passed is not a dictionary', or_none=True)
        
        if attributes is None:
            self._attributes = {}
        else:
            self._attributes = attributes

        # update the proteome if this is a novel type of site
        protein.proteome.__update_site_types(self._site_type)


    ## ------------------------------------------------------------------------
    ##
    @property
    def residue(self):
        """
        Returns the amino acid residue associated with the site position as
        a string.
        """
        return self._protein.residue(self._position)


    ## ------------------------------------------------------------------------
    ##
    @property
    def position(self):
        """
        Returns the actual sequence indexed position as an int (recall protein
        indexing starts at 1).
        """
        return self._position


    ## ------------------------------------------------------------------------
    ##
    @property
    def protein(self):
        """
        Return the Protein object this site is found within
        """
        return self._protein


    ## ------------------------------------------------------------------------
    ##
    @property
    def site_type(self):
        """
        Returns the site type (string)
        """
        return self._site_type


    ## ------------------------------------------------------------------------
    ##
    @property
    def symbol(self):
        """
        Returns the symbol associated with this site. Note a symbol is either
        None or a str type.
        """
        return self._symbol


    ## ------------------------------------------------------------------------
    ##
    @property
    def value(self):
        """
        Returns the value associated with this site. Note a value is either
        None or a float type.
        """
        return self._value

    ## ------------------------------------------------------------------------
    ##
    def update_site_value(self, new_value):
        """
        Function that updates the site value. Values must be
        numerical (float) or None.
        
        Parameters
        -----------
        new_value : float (or None)
            Updated value

        Returns
        -----------
        None
            Nothing but sets the value to be the new value

        """
        if new_value is None:
            self._value = None
        else:
            self._value = float(new_value)

    ## ------------------------------------------------------------------------
    ##
    def update_site_symbol(self, new_symbol):
        """
        Function that updates the site_symbol. The site tyoe must be
        a string. Note this function also updates the proteome list of
        non-redudant sites
        
        Parameters
        -----------
        new_symbol : str (or None)
            Updated symbol

        Returns
        -----------
        None
            Nothing but sets the symbol to be the new symbol

        """
        if new_symbol is None:
            self._symbol = None
        else:
            self._symbol = str(new_symbol)





     ###################################
    ##                               ##
    ##     ATTRIBUTE FUNCTIONS       ##
    ##                               ##
    ###################################


    ## ------------------------------------------------------------------------
    ##
    @property
    def attributes(self):
        """
        **[Property]**: Provides a list of the keys associated with every 
        attribute associated with this Site.        

        Returns
        -------
        list
            returns a list of the attribute keys associated with the protein. 
        """

        return list(self._attributes.keys())



    ## ------------------------------------------------------------------------
    ##
    def attribute(self, name, safe=True):

        """
        Function that returns a specific attribute as defined by the name. 

        Recall that attributes are name : value pairs, where the 'value' 
        can be anything and is user defined. This function will return 
        the value associated with a given name.
        
        

        Parameters
        ----------------
        name : str
             The attribute name. A list of valid names can be found by 
             calling the ``<Site>.attributes()`` (which returns a list 
             of the valid names)
             
        safe : bool (default = True)
            Flag which if true with throw an exception if an attribute 
            with the same name already exists.
                        
        Returns
        ---------
        Unknown 
            Will either return whatever was associated with that attribute 
            (which could be anything) or None if that attribute is missing.
                    
        """

        # if name is in the _atributes dictionary the  return
        if name in self._attributes:
            return self._attributes[name]
        else:

            # else if safe was passed raise an exception if that attribute was missing
            if safe:
                raise SiteException('Requesting attribute [%s] from Site [%s] but this attribute has not been assigned' % (name, str(self))) 

            # if safe not passed just return None
            else:
                return None
                


    ## ------------------------------------------------------------------------
    ##
    def add_attribute(self, name, val, safe=True):
        """
        Function that adds an attribute. Note that if safe is true, this 
        function will raise an exception if the attribute is already present. 
        If safe=False, then an existing value will be overwritten.

        Parameters
        ----------------

        name : str
            The parameter name that will be used to identify it

        val : <anything>
            An object or primitive we wish to associate with this attribute

        safe : bool (default = True)
            Flag which if True with throw an exception if an attribute with 
            the same name already exists, otherwise the newly introduced 
            attribute will overwrite the previous one.

        Returns
        ---------
            None - but adds an attribute to the calling object

        """

        if safe:
            if name in self._attributes:
                raise SiteException("Trying to add attribute [%s=%s] to Site [%s] but this attribute is already set.\nPossible options are: %s" %(name,val, str(self), str(self._attributes.keys())))
                
        self._attributes[name] = val


    ## ------------------------------------------------------------------------
    ##
    def remove_attribute(self, name, safe=True):
        """
        Function that removes a given attribute from the Site based on the 
        passed attribute name. If the passed attribute does not exist or is not 
        associate with the Site then this will trigger an exception 
        unless safe=False.

        Parameters
        ----------------

        name : str
            The attribute name that will be used to identify it

        safe : bool (default = True)
            Flag which if True with throw an exception if an 
            attribute this name does not exists. If set to
            False then if an attribute is not found it is simply
            ignored
            
        Returns
        ---------
        None
            No return type but will remove an attribute from the 
            protein if present.
            
        """

        if name not in self._attributes:
            if safe:
                raise ProteinException(f'Passed attribute [{name}] not found in {self}')
        else:
            del self._attributes[name]
        
        

    ## ------------------------------------------------------------------------
    ##
    def get_local_sequence_context(self, offset = 5):
        """
        Returns the local amino acid context around a residue +/- the offset
        provided.
        
        Note that the offset extends to the start/end of sequence and then 
        silently truncates.
        
        Parameters
        -----------
        offset : int
            Defines the +/- region around the position which is used to 
            define the local sequence context.
            

        Returns
        --------
        str
            Returns an amino acid sequence that corresponds to the local 
            sequence context around the site of interest
        """

        return self._protein.get_sequence_context(self.position, offset)    



    #######################################
    ##                                   ##
    ##        SITE domain FUNCTIONS      ##
    ##                                   ##
    #######################################

    ## ------------------------------------------------------------------------
    ##      
    def get_domains(self, offset=0, safe=True):
        """
        Function that returns the set of domains that the site lies within. The
        oofset parameter defines the wiggle room +/- that is tolerated, but
        defaults to 0.         

        Parameters
        --------------

        offset : int (default = 0)
            +/- values around the site from which regions are taken.

        safe : bool (default = True)
            If set to True, missing tracks trigger an exception, else they 
            just return None.

        Returns
        ----------
        list
            Returns a list of domain objects for which this site can be found
            in or near

        """
        
        valid_domains = []

        for domain in self.protein.domains:

            # get start and end positions with offset
            (p1, p2) = sequence_utilities.get_bounding_sites(domain.start, offset, self._protein._len)
            start = p1
            
            (p1, p2) = sequence_utilities.get_bounding_sites(domain.end, offset, self._protein._len)
            end = p2

            if start <= self.position and end >= self.position:
                valid_domains.append(domain)

        return valid_domains



    #######################################
    ##                                   ##
    ##        SITE TRACK FUNCTIONS       ##
    ##                                   ##
    #######################################

    ## ------------------------------------------------------------------------
    ##      
    def get_track_values(self, name, offset=0, safe=True):
        """
        Function that returns the region of a protein's values- track 
        associated with this site, +/- some offset.
        
        If the track name is missing and safe is True, this will throw 
        an exception, otherwise (if safe=False) then if the track is 
        missing the function returns None.


        Parameters
        --------------
        
        name : str
            Track name

        offset : int (default = 0)
            +/- values around the site from which regions are taken

        safe : bool (default = True)
            If set to True, missing tracks trigger an exception, else 
            they just return None


        Returns
        ----------
        list
            Returns a list of floats that corresponds to the set of 
            residues associated with the domain of interest
            

        """
        
        (p1, p2) = sequence_utilities.get_bounding_sites(self._position, offset, self._protein._len)
        

        # because calling values_region only makes sense IF the track exists, we have to split
        # this into two operations. Note this throws an exception if safe=True and the track
        # does not exist
        t = self._protein.track(name, safe)

        if t is not None:
            return t.values_region(p1, p2)
        else:
            return None


    ## ------------------------------------------------------------------------
    ##      
    def get_track_value(self, name, safe=True):
        """
        Function that returns the value associated with the track at
        the residue position associated with this site.

        If the track name is missing and safe is True, this will throw 
        an exception, otherwise (if safe=False) then if the track is 
        missing the function returns None.


        Parameters
        --------------
        
        name : str
            Track name

        safe : bool (default = True)
            If set to True, missing tracks trigger an exception, else 
            they just return None

        Returns
        ----------
        float or int
            Returns the value associated with the track of interest
            at this site.

        """
        return self.get_track_values(name, offset=0, safe=safe)[0]


    ## ------------------------------------------------------------------------
    ##      
    def get_track_symbols(self, name, offset=0, safe=True):
        """
        Function that returns the region of a protein's symbols track 
        associated with this site.
                
        If a Track of this name is not associated with the underlying
        protein and safe is True, this will throw an exception, 
        otherwise (if safe=False) then if the track is missing the 
        function returns None.

        Parameters
        --------------
        
        name : str
            Track name

        offset : int (default = 0)
            +/- values around the site from which regions are taken

        safe : bool (default = True)
            If set to True, missing tracks trigger an exception, 
            else they just return None


        Returns
        ----------
        list
            Returns a list of strs that corresponds to the set of 
            residues associated with the domain of interest.
        """

        (p1, p2) = sequence_utilities.get_bounding_sites(self._position, offset, self._protein._len)

        # because calling symbols_region only makes sense IF the track exists, we have to split
        # this into two operations. Note this throws an exception if safe=True and the
        # track does not exist
        t = self._protein.track(name, safe)

        if t is not None:
            return t.symbols_region(p1,p2)
        else:
            return None


    ## ------------------------------------------------------------------------
    ##      
    def get_track_symbol(self, name, safe=True):
        """
        Function that returns the symbol associated with the track at
        the residue position associated with this site.
                
        If a Track of this name is not associated with the underlying
        protein and safe is True, this will throw an exception, 
        otherwise (if safe=False) then if the track is missing the 
        function returns None.

        Parameters
        --------------
        
        name : str
            Track name

        safe : bool (default = True)
            If set to True, missing tracks trigger an exception, 
            else they just return None

        Returns
        ----------
        str
            Returns the string associated with the symbol at this site

        """

        return self.get_track_symbols(name, offset=0, safe=safe)[0]


    ## ------------------------------------------------------------------------
    ##
    def __repr__(self):             
        return "|Site: %s @ %i in protein %s" % (self._site_type, self.position, self.protein.unique_ID)
    
    
