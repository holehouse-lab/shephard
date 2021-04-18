"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

from . import general_utilities
from .exceptions import ProteinException, SiteException


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Class that defines a Site in sequence
#
class Site:
    """
    Sites are defined sub-positions within a protein that map to a specific residue.
    
    Proteins contain a list of 0 or more sites, and each site is associated 
    with the protein it originates from via the linking protein object.


    Parameters
    -------------

    position : int
        Position in sequence associated with this site

    site_type :  str
        Identifier for the site type

    protein : Protein
        Protein object for which this site is part of

    symbol : str
        Symbol associated with the site - a string-based representation of something specific
        to the site. For a mutation that could be the residue the native residue mutates to, 
        for example. Is not required. Default = None.

    value : float
        Value associated with the site - a numerical value (cast to a float). Is not required.
        Default = None.

    attributes : dict (optional)
        The attributes dictionary provides a key-value pairing for arbitrary information.
        This could include different types of identifies, track generator functions,
        a set of Site partners, or anything else one might wish to associated with the
        track as a whole. Default = None
    
    """
    
    def __init__(self, position, site_type, protein, symbol=None, value=None, attributes=None):

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
        general_utilities.variable_is_dictionary(attributes, SiteException, 'attributes argument passed to site %i in protein %s is not a dictionary' %(self._position, self._protein), or_none=True)
        
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
        return self._protein.get_residue(self._position)


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
        **[Property]**: Provides a list of the keys associated with every attribute associated
        with this Site.

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

        Recall that attributes are name : value pairs, where the 'value' can be 
        anything and is user defined. This function will return the value associated 
        with a given name.

        Parameters
        ----------------
        name : str
             The attribute name. A list of valid names can be found by calling the
             ``<Site>.attributes()`` (which returns a list of the valid names)

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
                raise SiteException('Requesting attribute [%s] from Site [%s] but this attribute has not been assigned' % (name, str(self))) 

            # if safe not passed just return None
            else:
                return None
                


    ## ------------------------------------------------------------------------
    ##
    def add_attribute(self, name, val, safe=True):
        """
        Function that adds an attribute. Note that if safe is true, this function will
        raise an exception if the attribute is already present. If safe=False, then
        an existing value will be overwritten.

        Parameters
        ----------------

        name : str
            The parameter name that will be used to identify it

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
                raise SiteException("Trying to add attribute [%s=%s] to Site [%s] but this attribute is already set.\nPossible options are: %s" %(name,val, str(self), str(self._attributes.keys())))
                
        self._attributes[name] = val

        
        

    ## ------------------------------------------------------------------------
    ##
    def get_local_sequence_context(self, offset = 5):
        """
        Returns the local amino acid context around a residue +/- the offset provided.
        
        Note that the offset extends to the start/end of sequence and then silently
        truncates.

        Parameters
        -----------
        offset : int
            Defines the +/- region around the position which is used to define
            the local sequence context.

        Returns
        --------
        str
            Returns an amino acid sequence that corresponds to the local sequence
            context around the site of interest

        """

        return self._protein.get_local_sequence_context(self.position, offset)    




    #######################################
    ##                                   ##
    ##        SITE TRACK FUNCTIONS       ##
    ##                                   ##
    #######################################

    ## ------------------------------------------------------------------------
    ##      
    def get_track_values(self, name, offset=5, safe=True):
        """
        Function that returns the region of a protein's values- track associated with
        this site, +/- some offset.
        
        If the track name is missing and safe is True, this will throw an exception,
        otherwise (if safe=False) then if the track is missing the function
        returns None

        Parameters
        --------------
        
        name : str
            Track name

        offset : int
            +/- values around the site from which regions are taken

        safe : bool
            If set to True, missing tracks trigger an exception, else they 
            just return None



        Returns
        ----------
        list
            Returns a list of floats that corresponds to the set of residues associated
            with the domain of interest

        """
        
        (p1, p2) = sequence_utilities.get_bounding_sites(self._position, offset, self._protein._len)
        

        # because calling values_region only makes sense IF the track exists, we have to split
        # this into two operations
        t = self._protein.track(name, safe)

        if t is not None:
            return t.values_region(p1, p2)
        else:
            return None


    ## ------------------------------------------------------------------------
    ##      
    def get_track_symbols(self, name, offset=5, safe=True):
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

        safe : bool
            If set to True, missing tracks trigger an exception, else they 
            just return None

        Returns
        ----------
        list
            Returns a list of strs that corresponds to the set of residues associated
            with the domain of interest.

        """

        (p1, p2) = sequence_utilities.get_bounding_sites(self._position, offset, self._protein._len)

        # because calling symbols_region only makes sense IF the track exists, we have to split
        # this into two operations
        t = self._protein.track(name, safe)

        if t is not None:
            return t.symbols_region(p1,p2)
        else:
            return None



    ## ------------------------------------------------------------------------
    ##
    def __repr__(self):             
        return "|Site: %s @ %i in protein %s" % (self._site_type, self.position, self.protein.unique_ID)
    
    
