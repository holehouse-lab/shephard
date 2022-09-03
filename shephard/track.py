"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

from . exceptions import TrackException

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Class that defines a Site in sequence
#
class Track:
    
    def __init__(self, name, protein, values=None, symbols=None, attribute_dictionary = None):

        """
        Tracks define information that maps along a protein sequence.

        A Track is, fundamentally, a vector which is the length of the 
        sequence. This could be an way to re-code the amino acid sequence, 
        or reflect some kind of sliding window analysis.
        
        Tracks can can either define a set of symbols that convert residues 
        to symbols (i.e. discrete classifications) or values (i.e. floating 
        number values associated with each position). Note that Tracks 
        CANNOT define both symbolic and numerical data (i.e. must be one or 
        the other).
    
        Parameters
        ------------
        
        name : string
            Defines the name of the track. This can be any value, but should
            be something that makes sense. The name can be used by analysis 
            routines.
                    
        protein : Protein object
            the protein from which the track is being added to
        
        values : iterable of numerical values (default = None)
            This iterable is passed over and convert into a list of floats. 
            Must be same length as the number of residues in the protein.
        

        symbols : iterable of strings (default = None)
            This iterable is directly assigned to the track.symbols         
            variable. Must be same length as the number of residues in 
            the protein.      
                
        attribute_dictionary : dict (default = None)
            The attribute_dictionary provides a key-value pairing for 
            arbitrary information. This could include different types of 
            identifies, track generator functions, a set of Track partners, 
            or anything else one might wish to associated with the track 
            as a whole.
        
        """

        # if values were provided for the track...
        if values is not None:

            # cannot have symbols and values!
            if symbols is not None:
                raise TrackException(f'Added tracks must be include either symbols or values but not both [Track={name}, Protein={protein}')
                
            # check the values provided is the same length as the number of residues - if not raise an exception
            if len(protein.sequence) != len(values):
                raise TrackException(f'Track length of %i does not match protein length of %i (values track)\b Track = %s\nProtein=%s' %(len(values), len(protein.sequence), name, str(protein)))

            # convert values to list of floats
            try:                
                values = [float(i) for i in values]
                
                # add leading zero for index purposes
                values = [0.0] + values
            except ValueError:
                raise TrackException(f'Unable to convert values passed into float64 numpy array [Track={name}, Protein={protein}')
            track_type = 'values'


        # if the symbols were provided
        elif symbols is not None:

            # check lengths match up
            if len(protein.sequence) != len(symbols):
                raise TrackException(f'Track length (symbols) does not match protein length\nTrack: {name}, length={len(symbols)}\nProtein: {protein}')

            # if we passed a list (which we now know is the right length) we're good!
            if isinstance(symbols, list):
                pass

            # if we passed a string convert to a list
            elif(symbols, str):
                symbols = list(symbols)
            
            else:
                raise TrackException(f'Unable to convert passed symbols track to a list of symbols. Symbols track should be either a list of symbols or a string. [Track={name}, Protein={protein}')

            # add leading ('-') for index purposes
            symbols = ['-'] + symbols 
            track_type = 'symbols'


        # if NEITHER symbols nor track were provided through an exception
        if symbols is None and values is None:
            raise TrackException('Empty track provided [Track=%s, Protein=%s' %(name, str(protein)))

            
        # set attribute dictionary IF a dictionary was passed
        if isinstance(attribute_dictionary, dict):
            self._attributes = attribute_dictionary

        # set dictionary to an empty dictionary if none was passed
        elif attribute_dictionary is None:
            self._attributes = {}

        else:
            raise exceptions.TrackException('[FATAL]: If provided, Track attribute must a dictionary')

        self._values  = values
        self._symbols = symbols
        self._name = name
        self._protein = protein
        self._track_type = track_type
        
        # update track name types
        protein.proteome.__update_track_names(self._name, self._track_type)


    ## ------------------------------------------------------------------------
    ##
    @property
    def name(self):
        """
        **[Property]**: Returns the track name
        """
        return self._name


    ## ------------------------------------------------------------------------
    ##
    @property
    def track_type(self):
        """
        **[Property]**: Returns the track type. Will always be one of 'values'
        or 'symbols'.
        """
        return self._track_type

    ## ------------------------------------------------------------------------
    ##
    @property
    def values(self):
        """
        **[Property]**: Returns a list that matches the complete set of values
        for this track. If no values are assigned returns None
        """
        if self._values is None:
            return None
        else:
            return self._values[1:]

    ## ------------------------------------------------------------------------
    ##
    @property
    def symbols(self):
        """
        **[Property]**: Returns a list that matches the complete set of symbols
        for this track. If no symbols are assigned returns None
        """
        if self._symbols is None:
            return None
        else:
            return self._symbols[1:]

    ## ------------------------------------------------------------------------
    ##
    @property
    def protein(self):
        """
        **[Property]**: Returns the Protein that this Track is associated with
        """
        return self._protein


    ## ------------------------------------------------------------------------
    ##
    def values_region(self, start, end=None):
        """
        Returns a single value or a subregion of values, depending on if
        a start and end position are provided or just a start position

        Parameters
        ----------
        start : int
            Starting position of interest

        end : int
            Ending position of interest. If not provided
            then only the 

        Returns
        --------
        list
            Returns a list of values that maps to the residues in the intervening region
            defined by start and end)

        """

        # this list comprehension checks start and end are valid options
        if end is not None:
            [self._protein._check_position_is_valid(i, helper_string=f'Invalid position [{start}] passed to track {str(self)}') for i in [start,end]]
            return self._values[start:end+1]
        else:
            self._protein._check_position_is_valid(start, helper_string = f'Invalid position [{start}] passed to track {str(self)}')
            return self._values[start:start+1][0]
            
        
        

    ## ------------------------------------------------------------------------
    ##
    def symbols_region(self, start, end=None):
        """
        Returns a single symbol or a subregion of symbols, depending on if
        a start and end position are provided or just a start position.

        Parameters
        ----------
        start : int
            Starting position of interest.

        end : int
            Ending position of interest.

        Returns
        --------
        list
            Returns a list of values that maps to the residues in the 
            intervening region defined by start and end).

        """


        # this list comprehension checks start and end are valid options
        if end is not None:
            [self._protein._check_position_is_valid(i, helper_string=f'Invalid position [{start}] passed to track {str(self)}') for i in [start,end]]
            return self._symbols[start:end+1]
        else:
            self._protein._check_position_is_valid(start, helper_string = f'Invalid position [{start}] passed to track {str(self)}')
            return self._symbols[start:start+1][0]



    ## ------------------------------------------------------------------------
    ##
    def value(self, position, safe=True):
        """
        Returns a single value from the passed position.

        Parameters
        ----------
        position : int
            Starting position of interest.

        safe : bool (default = True)
            Flag which if true with throw an exception if a
            value is requested from a symbol track.

        Returns
        --------
        float
            Returns a value associated with the passed position

        """
        if self.values is None:
            if safe is True:
                raise TrackException('Requesting value from a symbols track')
            else:
                return None

        # this list comprehension checks start and end are valid options
        self._protein._check_position_is_valid(position, helper_string = f'Invalid position [{position}] passed to track {str(self)}')
        return self._values[position]


    ## ------------------------------------------------------------------------
    ##
    def symbol(self, position, safe=True):
        """
        Returns a single symbol from the passed position
        Parameters
        ----------
        position : int
            Starting position of interest.

        safe : bool (default = True)
            Flag which if true with throw an exception if a
            symbol is requested from a symbol track.

        Returns
        --------
        str
            Returns a symbol associated with the passed position

        """
        if self.symbols is None:
            if safe is True:
                raise TrackException('Requesting symbol from a values track')
            else:
                return None

        # this list comprehension checks start and end are valid options
        self._protein._check_position_is_valid(position, helper_string = f'Invalid position [{position}] passed to track {str(self)}')
        return self._symbols[position]


      
    ## ------------------------------------------------------------------------
    ##      
    def __repr__(self):             
        return "Track [name: %s] associated with protein %s" % (self.name, self.protein)

    ## ------------------------------------------------------------------------
    ##      
    def __len__(self):
        return len(self._protein)

    
    
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
        attribute associated with this Track.

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

        NOTE: Track attributes cannot be loaded or saved to file when
        Tracks are read/written via interfaces.si_track.

        
        Parameters
        ----------------
        name : str
            The attribute name. A list of valid names can be found by 
            calling the  ``<Track>.attributes()`` (which returns a list 
            of the valid names).
            

        safe : bool (default = True)
            Flag which if true with throw an exception if an attribute 
            with the same name  already exists.
                        
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
                raise TrackException('Requesting attribute [%s] from protein [%s] but this attribute has not been assigned' % (name, str(self))) 

            # if safe not passed just return None
            else:
                return None
                


    ## ------------------------------------------------------------------------
    ##
    def add_attribute(self, name, val, safe=True):
        """
        Function that adds an attribute. Note that if safe is true, this 
        function will raise an exception if the attribute is already 
        present. If safe=False, then an existing value will be 
        overwritten.

        NOTE: Track attributes cannot be loaded or saved to file when
        Tracks are read/written via interfaces.si_track.

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
            None - but adds an attribute to the calling object.

        """

        if safe:
            if name in self._attributes:
                raise TrackException("Trying to add attribute [%s=%s] to Track [%s] but this attribute is already set.\nPossible options are: %s" %(name,val, str(self), str(self._attributes.keys())))
                
        self._attributes[name] = val


    ## ------------------------------------------------------------------------
    ##
    def remove_attribute(self, name, safe=True):
        """
        Function that removes a given attribute from the Track based on the 
        passed attribute name. If the passed attribute does not exist or is not 
        associate with the Track then this will trigger an exception 
        unless safe=False.

        Parameters
        ----------------

        name : str
            The parameter name that will be used to identify it

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
                raise TrackException(f'Passed attribute [{name}] not found in {self}')
        else:
            del self._attributes[name]
