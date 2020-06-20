"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

import numpy as np
from . import exceptions
from . import sequence_utilities
from .domain import Domain 
from .site import Site
from .track import Track
from .exceptions import ProteinException
from .import general_utilities


STANDARD_AAs = 'ACDEFGHIKLMNPQRSTVWY'

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Class that defines a Entry
#
class Protein:
    """
    Protein objects are the parent object to all sequence-based information. Protein objects 
    are explicitly associated with several different types of objects:

    * **tracks** - Vectorial information that maps to each residue and contains some set of information. A protein can have multiple tracks, but there must be a one-to-one mapping for sequence to track.

    * **domains** - Information on a single contiguous region in the protein. A protein can have multiple domains.
    
    * **sites** - Information associated with a single amino acid site. A protein can have multiple sites
   
    * **attributes** - Protein-specific information associated 
                     
    Parameters
    ------------

    seq : str
        Amino acid sequence for the protein. No validation is performed. 

    name : str
        Some sort of name-based identifier for the protein. Can be anything - is not used 
        internally so no restraints, but could be used by other bits of analysis. 

    unique_ID : str 
        The unique_ID should be a short unique identifier. When added to a Proteome the 
        Proteome object ensures the unique_ID is unique with respect to that Proteome. 

        We HIGHLY recommend using the uniprot accession number, as this meets the 
        requirement of a unique ID as well as allowing effective cross-refering from 
        other databases.

    attribute_dictionary : dict (optional)
        The attribute_dictionary provides a key-value pairing for arbitrary information.
        This could include gene names, different types of identifies, protein copy number,
        a set of protein partners, or anything else one might wish to associated with the
        protein as a whole.

    """
    
    def __init__(self, seq, name, proteome, unique_ID = None, attribute_dictionary = None):
        """
        """
        
        # define internal attributes that are then accessed via @properties 
        self._name     = name
        self._sequence = "-" + seq # the '-' at the start fixes our indexing woes
        self._unique_ID = unique_ID
        self._proteome = proteome

        self._len = len(seq) # protein length 
        self._true_len = len(self._sequence) # length of string the protein is in

        # set attribute dictionary IF a dictionary was passed
        if isinstance(attribute_dictionary, dict):
            self._attributes = attribute_dictionary

        # set dictionary to an empty dictionary if none was passed
        elif attribute_dictionary is None:
            self._attributes = {}

        else:
            raise exceptions.ProteinException('[FATAL]: If provided, protein attribute must a dictionary')
        
        # initialize the empty dictionaries for the set of sites, domains and tracks
        self._sites   = {}
        self._domains = {}    
        self._tracks  = {}

        # the domains by type and sites by type dictionaries are only built IF we request
        # domains or sites by type. This provides a mode of conditional memoization, so at 
        # least within a single session we do not have to search through domains and sites
        # multiple times to find a specific type
        self._domains_by_type ={}
        self._sites_by_type ={}




    ## ------------------------------------------------------------------------
    ##
    @property
    def unique_ID(self):
        """
        [Property]: Returns the protein's unique_ID

        Returns
        ---------------
        str
            Returns the protein's unique_ID

        """
        return self._unique_ID


    ## ------------------------------------------------------------------------
    ##
    @property
    def name(self):
        """
        [Property]: Returns the protein name.

        Returns
        ---------------
        str
            Returns a string that corresponds to the region of interest
        """

        return self._name


    ## ------------------------------------------------------------------------
    ##
    @property
    def proteome(self):
        """
        [Property]: Returns the ``Proteome`` object this protein is associated with.

        Returns
        --------
        Proteome 
            Returns a Proteome object that contains this Protein.

        """
        return self._proteome


    ###################################
    ##                               ##
    ##      SEQUENCE FUNCTIONS       ##
    ##                               ##
    ###################################

    ## ------------------------------------------------------------------------
    ##
    def get_residue(self, position):
        """
        Function that returns the natural residue found at a given position

        Parameters
        ----------
        position : int
            Position of interest

        Returns
        ----------
        str
            Returns a single character that corresponds to the string of interest

        """

        # cast to integer incase...
        return self._sequence[int(position)]



    ## ------------------------------------------------------------------------
    ##
    @property
    def sequence(self):
        """
        [Property]: Returns the protein amino acid sequence as a string. Recall 
        that in strings indexing occurs from 0 and is non-inclusive. For 
        proteins/biology indexing is from 1 and is inclusive.

        i.e. for sequence 'MAPSTA...' real/biological indexing of region 1-3 would give
        you 'MAP' while Python's indexing would give you 'AP'.

        As a result BEWARE if using the raw sequence for analysis! The Protein class
        provides a get_sequence_region(), get_sequence_context() and analogous functions
        for tracks that allow you to use normal indexing to select ranges or regions
        around a specific point. We suggest this is a safer way to extract vectorial
        information.

        Returns
        --------
        str 
            Amino acid sequence associated with the protein.

        """
        return self._sequence[1:]


    ## ------------------------------------------------------------------------
    ##
    def get_sequence_region(self, start, end):
        """
        Function that allows a region of the sequence to be extracted out.

        Parameters
        ---------------
        start : int
            Start position for region

        end : int
            End position for region (note this is inclusive)

        Returns
        ---------------
        str
            Returns a string that corresponds to the region of interest

        """

        self._check_position_is_valid(start, helper_string='Sequence region cannot start below 1 [%i]'%(start))
        self._check_position_is_valid(end, 'Sequence region cannot end after the sequence length (%i) [%i]'%(self._len, end))
            
        # note +1 because we're inclusive with positions
        return self._sequence[start:end + 1]


    ## ------------------------------------------------------------------------
    ##
    def get_sequence_context(self, position, offset=5):
        """
        Function that allows a local region of the sequence centered on a specific position
        to be extracted, including +/- an offset border that intelligently truncates if the 
        offset would extend outside the sequence region.

        Parameters
        ---------------
        position : int
            Position for which we'll interrogate the local sequence

        offset : int
            [**Default = 5**] Plus/Minus offset used to investigate the region around the position. Note that
            an offset is symmetrical around the position. 

        Returns
        ---------------
        str
            Returns a string that corresponds to the region of interest

        """
        
        # sanity check position input
        self._check_position_is_valid(position, helper_string='Sequence position is outside below 1 [%i]')
        
        # compute start/end of context according to the offset
        (p1, p2) = sequence_utilities.get_bounding_sites(position, offset, self._len)

        # note +1 because we're inclusive with positioning here (and index from 1)
        return self._sequence[p1:p2 + 1]




    ## ------------------------------------------------------------------------
    ##
    def check_sequence_is_valid(self):
        """
        Function that checks if the current protein sequence is valid (i.e. consists 
        of only the standard 20 amino acids).
        
        Returns
        ---------------
        bool
            Returns True if all residues are in the standard 20 amino acids, and
            False if not
        
        """

        for i in self._sequence:
            if i not in STANDARD_AAs:
                return False
        return True


    ## ------------------------------------------------------------------------
    ##
    def convert_to_valid(self, copy=False, safe=True):
        """
        Function that converts non-standard amino acid residues to standard ones and
        applies this version to the Protein's sequence.

        Specifically:

        ``B -> N``

        ``U -> C``

        ``X -> G``

        ``Z -> Q``

        ``* -> <empty string>``

        ``- -> <empty string>``

        By default this alters the underlying sequence. If you wish to return a copy
        of the altered sequence instead set copy=True. Otherwise the underlying sequence
        is changed. Note that removing the ``*`` and ``-`` characters will change the 
        sequence length which could cause major issues as none of the internal position-
        specific references will automatically update. Note that if safe=True such changes
        will trigger an exception.

        Parameters
        ---------------
        copy : bool (default = False)
            Boolean flag - if set to true a copy of the updated sequence is returned, 
            if False then the function returns None. In both cases the associated 
            protein's sequence is  altered.
        
        safe : bool (default = True)
            Boolean flag that defines how to respond if an update changes the sequence
            length. If set to true, a change that alters the sequence length will
            trigger an exception, if False it will continue unannounced.

        Returns
        ---------------
        Unknown
            Depends on if copy is set to True or not. If copy=True 
        
        """

        if copy is True:
            # create a copy, such that within the protein the underlying sequence is
            # unaltered
            s = self._sequence[:]            
        else:
            # create a view, which means the protein sequence is altered
            s = self._sequence

        old_len=len(s)
        
        # systematically replace common  'non-canonical' one-letter codes
        # with acceptable codes. Code explanations from
        # https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp
        s = s.replace('B','N') # B = N/E
        s = s.replace('U','C') # U = selenocysteine
        s = s.replace('X','G') # X = Any
        s = s.replace('Z','Q') # Z = Q/D
        s = s.replace('*','')  # * = stop
        s = s.replace('-','')  # - = gap

        if safe is True:
            if len(s) != old_len:
                raise ProteinException('When altering the sequence to remove invalid characters the sequence length changed. This will invalidate positional attributes, such as sites and domains, as these are not automatically updated!')

        if copy is True:
            return s
        else:
            return None


    ## ------------------------------------------------------------------------
    ##
    def _check_position_is_valid(self, position, helper_string=None):
        """
        Internal function that ensures a passed position falls inside the sequence.
        The helper-string allows the exception raised to be customized
        """

        # recal we're operating in a 1-indexed space
        if sequence_utilities.inside_region(1, self._len, position):
            return None
        else:
            if helper_string:                
                raise ProteinException(helper_string)
            else:
                raise ProteinException('Position %i falls outside of sequence'%(position))


    
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
        with this protein.

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
             ``<Protein>.attributes()`` (which returns a list of the valid names)

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
                raise ProteinException('Requesting attribute [%s] from protein [%s] but this attribute has not been assigned' % (name, str(self))) 

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
                raise ProteinException("Trying to add attribute [%s=%s] to protein [%s] but this attribute is already set.\nPossible options are: %s" %(name,val, str(self), str(self._attributes.keys())))
                
        self._attributes[name] = val




    ###############################
    ##                           ##
    ##     TRACK FUNCTIONS       ##
    ##                           ##
    ###############################

    ## ------------------------------------------------------------------------
    ##
    @property
    def tracks(self):
        """
        **[Property]**: Provides a list of the keys associated with each track 
        associated with this protein.
        
        Returns
        -------
        list
            returns a list of the track keys associated with the protein. 


        """
        return list(self._tracks.keys())


    ## ------------------------------------------------------------------------
    ##
    def track(self, name, safe=True):

        """
        Function that returns a specific track as defined by the name. 

        Recall that tracks are defined by a name. If a track by this name exists
        this function returns the actual ``Track`` object, NOT the **values** or **symbols**
        associated with the track. If a track by this name does *not* exist then if safe=`True`
        an exception will be raised, otherwise the function returns None.
        
        For direct access to values and symbols, use the ``<Protein>.get_track_values(<track_name>)`` and 
        ``<Protein>.get_track_symbols(<track_name>)``.

        Parameters
        ----------------
        name : str
             The track name. A list of valid names can be found by calling the
             ``<Protein>.tracks()`` (which returns a list of the valid track names)

        Returns
        ---------
        Unknown 
            Will either return the Track object associated with the name, OR
            will return None if safe=False and there was no Track object that
            matched the name.
        
        """

        if name in self._tracks:
            return self._tracks[name]

        elif safe:
            raise exceptions.ProteinException('No track named [%s] in protein %s\n\nAvailable options are: %s' %(name, str(self), self.tracks))



    ## ------------------------------------------------------------------------
    ##
    def get_track_values(self, name, start=None, end=None, safe=True):
        """
        Function that returns the values associated with a specific track, as 
        defined by the name.

        Recall that tracks are defined by a name. If a track by this name exists
        this function returns the values IF these are associated with the track.
        If no values are associated then the function will throw an exception 
        unless safe is set to False, in which case it will return None.
        
        Parameters
        ----------------
        name : string
             The track name. A list of valid names can be found by calling the
             <Protein>.tracks (which returns a list of the valid track names)

        start : int (default none)
            If provided defines the start position along the track. Note
            if start is provided end must also be provided.

        end : int (default none)
            If provided defines the start position along the track. Note
            if end is provided start must also be provided.

        safe : bool (default = True)
            Flag which if true with throw an exception if a track that matches
            the passed name does not already exist.
            
        Returns
        ---------
        Unknown 
            Will either return the values associated with the track, OR
            will return None if safe=False and there was no Track that
            matched the name.
        
        """

        (_start, _end) = self._build_start_end(start,end)
        
        # call internal function
        return self.__get_track_info(name, safe, _start, _end, 'values')



    ## ------------------------------------------------------------------------
    ##
    def get_track_symbols(self, name, start=None, end=None, safe=True):
        """
        Function that returns the symbols associated with a specific track, as 
        defined by the name.

        Recall that tracks are defined by a name. If a track by this name exists
        this function returns the symbols IF these are associated with the track.
        If no symbols are associated then the function will throw an exception 
        unless safe is set to False, in which case it will return None.

        Parameters
        ----------------
        name : string
             The track name. A list of valid names can be found by calling the
             <Protein>.tracks (which returns a list of the valid track names)

        start : int (default none)
            If provided defines the start position along the track. Note
            if start is provided end must also be provided.

        end : int (default none)
            If provided defines the start position along the track. Note
            if end is provided start must also be provided.

        safe : bool (default = True)
            Flag which if true with throw an exception if a track that matches
            the passed name does not already exist.
            
        Returns
        ---------
        Unknown 
            Will either return the values associated with the track, OR
            will return None if safe=False and there was no Track that
            matched the name.
        
        """

        # build the start and end position
        (_start, _end) = self._build_start_end(start,end)

        # call internal function
        return self.__get_track_info(name, safe, _start, _end, 'symbols')



    ## ------------------------------------------------------------------------
    ##
    def __get_track_info(self, name, safe, start, end, mode):
        """
        Internal function that follows the exact same logic as the public-facing
        get_track_symbols or get_track_values

        """
        
        t = self.get_track(name, safe)

        if t is None:
            # note - technically as the code is written now we don't need this, (the safety is dealt in get_track())
            # but I'm including it for best practice to avoid implicit dependencies in the code
            if safe:
                raise exceptions.ProteinException('No track named [%s] in protein %s\n\nAvailable options are: %s' %(name, str(self), self.tracks))
            else:
                return None

        # try and get values
        if mode == 'values':
            v = t.values_region(start, end)
        elif mode == 'symbols':
            v = t.symbols_regions(start, end)

        # if v is a value or safe is False just return v (will either be values
        # or None)
        if v is not None or not safe:
            return v

        # we only get here if v is None and safe is True
        else:
            raise exceptions.ProteinException('Requested track values for track [%s] in protein [%s] but no values available' %(t.name, self))
        


    ## ------------------------------------------------------------------------
    ##
    def add_track(self, name, values=None, symbols=None, safe=True):
        """
        Function that adds a track to this protein. For more information on
        Tracks see the relevant documentation. However, some general guidelines
        are provided below for convenience.

        (1) A values track should be a list/array of numerical values

        (2) A symbols track should be a list or string of symbolic characters

        (3) In either case, the iterable should have a 1:1 mapping with the sequence

        (4) Tracks can have both a value and a symbol, although in general it probably
            makes sense to use multiple tracks. 
        
        Parameters
        ---------------

        name : string
            Name for track. NOTE that this is a unique identifier, and each track within
            a given protein should must have a unique name. 

        values : list or np.array (default None)
            A numerical iterable collection of values, where each value maps to a specific 
            residue in the sequence. 

        symbols : list or string (default None)
            A symbolic collection of characters, where each symbol maps to a specific 
           residue in the sequence. 
        
        safe : boolean (default = True)
            If set to True over-writing tracks will raise an exception, otherwise overwriting
            a track will simply over-write it.

        Returns
        ----------

            Nothing, but adds a track to the calling object.

        """

        if name in self.tracks:

            if safe is True:
                raise exceptions.ProteinException('Trying to add Track [%s] in protein [%s] but Track already exists' % (name, self.name))
                
        self._tracks[name] = Track(name, self, values, symbols)


    ## ------------------------------------------------------------------------
    ##
    def build_track_values_from_sequence(self, name, trackfunction, input_dictionary={}, safe=True):
        """
        Tracks can be added as pre-loaded values. However, sometimes you want to build a track based on some
        analysis of the sequence on the fly. This function allows you to pass in your own function (with keyword
        arguments in the keywords dictionary) that will take in the protein sequence, generate a new track, and 
        add that track to the protein.

        build_track_values allows you to define a function that converts amino acid sequence into a numerical list
        or np.array, which gets written as a values track. If you want a symbols track, use build_track_symbols().

        Specifically, the argument trackfunction must be a user-defined function. This function can be defined anywhere, 
        but should take either one or two arguments:

        (1) The first/only argument should be an amino acid sequence.
        (2) The second argument a dictionary of key-value pairs.

        When build_track_values_from_sequence is called, the sequence of the protein is passed as the first
        argument into the trackfunction, and  - if present - the input_dictionary is passed as the second argument.

        In this way a new track is defined internally, with the trackfunction using the proteins sequence and 
        any/all pass input_dictionary to convert the sequence into some numerical representation.
        

        Parameters
        ------------

        name : string
            Name of the track to be used. Should be unique and will always overwrite an existing track with the same
            name (no safe keyword provided here).

        trackfunction : function
            A user define function that has the following properties:
        
            (1) First argument is expected to be amino acid sequence
            (2) Second argument (if provided) should be a dictionary which is passed (untouched) THROUGH build_track_values 
                from sequence to the trackfunction at runtime

        function_keywords : dictionary
            This is a dictionary that will be passed to the trackfunction as the second argument IF it is provided.
            In this way, the user can pass an arbitrarily complex set of arguments to the trackfunction each time 
            the build_track_values_from_sequence is called.

        safe : boolean (default = True)
            If set to True over-writing tracks will raise an exception, otherwise overwriting
            a track will simply over-write it.

        Example
        --------
        
        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        # define a function that takes in a sequence and converts it into some other numerical
        # list. Note this is INLINE with the code, or could be elsewhere. This function MUST
        # take either ONE argument (sequence) or TWO arguments (sequence and input_dictionary).
        # Also the names of these arguments does not matter, but the order does (i.e. first argument
        # will always get the sequence).

        def trackbuilder(seq, input_dictionary):
            ''' 
                This function takes in a sequence (seq) as first argument, and the v1 and v2
                as additional arguments. See below for what it's doing (pretty simple).
            '''
            newseq=[]  
            
            # we are extracting out the 'values' from the input dictionary
            # for the sake of code clarity
            v1 = input_dictionary['v1']
            v2 = input_dictionary['v2']

            # for each residue in the sequence
            for i in seq:

                # is that residue in v1 (append 1) or v2 (append -1)? If neither append 0
                if i in v1:
                    newseq.append(1)
                elif i in v2:
                    newseq.append(-1)
                else:
                    newseq.append(0)
        
            return newseq
        
        # define the input_dictionary (note again that the variable names here do not matter)
        input_dictionary = {'v1':['K','R'], 'v2':['E','D']}

        # now assuming ProtOb is a Protein object, this will add a new track
        ProtOb.build_track_values('charge_vector', trackbuilder, function_dictionary=input_dictionary)

        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


        In this example we defined a function that converts an amino acid string into a numerical list where
        positively charged residues = +1 and negatively charged residues = -1. We applied this function to generate
        a 'charge_vector' track.


        Note this is analagous to defining our function and then running:

        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        s = ProtOb.sequence
        newtrack = trackbuilder(s, ['K','R'], ['E',D'])
        ProbOb.add_track('charge_vector', values=newtrack)

        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


        Some FAQs:

        Do I need to pass an input_dictionary to the custom function?
            No!

        Does the name of the custom function matter?
            No!

        Does the custom function have to accepted the amino acid sequence as the first argument?
            Yes!

        """

        # if this will overwrite an existing track and safe is on...
        if name in self.tracks:
            if safe is True:
                raise exceptions.ProteinException('Trying to add Track [%s] in protein [%s] but Track already exists' % (name, self.name))

        # build the new track with the trackfunction, correctly handling between 0 and n additional
        # arguments to be passed to the trackfunction
        if len(input_dictionary) == 0:
            built_track = trackfunction(self.sequence)
        else:
            built_track = trackfunction(self.sequence, input_dictionary)
            
        # finally add the track
        self._tracks[name] = Track(name, self, values=built_track, symbols=None)



    ## ------------------------------------------------------------------------
    ##
    def build_track_symbols_from_sequence(self, name, trackfunction, input_dictionary={}, safe=True):
        """
        Tracks can be added as pre-loaded values. However, sometimes you want to build a track based on some
        analysis of the sequence on the fly. This function allows you to pass in your own function (with keyword
        arguments) that will take in the protein sequence, generate a new track, and add that track to the Protein.

        build_track_symbols allows you to define a function that converts amino acid sequence into a symbolic list
        or string, which gets written as a symbols track. If you want a values track, use build_track_values().

        Specifically, the argument trackfunction must be a user-defined function. This function can be defined anywhere, 
        but should take either one or two arguments:

        (1) The first/only argument should be an amino acid sequence.
        (2) The second argument a dictionary of key-value pairs.

        When build_track_symbols_from_sequence is called, the sequence of the protein is passed as the first
        argument into the trackfunction, and  - if present - the input_dictionary is passed as the second argument.

        In this way a new track is defined internally, with the trackfunction using the proteins sequence and 
        any/all pass input_dictionary to convert the sequence into some other symbolic representation.
        

        Parameters
        ------------

        name : string
            Name of the track to be used. Should be unique and will always overwrite an existing track with the same
            name (no safe keyword provided here).


        trackfunction : function
            A user define function that has the following properties:
        
            (1) First argument is expected to be amino acid sequence
            (2) Second argument (if provided) should be a dictionary which is passed (untouched) THROUGH build_track_values 
                from sequence to the trackfunction at runtime

        function_keywords : dictionary
            This is a dictionary that will be passed to the trackfunction as the second argument IF it is provided.
            In this way, the user can pass an arbitrarily complex set of arguments to the trackfunction each time 
            the build_track_symbols_from_sequence is called.

        safe : boolean (default = True)
            If set to True over-writing tracks will raise an exception, otherwise overwriting
            a track will simply over-write it.


        Example
        --------        
        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        # define a function that takes in a sequence and converts it into some other symbolic representation
        # as a string. Note this is INLINE with the code, or could be elsewhere. This function MUST
        # take either ONE argument (sequence) or TWO arguments (sequence and input_dictionary).
        # Also the names of these arguments does not matter, but the order does (i.e. first argument
        # will always get the sequence).

        def trackbuilder(seq, input_dictionary):
            ''' 
                This function takes in a sequence (seq) as first argument, and the v1 and v2
                as additional arguments. See below for what it's doing (pretty simple).
            '''
            new_string_list=[]  
            
            # we are extracting out the 'values' from the input dictionary
            # for the sake of code clarity
            v1 = input_dictionary['v1']
            v2 = input_dictionary['v2']

            # for each residue in the sequence
            for i in seq:

                # is that residue in v1 (append 1) or v2 (append -1)? If neither append 0
                if i in v1:
                    new_string_list.append('+')
                elif i in v2:
                    new_string_list.append('-')
                else:
                    new_string_list.append('0')
        
            # convert the list into a string
            newstring = "".join(new_string_list)
            return newstring
        
        # define the input_dictionary (note again that the variable names here do not matter)
        input_dictionary = {'v1':['K','R'], 'v2':['E','D']}

        # now assuming ProtOb is a Protein object, this will add a new track
        ProtOb.build_track_values('charge_string', trackbuilder, function_dictionary=input_dictionary)

        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        In this example we defined a function that converts an amino acid string into a coarse-grained string 
        representation where positive residues are "+", negative are "-" and neutral are "0".

        Note this is analagous to defining our function and then running:

        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        s = ProtOb.sequence
        newtrack = trackbuilder(s, ['K','R'], ['E',D'])
        ProbOb.add_track('charge_vector', values=newtrack)

        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        Some FAQs:

        Do I need to pass an input_dictionary to the custom function?
            No!

        Does the name of the custom function matter?
            No!

        Does the custom function have to accepted the amino acid sequence as the first argument?
            Yes!

        """

        # if this will overwrite an existing track and safe is on...
        if name in self.tracks:
            if safe is True:
                raise exceptions.ProteinException('Trying to add Track [%s] in protein [%s] but Track already exists' % (name, self.name))

        # build the new track with the trackfunction, correctly handling between 0 and n additional
        # arguments to be passed to the trackfunction
        if len(input_dictionary) == 0:
            built_track = trackfunction(self.sequence)
        else:
            built_track = trackfunction(self.sequence, input_dictionary)
            
        # finally add the track
        self._tracks[name] = Track(name, self, values=None, symbols=built_track)



    ## ------------------------------------------------------------------------
    ##
    def build_track(self, name, input_data, track_definition_function, safe=True):
        """
        Function that constructs a track using a given track_definition_function and 
        a user provided input_data object. Very little constraint is set here, other
        than the fact the name should be a string and track_definition function should
        return a dictionary with (at least) two key:value pairings: `symbols` and `values`,
        where the corresponding value for each is bona-fide track input data.
        
        Parameters
        ------------

        name : string
            Name of the track to be used. Should be unique and will always overwrite an existing track with the same
            name (no safe keyword provided here).

        input_data : ?
            Some kind of data that will be passed to the track_definition_function 
        
        track_definition_function : function
            Function that takes in `input_data` and returns a dictionary with a 'values' and a 'symbols' key and value
            pairing. The values that map to 'values' and 'symbols' will be added as a single new track defined by name.

        safe : boolean (default = True)
            If set to True over-writing tracks will raise an exception, otherwise overwriting
            a track will simply over-write it.


        """

        if name in self.tracks:
            if safe is True:
                raise exceptions.ProteinException('Trying to add Track [%s] in protein [%s] but Track already exists' % (name, self.name))


        track_out = track_definition_function(input_data)

        values = track_out['values']
        symbols = track_out['symbols']

        self._tracks[name] = Track(name, self, values=values, symbols=symbols)

        



    ###############################
    ##                           ##
    ##     DOMAIN FUNCTIONS      ##
    ##                           ##
    ###############################

    ## ------------------------------------------------------------------------
    ##
    @property
    def domains(self):
        """
        **[Property]**: Returns a list of the domain names associated with this protein

        """
        return list(self._domains.keys())


    ## ------------------------------------------------------------------------
    ##
    @property
    def domain_types(self):
        """
        **[Property]**: Returns a list of the unique domain types associated with 
        this protein. There will be no duplicates here.

        """

        # define an empty set
        domain_types = set([])

        # cycle through the domains and add the domain type to the set
        for domain in self._domains:
            domain_types.add(self._domains[domain].domain_type)

        # convert the set to a list and return
        return list(domain_types)


    ## ------------------------------------------------------------------------
    ##
    def domain(self, name, safe=True):

        """
        Function that returns a specific domain as defined by the name. Note it is
        often more useful to request a domain by type rather than by the name, in
        which case get_domain_by_type(<domain_type>) is the relevant syntax.
        
        Note domains can also be requested based on position (get_domain_by_position).

        Parameters
        ----------------
        name : string
             The Domain name. A list of valid names can be found by calling the
             <Protein>.domains (which returns a list of the valid track names)

        safe : bool (default = True)
            Flag which if true with throw an exception if no domain exists with
            this name. If false function will return None instead.

        Returns
        ---------
        Unknown 
            Will either return the Domain object associated with the name, OR
            will return None if safe=False and there was no Domain object that
            matched the name.
        
        """

        if name in self._domains:
            return self._domains[name]
        elif safe:
            raise exceptions.ProteinException('No domains named [%s] in protein %s\n\nAvailable options are: %s' %(name, self, str(self.domains)))
      
      
    ## ------------------------------------------------------------------------
    ##
    def add_domains(self, list_of_domains, safe=True, autoname=False):
        """

        Function that takes a list of domain dictionaries and adds those domains to the protein.

        Each domain dictionary within the list must have a key-value pair that defines the following info:
                    
        start : domain start position (in real sequence, not i0 indexing)
        end   : domain end position (in real sequence, not i0 indexing)
        domain_type : type of the domain (string)
        attribute_dictionary : a dictionary of attributes to associated with the domain (optional)
        
        Note that in start, end, and domain_type are the only required key-value pairs required in
        the dictionary.

        Parameters
        -------------

        list of domains : list of domain dictionaries
            A "domain dictionary" is defined above, but in short is a dictionary with the following
            key-value pairs:

            REQUIRED:
            start                : int (domain start position)
            end                  : int (domain end position)
            domain_type          : string (domain type)

            OPTIONAL:
            attribute_dictionary : dictionary of arbitrary key-value pairs 
                                   that will be associated with the domain
            

        safe : boolean 
            If set to True over-writing domains will raise an exception. If False, overwriting
            a domain will silently over-write. Default = True.

        autoname : boolean
            If autoname is set to true, this function ensures each domain ALWAYS has a unique
            name - i.e. the allows for multiple domains to be perfectly overlapping in position
            and type. This is generally not going to be required and/or make sense, but having
            this feature in place is useful. In general we want to avoid this as it makes it 
            easy to include duplicates which by default are prevented when autoname=False.


        Returns
        -------

        No return value, but will add the passed domains to the protein or throw an exception
        if something goes wrong!

        """

        # for each domain definition check we can extract the relevant info 
        for new_domain in list_of_domains:

            try:
                # note we split these operations into two blocks to aid with
                # debugging if someone makes an invalid domain dictionary
                start = new_domain['start']
                end = new_domain['end']
                domain_type = new_domain['domain_type']
                
                # cast to integers
                start=int(start)
                end = int(end)

            except Exception:
                raise ProteinException("When trying to add domains with the add_domains function one or more dictionaries in the return list lacked a 'start', 'end', or 'domain_type' key, and/or the start/end values could not be correctly coverted to integers")
                
            if 'attribute_dictionary' in new_domain:
                ad = new_domain
            else:
                ad = {}

            self.add_domain(start, end, domain_type, ad, safe, autoname)


        
    ## ------------------------------------------------------------------------
    ##
    def add_domain(self, start, end, domain_type, attribute_dictionary=None, safe=True, autoname=False):
        """
        Function that adds a domain, automatically generating a unique name if 
        none is provided. Domain type can be used to assign a specific type if
        we want to retrieve domains of a specific type at some point. Position
        indexing is done for 1 - i.e. the first residue in a protein is 1, not
        0.

        Allows a domain at a specific position to be 

        Parameters
        -----------

        start : int
            Position of the start of the domain, inclusive.

        end : int 
            Position of the end of the domain (not inclusive). i.e. if we had a
            domain that ran from start=10 end=20, it would be 10 residues long
            and include residues [10, 11, 12, 13, 14, 15, 16, 17, 18, 19].

        domain_type : string 
            None unique string that allows a type identifier to be associated 
            with a domain. 

        attribute_dictionary : dictionary [optional]
            Optional dictionary which allows an arbitrary set of attributes to be
            associated with a domain, in much the same way that they can be associated
            with a protein. Default = None.

        safe : boolean 
            If set to True over-writing tracks will raise an exception, otherwise overwriting
            a track will simply over-write it. Default = True.

        autoname : boolean
            If autoname is set to true, this function ensures each domain ALWAYS has a unique
            name - i.e. the allows for multiple domains to be perfectly overlapping in position
            and type. This is generally not going to be required and/or make sense, but having
            this feature in place is useful. In general we want to avoid this as it makes it 
            easy to include duplicates which by default are prevented when autoname=False.

        """

        # append start and end position to name. 
        full_name = "%s_%i_%i"%(domain_type, start, end)    
        
        # if this domain name was already found...
        if full_name in self.domains:

            # if we're in autoname mode create a new unique name. This acts to add an incrementer to the
            # end and cycles through until a unique domaintype_star_end_incrementer name is found, where 
            # incremementer keeps being incremented
            if autoname:
                increment = 0
                found = False

                while found is False:

                    increment = increment + 1
                    newname = "%s_%i_%i_%i"%(domain_type, start, end, increment)    

                    if newname not in self.domains:
                        found = True                    
                full_name = newname
            elif safe:
                raise exceptions.ProteinException('Domain [%s] already found in proteins %s' % (full_name, self.name))
            
        self._domains[full_name] = Domain(start, end, self, domain_type, attribute_dictionary=attribute_dictionary)



    ## ------------------------------------------------------------------------
    ##
    def build_domain(self, input_data, domain_definition_function,  safe=True, autoname=False):
        """
        Function that is somewhat analogous to build_tracks, but allows the user to define a custom function 
        (domain_definition_function) that takes input_data and returns domain information, and then assigns
        that domain information to a new domain


        Parameters
        ----------
        input_data : anything
            Any input that makes sense when passed to the domain_definition_function

        domain_definition_function : function that takes a single argument (input_data)
            Function that takes input_data as it's input argument and returns a list of 0 or more dictionaries. 
            Each dictionary within the list has a key-value pair that defines the following info:
                    
            start : domain start position
            end   : domain end position
            domain_type : type of the domain
            attribute_dictionary : a dictionary of attributes to associated with the domain (optional)
            
            Note that in principle only start and end are required, although we highly recommend  one/both of domain name and 
            domain type.

            Some important requirements to consider
            (1) domain_definition_function must return a list of zero or more dictionaries
            

        safe : bool (default = True)
            Flag which if true with throw an exception of a domain with the same name already exists

        autoname : boolean
            If autoname is set to True, this function ensures each domain ALWAYS has a unique
            name - i.e. the allows for multiple domains to be perfectly overlapping in position
            and type. This is generally not going to be required and/or make sense, but having
            this feature in place is useful. In general we want to avoid this as it makes it 
            easy to include duplicates which by default are prevented when autoname=False.


        """

        # build our domain definitions. Should probably add code that ensyres
        # domain_definitions is a list
        domain_definitions = domain_definition_function(input_data)

        self.add_domains(domain_definitions)
            
                
    """
    ## To Do!
    ## ------------------------------------------------------------------------
    ##
    def get_domains_by_position(self, position, wiggle = 0):
        pass


    ## ------------------------------------------------------------------------
    ##
    def get_domains_by_position_range(self, start, end, wiggle = 0):
        pass
    """

    ## ------------------------------------------------------------------------
    ##
    def get_domains_by_type(self, domain_type, perfect_match=True):
        """
        Function that returns a dictionary of domains that match a specific type
        of domain. 
        
        Parameters
        ------------
        domain_type : string
            String associated domain_type that you want to search for

        perfect_match : boolean
            Flag that identifies if the domain names should be a perfect match (=True) 
            or if the string passed should just appear somewhere in the domain_type 


        Returns
        -----------
        dict
            Returns a dictionary where keys are domain_id and values are the actual domain objects
            that match the request
            
        """
        if perfect_match:
            def selection(t):
                if t  == domain_type: 
                    return True
                else:
                    return False
        else:
            def selection(t):
                if t.find(domain_type) > -1:
                    return True
                else:
                    return False

        return_dict = {}
        for domain_id in self.domains:
            d = self.domain(domain_id)
            if selection(d.domain_type):
                return_dict[domain_id] = d

        return return_dict
            
        

    ###############################
    ##                           ##
    ##     SITE FUNCTIONS        ##
    ##                           ##
    ###############################

    ## ------------------------------------------------------------------------
    ##
    @property
    def sites(self):
        """
        **[Property]**: Provides a list of the keys associated with every site on 
        the protein. Recall sites are indexed simply by their position, so this returns
        a list of positions associated with the protein where sites are found.

        """
        return list(self._sites.keys())


    ## ------------------------------------------------------------------------
    ##
    def site(self, position, safe=True):
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

        return self._sites[int(position)]
        

    ## ------------------------------------------------------------------------
    ##
    def add_site(self, position, site_type, symbol=None, value=None, attributes={}):
        """
        Function that adds a site to a specific position in the sequence. Sites are
        indexed by residue position, and multiple sites can co-exist on the same site,
        so no name is required (unlike Proteins, Tracks or Domains).
       
        site_type is a non-unique identifier that allows sites to be specifically
        identified/selected.

        Sites can be associated with a numerical value, a symbol, or both. Sites can
        also have attributes associated with them.
        
        Parameters
        -----------

        position : int
            Position of site (recall we index from 1 - i.e. the first residue in 
            a protein = 1, not 0. Note that this value is cast to int.

        site_type : string 
            Non-unique string that allows a type identifier to be associated 
            with a site.

        symbol : string
            Symbol associated with a site. Symbols are string-based - will often
            be a single character but could be multiple characters. Default = None.

        value : float64 
            Numerical value associated with a site. Note that the value is cast
            to a float64. Default = None

        attributes : dictionary [optional]
            Optional dictionary which allows an arbitrary set of attributes to be
            associated with a domain, in much the same way that they can be associated
            with a protein. Default = empty dictionary.
        
        """
        
        # recal inside_regions is inclusive
        if not sequence_utilities.inside_region(1, self._len, position):
            raise ProteinException("Trying to add site to protein [%s] at positions [%i] - this falls outside the protein's dimensions [%i-%i]" %(str(self), position, 1, self._len))

        # cast the position to an int and if there are no sites at that position create an empty list there
        position = int(position)        
        if position not in self._sites:
            self._sites[position] = []

        # add the site!
        self._sites[position].append(Site(position, site_type, self, symbol, value, attributes))


    ## ------------------------------------------------------------------------
    ##
    def get_sites_by_position(self, position, wiggle = 0):
        """
        Get all sites at a specific position

        Parameters
        ---------------
        position : int
            Residue position of interest (position in sequence)

        wiggle : int (default = 0)
            Value +/- the position (i.e. lets you look at sites around a 
            specific position)

        Returns
        -----------
        dictionary

            Returns a dictionary where key-value pairs are

            keys: positions along the sequence

            values: list of one or more sites found at that position

        """

        return self.get_sites_by_range(position, position, wiggle)


    ## ------------------------------------------------------------------------
    ##
    def get_sites_by_range(self, start, end, wiggle = 0): 
        """
        Get all sites within a certain range.

        Parameters
        ---------------
        start : int
            Start of region of interest (position in sequence)

        end : int
            End of region of interest (position in sequence)

        wiggle : int (default = 0)
            Value +/- at the edges that are included. 


        Returns
        -----------
        dictionary

            Returns a dictionary where key-value pairs are

            keys: positions along the sequence

            values: list of one or more sites found at that position

        """
       
        return_dict = {}

        # recal p1 and p2 should be in real-world indices 
        p1 = max(1, start - wiggle)
        p2 = min(end + wiggle, self._len)

        # recall we need +1 offset so we go to the end - positions/ranges are inclusive
        # when talking about proteins
        for j in range(p1, p2+1):
            if j in self._sites:
                return_dict[j] = self._sites[j]
            
        return return_dict


    ## ------------------------------------------------------------------------
    ##
    def get_sites_by_type(self, site_types):
        """
        Get a set of sites that match a specified site-type.

        Parameters
        ------------------
    
        site_types : string or list of strings
            One or more possible site_types that may be found in the protein. Either
            a single string or a list of strings can be passed, allowing for one or
            more sites to be grouped together

        Returns 
        ----------
        dict 
            Returns a dictionary where the key is a position (site) and the value
            is a list of one or more sites at that position that match the site
            type of interest. 
        
        """

        return self.__site_by_type_internal(self._sites, site_types)
        
        

    ## ------------------------------------------------------------------------
    ##
    def get_sites_by_type_and_range(self, site_types, start, end, wiggle=0):
        """
        Returns a set of sites that match both a type of interest and are found
        in the range provided. 

    
        Parameters
        ------------------
    
        site_types : string or list of strings
            One or more possible site_types that may be found in the protein. Either
            a single string or a list of strings can be passed, allowing for one or
            more sites to be grouped together

        start : int
            Start residue that defines start of region to be examined

        end : int
            End reidue that defines end of region to be examined

        wiggle : int
            Value that adds a bit of slack to the start/end positions. Default
            = 0.

        Returns 
        ----------
        dict 
            Returns a dictionary where the key is a position (site) and the value
            is a list of one or more sites at that position that match the site
            type of interest. 
      
        """

        # first get sites within the range
        initial_dict = self.get_sites_by_range(start, end, wiggle)
        
        # and then subselect sites of the right type
        return self.__site_by_type_internal(initial_dict, site_types)



    ## ------------------------------------------------------------------------
    ##
    def __site_by_type_internal(self, indict, site_types):
        """
        Internal function that allows a subset of sites to be selected based 
        on the passed site_type(s).

        Parameters
        ------------------
        site_types : string or list of strings
            One or more possible site_types that may be found in the protein. Either
            a single string or a list of strings can be passed, allowing for one or
            more sites to be grouped together

        Returns
        -----------
        dict 
            Returns a dictionary where the key is a position (site) and the value
            is a list of one or more sites at that position that match the site
            type of interest. This is exactly the same structure as the 
            self._sites dictionary, just filtered for a specific site_type.

        """

        # function that allows site_types to be either a string or a list
        # of strings so one or more sity_types can be passed
        site_types = general_utilities.string_to_list_of_strings(site_types)

        return_dict = {}

        # for each key (which reflects a site position in the passed dictionary)
        for i in indict:

            # for each site found in the list associated with that position
            for site_object in indict[i]:

                # for the one or more site types in the site_types list
                for ST in site_types:

                    # if that site type matches the target type
                    if site_object.site_type == ST:

                        # add that site to a new dictionary
                        if i in return_dict:
                            return_dict[i].append(site_object)
                        else:
                            return_dict[i] = [site_object]

        return return_dict

        

    ## ------------------------------------------------------------------------
    ##
    def __repr__(self):             
        return "| Protein: %s - L=%i, #t=%i, #d=%i, #s=%i, #a=%i |" %(self.unique_ID, self._len, len(self.tracks), len(self.domains), len(self.sites), len(self.attributes))
        


    ## ------------------------------------------------------------------------
    ##
    def __len__(self):             
        return self._len
        
