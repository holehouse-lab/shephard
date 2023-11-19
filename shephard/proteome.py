"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

from . import general_utilities
from .exceptions import ProteomeException
from .protein import Protein
from itertools import islice
import copy

class Proteome:
    """
    The Proteome object is the main unit for information storage in 
    SHEPHARD.

    There are a few ways that new Proteomes can be generated:

    * By reading in a FASTA file (using shephard.interfaces.apis.fasta)

    The Proteome constructor takes a single argument, which is a 
    list of protein dictionaries or a list of Protein objects. This 
    means Proteome objects can be generated directly (see below for 
    a definition of protein dictionaries). However, it is often more 
    convenient to build Proteomes from FASTA files. For more information 
    on this see the ``api`` documentation.
    
        
    Protein dictionaries are dictionaries that must contain four elements 
    (others are ignored).
    
    * ``sequence`` : *str* - Amino acid sequence of the protein.        

    * ``name`` : *str* - Name of the protein (this can be anything, it is not used internally so no constraints on what this is.

    * ``unique_ID`` :  *str* - This must be unique with respect to all other unique_IDs in the set of proteins in the input list.
                                          
    * ``attributes`` : *dict* - Dictionary of one or more attributes to apply to this protein. Key/value pairs in this dictionary can be arbitrary and are user defined.
    
    

    As an example:

    >>> protein_dictionary_example = {'sequence':'ALAPSLLPAMPALSPALSP', 
                                      'name': 'my protein fragment', 
                                      'unique_ID':'UXX01', 'attributes':{}}
    >>> dictionary_list = []
    >>> dictionary_list.append(protein_dictionary_example)
    >>> P = Proteome(dictionary_list)
                                          
    Note that ``sequence``, ``name`` and ``unique_ID`` are cast to *str* by 
    the function, so if numerical values are passed for any these will be 
    converted to strings.
    
    **Notes**

    * NOTE that ALL FOUR of these are required for EACH protein, even if the 
      attributes dictionary is empty.

    * The unique_ID is checked for uniqueness against all others in the Proteomes 
      and will throw and exception if it is, in fact, not unique.

    * Additional proteins can be added using the `.add_protein()` or 
      `.add_proteins() function.

    """

    ## ------------------------------------------------------------------------
    ##
    def __init__(self, input_list = None, attributes = None, force_overwrite=False):
        # See the Proteome class documentation for constructor info
        """
        Constructor that generates a new Proteome object. This includes 
        taking a list of Protein objects or protein dictionaries (input_list) 
        and, optionally an attributes dictionary for the proteome itself.
                
        In addition, force_overwrite can be used to deal with duplicate 
        entries in the input_list.

        There are two types of lists that can be tolerated when passed to
        the Proteome constructor: a list of protein dictionaries and a 
        list of Protein objects. Note this simply is using the 
        .add_proteins() function to add the passed input_list.

        **Protein dictionaries**        
        One mode of adding multiple proteins is by passing a list of 
        Protein dictionaries.

        Protein dictionaries are dictionaries that posses the following 
        key-value pairs ::

            'sequence'   : amino acid sequence (str)
            'name'       : protein name (str)
            'unique_ID'  : The unique identification number used for the 
                           protein (str)
            'attributes' : A dictionary of arbitrary key-value pairs to 
                           associate with the protein (dict or None)
        
        Additional keys/value pairs are ignored and ALL four of these must
        be included. If any are missing for any protein entry this function 
        raises a ProteomeException.

        **Protein objects**
        A second mode of adding multiple proteins is by passing a list of
        Protein objects. This is useful if you're creating a new Proteome
        based on a subset of proteins from an existing Proteome.

        In both cases, the function automatically determines the type of 
        the passed list, and adds dictionaries accordingly. Note that 
        in both cases proteins are added by value - i.e. a new Protein
        objects are generated, and changes to Proteomes in the new 
        Proteome will not affect the original Proteome.
        
        Parameters
        -----------
        input_list : list (default = None)
            List of Protein dictionaries

        attributes : dict (default = None)
            A an arbitrary set of key-value pairs to annotate the proteome 
            with metadata.
        
        force_overwrite : bool (default = False)
            If set to False and there are duplicate unique_IDs in protein 
            dictionaries in the input_list this will trigger an exception. 
            However, if set to True then the 'last' entry overwrites a             
            more recent one in the case of duplicates.

        """

        # initiallize book-keeping instruments
        self._records = {}
        self._unique_domain_types = {}
        self._unique_site_types = {}         # dict that mapes a Site type to a count
        self._unique_track_names = {}        # dict that mapes Track name to count
        self._track_name_to_track_type = {}  # dict that maps Track name to track type ('values' or 'symbols')
        
        # check attributs dictionary
        general_utilities.variable_is_dictionary(attributes, ProteomeException, 'attributes argument passed to proteome is not a dictionary', or_none = True)

        if attributes is None:            
            self._attributes = {}
        else:
            self._attributes = attributes

        # if no/empty input provided then we're done
        if input_list is None or len(input_list) == 0:
            return

        # else try and add proteins - probably could be more soph
        self.add_proteins(input_list)



    ## ------------------------------------------------------------------------
    ##                
    def check_unique_ID(self, unique_id):
        """
        Function that checks if a given unique ID is found. Note that this 
        function is not needed for testing if a unique_ID is present if the 
        goal is to request Protein Objects (or not). Instead, one can use 
        the .protein(<unique_ID>, safe=False). By setting safe=False if the
        unique_ID is not found then this function will simply return None.

        Parameters
        -----------
        unique_id : string
            String corresponding to a unique_ID associated with some protein

        Returns
        ----------
        bool
            Returns True if the passed ID is present, or False if not.

        """
        if unique_id in self._records.keys():
            return True
        else:
            return False



    ## ------------------------------------------------------------------------
    ##
    @property
    def proteins(self):
        """
        Returns a list of unique_IDs that correspond to the proteins in 
        this Proteome. NOTE this returns a list of the IDs, not the actual 
        Protein objects. To get the corresponding protein object one must 
        use the ``.protein(<unique_ID>)`` notation.
                
        Returns
        --------

        ``list`` of ``str``
            Returns a list of unique IDs

        """

        return list(self._records.keys())



    ## ------------------------------------------------------------------------
    ##
    def protein(self, unique_ID, safe=True):
        """
        Returns the ``Protein`` object associated with the passed unique_ID. 
        If there is no Protein associated with the provided unique_ID then 
        if ``safe=True`` (default) an exception is raised, while if 
        ``safe=False`` then ``None`` is returned.

        Parameters
        -----------
        unique_id : string
            String corresponding to a unique_ID associated with some protein

        safe : bool (default = True)
            If set to True then a missing unique_ID  will raise an exception. 
            If ``False`` then a missing unique_ID will simply return None
            
        Returns
        --------
        
        Protein Object, None
            Depending on if the passed unique_ID is found in the ``Proteome``, 
            a ``Protein`` object or None will be returned.            

        """

        # convert unique IDs to strings because this typing also happens 
        # when new proteins are added
        unique_ID_str = str(unique_ID)

        try:
            return self._records[unique_ID_str]
        except KeyError:
            if safe:
                raise ProteomeException("unique_ID '%s' not found in proteome" % (unique_ID))
            else:
                return None



    ## ------------------------------------------------------------------------
    ##
    def add_protein(self, sequence, name, unique_ID, attributes=None, force_overwrite=False):
        """
        Function that allows the user to add a new protein to a Proteomes 
        in an ad-hoc fashion. In general most of the time it will make 
        sense to add proteins all at once from some input source, but the 
        ability to add proteins one at a time is also useful. 

        If a duplicate unique_ID is passed an exception (ProteomeException) 
        is raised.
        
        Parameters
        -----------
        sequence : string
            Amino acid sequence of the protein. Note  - no sanity check of 
            the sequence is performed.

        name : string
            String reflecting the protein name. Again this can be 
            anything.
        
        unique_id : string
            String corresponding to a unique_ID associated with some 
            protein.

        attributes : dict (default = None)
            The attributes dictionary provides a key-value pairing for 
            arbitrary information. This could include gene names, different 
            types of identifies, protein copy number, a set of protein 
            partners, or anything else one might wish to associated with the
            protein as a whole. Default is None.

        force_overwrite : Bool (default = False)
            If set to False and a unique_ID is included that already is 
            found then this function will raise an exception. However, if 
            set to True it will automatically overwrite the pre-existing 
            entry. (Default = False).
           
        Returns
        --------        
        None
            No return status, but valid proteins included in the input_list 
            will be added to to the underlying proteome.

        """
        
        unique_ID_str = str(unique_ID)

        if unique_ID_str in self._records:
            if force_overwrite is False:
                raise ProteomeException('Non-unique unitque_ID passed [%s]' % (unique_ID_str))

        self._records[unique_ID_str] = Protein(sequence, name, self, unique_ID_str, attributes)



    ## ------------------------------------------------------------------------
    ##        
    def add_proteins(self, input_list, force_overwrite=False):
        r"""
        Function that allows the user to add a multiple new proteins using 
        either a list of protein dictionaries (described below) or a list 
        of Protein objects.        

        **Protein dictionaries**
        
        One mode of adding multiple proteins is by passing a list of 
        Protein dictionaries.

        Protein dictionaries are dictionaries that posses the following 
        key-value pairs ::

            'sequence'   : amino acid sequence (str)
            'name'       : protein name (str)
            'unique_ID'  : The unique identification number used for the 
                           protein (str)
            'attributes' : A dictionary of arbitrary key-value pairs to 
                           associate with the protein (dict or None)
        
        Additional keys/value pairs are ignored and ALL four of these must
        be included. If any are missing for any protein entry this function 
        raises a ProteomeException.

        **Protein objects**
        A second mode of adding multiple proteins is by passing a list of
        Protein objects

        In both cases, the function automatically determines the type of 
        the passed list, and adds dictionaries accordingly. Note that 
        in both cases proteins are added by value - i.e. a new Protein
        object is generated.
        
        Parameters
        -----------
        input_list : list
            List of Protein dictionaries or list of Protein objects
        
        force_overwrite : bool (default = False)
            If set to False and a unique_ID is included that already is 
            found then this function will raise an exception. However, if 
            set to True it will automatically overwrite the pre-existing 
            entry. 

        Returns
        --------
        
        None
            No return status, but valid proteins included in the input_list 
            will be added to to the underlying proteome.
            

        """

        # cycles over every element in the input list and builds a new
        # list where each 
        type_list = list(set([type(i) for i in input_list]))

        # checks if only one type of object is found here
        if len(type_list) > 1:
            raise ProteomeException(f'Trying to add Proteins to a Proteome and the input_list contains more than one type {type_list}')

        # if we're using a list of protein dictionaries
        if type_list[0] == dict:
            self._add_proteins_dict(input_list, force_overwrite)

        # if we're using a list of Protein objects
        elif type_list[0] == Protein:
            self._add_proteins_Protein(input_list, force_overwrite)

        

    ## ------------------------------------------------------------------------
    ##        
    def _add_proteins_dict(self, input_list, force_overwrite=False):
        """
        Internal function that mirrors add_proteins() but operates if every
        element in the input_list is a dictionary.

        Importantly, this works by create a NEW proteins and populating
        based on the key-value mapping in the protein dictionary. 

        Protein dictionaries are dictionaries that posses the following 
        key-value pairs:

            'sequence'   : amino acid sequence (str)
            'name'       : protein name (str)
            'unique_ID'  : The unique identification number used for the 
                           protein (str)
            'attributes' : A dictionary of arbitrary key-value pairs to 
                           associate with the protein (dict or None)


        This function is not be called directly, but instead used by 
        add_proteins() function

        Parameters
        --------------
        input_list : list 
            List of dictionaries  

        force_overwrite : bool, default=False
            Flag which if set to True will mean a Protein in the input_list
            would overwrite an existing protein with the same unique_ID.

        Returns
        --------
        None
            No return type, but will add the protein dictoinaries in the 
            input_list into this Proteome.

        
        """
        
        # for each protein entry in the input list
        for entry in input_list:
        
            try:
                sequence       = str(entry['sequence'])
                name           = str(entry['name'])
                unique_ID      = str(entry['unique_ID'])
                attributes     = entry['attributes']
            except KeyError:
                # if something goes wrong while extracting the four required attributes we build a 
                # diagnosis string and then print this as we raise an exception. The goal here is to
                # try and provide the user with as much info as possible to diagnose the problem
            
                diagnosis_string = self.__build_diagnosis_string_proteome_construction(entry)
                raise ProteomeException('%s'%(diagnosis_string))
            
            if unique_ID in self._records:
                if force_overwrite is False:
                    raise ProteomeException('Non-unique unique_ID passed [%s]' % (unique_ID))

            # add in a new protein
            self._records[unique_ID] = Protein(sequence, name, self, unique_ID, attributes)


    ## ------------------------------------------------------------------------
    ##        
    def _add_proteins_Protein(self, input_list, force_overwrite=False):
        """
        Internal function that mirrors add_proteins() but operates if every
        element in the input_list is a Protein object.

        Importantly, this works by create a NEW protein, and ensure all 
        complex datatypes associated with that new protein are copied. 

        This function is not be called directly, but instead used by 
        add_proteins() function

        Parameters
        --------------
        input_list : list of shephard.protein.Protein objects
            List of Protein objects to be copied into this Proteome

        force_overwrite : bool, default=False
            Flag which if set to True will mean a Protein in the input_list
            would overwrite an existing protein with the same unique_ID.

        Returns
        --------
        None
            No return type, but will add the Protein objects in the input_list
            into this Proteome.
        
        """
        
        # for each protein entry in the input list
        for entry in input_list:
            
            # get the unique ID and, if this ID is already found in the Proteome
            # raise an exception UNLESS force_overwrite is True
            unique_ID = entry.unique_ID
            if unique_ID in self._records:
                if force_overwrite is False:
                    raise ProteomeException('Non-unique unique_ID passed [%s]' % (unique_ID))


            ##
            ## New protein is fully created and all complex data types are copied, so this
            ## new protein is a completely distinct entity to the protein in the original
            ## input list. This avoids any possible issues with cross-referencing back against
            ## old data structures and keeps things clean. Also ensures that the new proteome
            ## has updated unique sites and unique domains lists in the usual way. Basically
            ## this is the most appropritae way to add an existing protein to a new proteome
            ##
            # create new protein
            new_protein = Protein(entry.sequence, entry.name, self, entry.unique_ID)

            # update attributes
            for a in entry.attributes:
                new_protein.add_attribute(a, copy.deepcopy(entry.attribute(a)))

            # update domains
            for d in entry.domains:
                new_protein.add_domain(d.start, d.end, d.domain_type, copy.deepcopy(d._attributes))

            # update sites
            for s in entry.sites:
                new_protein.add_site(s.position. s.site_type, s.symbol, s.value, copy.deepcopy(s._attributes))

            # update tracks
            for t in entry.tracks:
                vals = t.values

                # if vals present then we're creating a values track. Note we can get away with
                # a shallow copy because we know these tracks will only have ints or chars in their
                # elements, which are appropriately copied by a shallow copy
                if vals:                    
                    new_protein.add_track(t.name, vals.copy(), None)
                else:
                    new_protein.add_track(t.name, None, t.symbols.copy())

            self._records[unique_ID] = new_protein
           
            


    ## ------------------------------------------------------------------------
    ##        
    def remove_protein(self, unique_ID, safe=True):
        """
        Function that removes a given protein from the Proteome based on the
        passed unique_ID. If the passed unique_ID does not exist then this 
        will trigger an exception unless safe=False.

        
        Parameters
        ------------
        unique_ID : str
            Unique ID that will be used to retrieve a given protein

        safe : bool (default = True)
            Flag that if set to True means if a passed unique_ID is missing 
            from the underlying proteome object an exception wll be raised 
            (ProteomeException). If set to False, a missing unique_ID is 
            ignored.


        Returns
        -----------
        None
            No return type but will remove an entry from the Proteome.
           
        """

        unique_ID_str = str(unique_ID)

        if unique_ID_str in self._records:
            del self._records[unique_ID_str]
        else:
            if safe:
                raise ProteomeException('Passed unique_ID [%s] not found in this proteome' % (unique_ID_str))
            


    ## ------------------------------------------------------------------------
    ##        
    def remove_proteins(self, input_list, safe=True):
        """
        Function that removes a given proteome from the Proteome based on
        the passed unique_ID. If the passed unique_ID does not exist then 
        this will trigger an exception unless safe = False.

        Parameters
        ------------
        input_list : list of str
            List that contains the unique IDs that will be used to select 
            proteins for deletion.

        safe : bool (default = True)
            Flag that if set to True means if a passed unique_ID is missing 
            from the underlying proteome object an exception wll be raised 
            (ProteomeException). If False a missing unique_ID is ignored.

        Returns
        -----------
        None
            No return type but will remove an entry from the proteome           

        """

        for unique_ID in input_list:
            self.remove_protein(unique_ID, safe=safe)
            


     
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
        Provides a list of the keys associated with every attribute associated
        with this protein.

        Returns
        -------
        list
            returns a list of the attribute keys associated with the Proteome. 


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
             calling the ``<Proteome>.attributes()`` (which returns a list 
             of the valid names).
             
        safe : bool (default = True)
            Flag which if true with throw an exception if an attribute with 
            the same name already exists.
            
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
                raise ProteomeException('Requesting attribute [%s] from Proteome [%s] but this attribute has not been assigned' % (name, str(self))) 

            # if safe not passed just return None
            else:
                return None
                


    ## ------------------------------------------------------------------------
    ##
    def add_attribute(self, name, val, safe=True):
        """
        Function that adds an attribute. Note that if safe is true, 
        this function will raise an exception if the attribute is 
        already present. If safe=False, then an exisiting value will 
        be overwritten.

        Parameters
        ----------------

        name : str
            The parameter name that will be used to identify it

        val : <anything>
            An object or primitive we wish to associate with this 
            attribute

        safe : bool (default = True)
            Flag which if True with throw an exception if an attribute 
            with the same name already exists, otherwise the newly 
            introduced attribute will overwrite the previous one.

        Returns
        ---------
            None - but adds an attribute to the calling object

        """

        if safe:
            if name in self._attributes:
                raise ProteomeException("Trying to add attribute [%s=%s] to Proteome [%s] but this attribute is already set.\nPossible options are: %s" %(name,val, str(self), str(self._attributes.keys())))
                
        self._attributes[name] = val        
        

    ## ------------------------------------------------------------------------
    ##
    def remove_attribute(self, name, safe=True):
        """
        Function that removes a given attribute from the Proteome based on the 
        passed attribute name. If the passed attribute does not exist or is not 
        associate with the protein then this will trigger an exception 
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
            Proteome if present.
            
        """

        if name not in self._attributes:
            if safe:
                raise ProteomeException(f'Passed attribute [{name}] not found in {self}')
        else:
            del self._attributes[name]



    ###################################
    ##                               ##
    ##       DOMAIN FUNCTIONS        ##
    ##                               ##
    ###################################

    ## ------------------------------------------------------------------------
    ##
    @property
    def domains(self):
        """
        Function that returns a list of all domain objects associated with
        the Proteome. 

        This function is useful if you wish to indiscriminately ask questions 
        of domains without considering the proteins they come from. However, 
        each Domain has a Protein object associated with it (via the .protein 
        operator), so one can always map a Domain back to a Protein.

        Returns
        --------------
        list of Domains
             A list of all the Domains from every protein in the Proteome
        
        """

        all_domains = []
        for prot in self:
            all_domains.extend(prot.domains)

        return all_domains


        
    ## ------------------------------------------------------------------------
    ##                
    @property
    def unique_domain_types(self):
        """
        Returns the list of unique Domain types associated with this Proteome.

        Return
        -------

        list of str
            Each element in the list is a string that corresponds to a Domain 
            type.

        """

        # Some description of what's going on here is in order. Every time a new domain
        # is added, the Domain constructor calls the function _Domain__update_domain_types
        # which checks if the domain_type of the domain being added is already in the
        # _unique_domain_types list. If yes, fine, if no, it gets added. This means
        # _unique_domain_types keeps track of a count of the complete number of unique
        # domains in the Proteome. An analogous setup holds true for the sites and tracks.
        
        return list(self._unique_domain_types.keys())


    ## ------------------------------------------------------------------------
    ##
    def get_domains_by_type(self, domain_type, perfect_match=True):
        """
        Function that returns a list of domains from all proteins that matched against
        a specific domain type name. 
        
        Parameters
        ------------
        domain_type : string
            String associated domain_type that you want to search for.

        perfect_match : bool (default = True)
            Flag that identifies if the domain names should be a perfect 
            match (=True) or if the string passed should just appear 
            somewhere in the domain_type string
            
        Returns
        -----------
        list
            Returns a list of Domain objects that match the requested type. 
            Objects are ordered by starting position in sequence.
                                
        """
        return_list = []
        for p in self:
            return_list.extend(p.get_domains_by_type(domain_type, perfect_match))

        return return_list



    ###################################
    ##                               ##
    ##        SITES FUNCTIONS        ##
    ##                               ##
    ###################################

    ## ------------------------------------------------------------------------
    ##                
    @property
    def sites(self):
        """
        Function that returns a list of all Site objects associated with 
        the Proteome. 

        This function is useful if you wish to indiscriminately ask questions 
        of sites without considering the proteins they come from. However, 
        each Site has a Protein object associated with it (via .protein 
        operator), so one can always map a Site back to a Protein.

        Returns
        --------------
        list of Sites
             A list of all the Sites from every protein in the Proteome
        
        """

        all_sites = []

        for prot in self:
            all_sites.extend(prot.sites)
        return all_sites
                    

    ## ------------------------------------------------------------------------
    ##                
    @property
    def unique_site_types(self):
        """
        Returns the list of unique Site types associated with this Proteome.

        Return
        -------

        list of str
            Each element in the list is a string that corresponds to a Site type

        """

        # Some description of what's going on here is in order. Every time a new site
        # is added, the Site constructor calls the function _Site__update_site_types
        # which checks if the domain_type of the domain being added is already in the
        # _unique_site_types list. If yes, fine, if no, it gets added. This means
        # _unique_site_types keeps track of a count of the complete number of unique
        # sites in the Proteome. An analogous setup holds true for domains and tracks.

        return list(self._unique_site_types.keys())
    

    ## ------------------------------------------------------------------------
    ##
    def get_sites_by_type(self, site_types):
        """
        Function that returns a list of sites from all proteins that matched against
        a specific site type name or set of site type names.
        
        Parameters
        ------------

        site_types : string or list of strings
            One or more possible site_types that may be found in the protein. 
            Either a single string or a list of strings can be passed, 
            allowing for one or more sites to be grouped together

            
        Returns
        -----------
        list
            Returns a list of Domain objects that match the requested type. 
            Objects are ordered by starting position in sequence.
                                
        """
        return_list = []
        for p in self:
            return_list.extend(p.get_sites_by_type(site_types, return_list=True))

        return return_list



    ###################################
    ##                               ##
    ##        TRACK FUNCTIONS        ##
    ##                               ##
    ###################################

    ## ------------------------------------------------------------------------
    ##                
    @property
    def unique_track_names(self):
        """
        Returns the list of unique Track names associated with this Proteome.

        Return
        -------

        list of strings
            Each element in the list is a string that corresponds to a 
            Track name found in one (or more) proteins

        """

        # Some description of what's going on here is in order. Every time a new track
        # is added, the Track constructor calls the function _Track__update_track_names
        # which checks if the track_name of the domain being added is already in the
        # _unique_track_types list. If yes, fine, if no, it gets added. This means
        # _unique_track_names keeps track of a count of the complete number of unique
        # domains in the Proteome. An analogous setup holds true for sites and domains
        
        
        return list(self._unique_track_names.keys())

    ## ------------------------------------------------------------------------
    ##                
    @property
    def track_names_to_track_type(self):
        """
        Returns a (copy of a) dictionry that maps track name to track type. 
        We return a copy so there's no way we can accidentally break the 
        internal book-keeping of the Proteome object.

        Return
        -------
        dict
            A dictionary that contains the unique track names and maps 
            each name to either values or type.
        
        """

        return dict(self._track_name_to_track_type)


    ####################################
    ##                                ##
    ##       INTERNAL FUNCTIONS       ##
    ##                                ##
    ####################################

    ## ------------------------------------------------------------------------
    ##                
    def __len__(self):
        """
        The length of the Proteome is defined as the number of proteins 
        in it.
        
        Returns
        -------
        int
            Returns an integer that reflects the number of proteins
        """

        # this function means when we call len(Proteome_Object) we get back
        # the number of proteins in it

        return len(self._records)



    ## ------------------------------------------------------------------------
    ##                
    def __repr__(self):
        """
        Provides a nice representation of the Proteome
        
        Returns
        -------
        string
            Formatted description of the Proteome.

        """

        # this function means when we print a Proteome object or cast it to 
        # a string we get a nice/informative representation, rather than the
        # id of the object

        return "[Proteome]: Sequence dataset with %i protein records" %(len(self))



    ## ------------------------------------------------------------------------
    ##                
    def __iter__(self):
        """
        Allows a Proteome object to act as a generator that yields actual 
        proteins, so the syntax
        
        .. code-block:: python

           for protein in ProteomeObject:
               print(protein.sequence)
        
        is be valid and would iterate through the proteins in the Proteome. 
        
        This makes performing some analysis over all proteins quite easy.

        """

        for i in self._records:
            yield self._records[i]

    ## ------------------------------------------------------------------------
    ##                
    def __contains__(self, m):
        """
        Enables the syntax X in Proteome to be used, where X can be 
        either a unique ID or a Proteome object.

        .. code-block:: python

           if protein.unique_ID in ProteomeObject:
               print(f'The protein {protein} is in the Proteome!')


        """

        if type(m) == str:
            if m in self._records.keys():
                return True
            else:
                return False

        elif type(m) == Protein:
            if m.unique_ID in self._records.keys():
                return True
            else:
                return False

    ## ------------------------------------------------------------------------
    ##                        
    def __getitem__(self, key):
        """
        Allows slicing index into Proteome to retrieve subsets of protein

        .. code-block:: python

           first_protein = ProteomeObject[0]
           print(f'The first protein is {first_protein}')

        """


        if isinstance(key, int) and key >= 0:
            return list(islice([self._records[i] for i in self._records], key, key+1))[0]

        elif isinstance(key, slice):
            return list(islice([self._records[i] for i in self._records], key.start, key.stop, key.step))
        else:
            raise KeyError("Key must be non-negative integer or slice, not {}"
                           .format(key))


    ## ------------------------------------------------------------------------
    ##                
    def _Domain__update_domain_types(self, domain_type):
        """
        INTERNAL FUNCTION (not for public API use)

        Note - we this function is named as __Domain_... so it can be 
        specifically and uniquely be called from a Domain object. This 
        function is ONLY called last thing in the Domain constructor 
        where it allows the Proteome object to keep track of the total 
        number of unique domain types in the Proteome.

        The function is (by default) called by the Domain constructor.
        
        Parameters
        ----------------
        domain_type : string
            String that defines a domain type

        Returns
        ---------------

        No return value, but will appropriately update the Proteome object

        """
        # because _unique_track_names is a dictionary this scales O(1) with number 
        # of track names
        if domain_type not in self._unique_domain_types:
            self._unique_domain_types[domain_type] = 1
        else:
            self._unique_domain_types[domain_type] = self._unique_domain_types[domain_type] + 1

    ## ------------------------------------------------------------------------
    ##                
    def _Protein__decrement_domain_types(self, domain_type):
        """
        INTERNAL FUNCTION (not for public API use)
        
        Note - we this function is named as __Protein_... so it 
        can be specifically and uniquely be called from a Protein 
        object. This function is ONLY called last thing when a 
        Protein object deletes a domain

        Parameters
        ----------------
        track_name : string
            String that defines a Domain type

        Returns
        ---------------

        No return value, but will appropriately update the Proteome 
        object

        """

        # if we can't find the domain name in the unique domain types... this is bad!
        if domain_type not in self._unique_domain_types:
            raise ProteomeException("Tried to remove a Domain type [{domain_type}] from the Proteome.unique_domain_types dictionary but the Domain type could not be found. This is a bug. Please report as a GitHub Issue.")
            
        # we are removing a unique domain_type name! Big event!
        elif self._unique_domain_types[domain_type] == 1:
            del self._unique_domain_types[domain_type]
        # else decrement one
        else:
            self._unique_domain_types[domain_type] = self._unique_domain_types[domain_type] - 1


    ## ------------------------------------------------------------------------
    ##                
    def _Track__update_track_names(self, track_name, track_type):
        """
        INTERNAL FUNCTION (not for public API use)
        
        Note - we this function is named as __Track_... so it can be 
        specifically and uniquely be called from a Track object. This 
        function is ONLY called last thing in the Track constructor 
        where it allows the Proteome object to keep track of the total 
        number of unique Track types in the Proteome.

        The function is (by default) called by the Track constructor

        Parameters
        ----------------
        track_name : string
            String that defines the Track name 

        Returns
        ---------------

        No return value, but will appropriately update the Proteome object

        """

        # because _unique_track_names is a dictionary this scales O(1) with number 
        # of track names
        if track_name not in self._unique_track_names:
            self._unique_track_names[track_name] = 1
            self._track_name_to_track_type[track_name] = track_type

        else:
            self._unique_track_names[track_name] = self._unique_track_names[track_name] + 1
            if self._track_name_to_track_type[track_name] != track_type:
                raise ProteomeException(f"Tried to assigned track name [{track_name}] as a [{track_type}] track, but this track was already assigned as a [{self._track_name_to_track_type[track_name]}] track. Cannot have two tracks with the same name but different types")
            

    ## ------------------------------------------------------------------------
    ##                
    def _Protein__decrement_track_names(self, track_name):
        """
        INTERNAL FUNCTION (not for public API use)
        
        Note - we this function is named as __Protein_... so it can be
        specifically and uniquely be called from a Protein object. 
        This function is ONLY called last thing when a Protein object 
        deletes a Track
                

        Parameters
        ----------------
        track_name : string
            String that defines the Track name 

        Returns
        ---------------

        No return value, but will appropriately update the Proteome 
        object.

        """

        # if we can't find the track name in the unique track names... this is bad!
        if track_name not in self._unique_track_names:
            raise ProteomeException("Tried to remove a Track name [{track_name}] from the Proteome.unique_track_names dictionary but the track could not be found. This is a bug. Please report as a GitHub Issue.")
            
        # we are removing a unique track name! Big event!
        elif self._unique_track_names[track_name] == 1:
            del self._unique_track_names[track_name]
            del self._track_name_to_track_type[track_name]

        # else decrement one
        else:
            self._unique_track_names[track_name] = self._unique_track_names[track_name] - 1


    ## ------------------------------------------------------------------------
    ##                
    def _Site__update_site_types(self, site_type):
        """
        INTERNAL FUNCTION (not for public API use).

        Note - we this function is named as __SITE_... so it can be 
        specifically and uniquely be called from a Site object. This 
        function is ONLY called last thing in the Site constructor where
        it allows the Proteome object to keep track of the total number 
        of unique Site types in the Proteome.

       
        The function is (by default) called by the Site constructor.

        Parameters
        ----------------
        site_type : string
            String that defines a site type

        Returns
        ---------------

        No return value, but will appropriately update the Proteome object.

        """

        # because _unique_track_names is a dictionary this scales O(1) with number 
        # of track names
        if site_type not in self._unique_site_types:
            self._unique_site_types[site_type] = 1

        else:
            self._unique_site_types[site_type] = self._unique_site_types[site_type] + 1


    ## ------------------------------------------------------------------------
    ##                
    def _Protein__decrement_site_types(self, site_type):
        """
        INTERNAL FUNCTION (not for public API use)

        
        Note - we this function is named as __Protein_... so it can be specifically
        and uniquely be called from a Protein object. This function is ONLY called
        last thing when a Protein object deletes a Site

        Parameters
        ----------------
        track_name : string
            String that defines a site type

        Returns
        ---------------

        No return value, but will appropriately update the Proteome object

        """

        # if we can't find the track name in the unique track names... this is bad!
        if site_type not in self._unique_site_types:
            raise ProteomeException("Tried to remove a Site type [{site_type}] from the Proteome.unique_site_types dictionary but the Site type could not be found. This is a bug. Please report as a GitHub Issue.")
            
        # we are removing a unique site_type name! Big event!
        elif self._unique_site_types[site_type] == 1:
            del self._unique_site_types[site_type]
        # else decrement one
        else:
            self._unique_site_types[site_type] = self._unique_site_types[site_type] - 1



    ## ------------------------------------------------------------------------
    ##                
    def __build_diagnosis_string_proteome_construction(self, entry):
        """
        INTERNAL FUNCTION (not for public API use)
        
        Function that builds a string to help in diagnosing what might go 
        wrong during Proteome construction. Called by the Proteome constructor.
        
        
        """

        ds = "Error building proteome when parsing the following entry:\n"


        # check the sequence
        try:
            s = str(entry['sequence'])
        except Exception:
            s = 'FAILED'

        ds = ds + "sequence: %s\n" %(s)


        # check the name
        try:
            s = str(entry['name'])
        except Exception:
            s = 'FAILED'

        ds = ds +"name: %s\n" %(s)


        # check the unique_ID
        try:
            s = str(entry['unique_ID'])
        except Exception:
            s = 'FAILED'

        ds = ds +"unique_ID: %s\n" %(s)


        # 
        try:
            s = str(entry['attributes'])
        except Exception:
            s = 'FAILED'

        ds = ds +"attributes: %s\n" %(s)

        return ds

