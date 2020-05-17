"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

from .exceptions import ProteomeException
from .protein import Protein
from itertools import islice

    

class Proteome:
    """
    The Proteome object is the main unit for information storage in **SHEPHARD**.

    The Proteome constructor takes a single argument, which is a list of protein dictionaries. This means
    Proteome objects can be generated directly (see below for a definition of protein dictionaries). However,
    it is often more convenient to build Proteomes from FASTA files. For more information on this see the ``api``
    documentation.
    
    Protein dictionaries are dictionaries that must contain four elements (others are ignored).
    
    ``sequence`` : *str* - Amino acid sequence of the protein
        
    ``name`` : *str* - Name of the protein (this can be anything, it is not used internally so no constraints on what this is 
                                      
    ``unique_ID`` :  *str* - This must be unique with respect to all other unique_IDs in the set of proteins in the input list. 
                                      
    ``attribute_dict`` : *dict* - Dictionary of one or more attributes to apply to this protein. Key/value pairs in this dictionary can be arbitrary and are user defined.

    As an example

    >>> protein_dictionary_example = {'sequence':'ALAPSLLPAMPALSPALSP', 'name': 'my protein fragment', 'unique_ID':'UXX01', 'attribute_dict':{}}
    >>> dictionary_list = []
    >>> dictionary_list.append(protein_dictionary_example)
    >>> P = Proteome(dictionary_list)
                                          
    Note that ``sequence``, ``name`` and ``unique_ID`` are cast to *str* by the function, so if numerical values are passed for any these will be converted to strings.
    
    **Notes**

    * NOTE that ALL FOUR of these are required for EACH protein, even if the attribute_dict is empty.
    * The unique_ID is checked for uniqueness against all others in the Proteomes and will throw and exception if it is, in fact, not unique.
    * Additional proteins can be added using the `.add_protein()` or `.add_proteins() function.

    """

    ## ------------------------------------------------------------------------
    ##
    def __init__(self, input_list, attribute_dictionary = None):
        # See the Proteome class documentation for constructor info
        """
        """

        # initiallize book-keeping instruments
        self._records = {}
        self._unique_domain_types = []
        self._unique_site_types = []
        
        # set attribute dictionary IF a dictionary was passed
        if isinstance(attribute_dictionary, dict):
            self._attributes = attribute_dictionary

        # set dictionary to an empty dictionary if none was passed
        elif attribute_dictionary is None:
            self._attributes = {}

        else:
            raise exceptions.ProteinException('[FATAL]: If provided, protein attribute must a dictionary')

        # for each entry in the input list
        for entry in input_list:
        
            try:
                sequence       = str(entry['sequence'])
                name           = str(entry['name'])
                unique_ID      = str(entry['unique_ID'])
                attribute_dict = entry['attribute_dictionary']
            except KeyError:
                # if something goes wrong while extracting the four required attributes we build a 
                # diagnosis string and then print this as we raise an exception. The goal here is to
                # try and provide the user with as much info as possible to diagnose the problem
            
                diagnosis_string = __build_diagnosis_string_proteome_construction(entry)
                raise ProteomeException('%s'%(diagnosis_string))
            
            if unique_ID in self._records:
                raise ProteomeException('Non-unique unique_ID passed [%s]' % (unique_ID))

            # add in a new protein
            self._records[unique_ID] = Protein(sequence, name, self, unique_ID, attribute_dict)



    ## ------------------------------------------------------------------------
    ##
    @property
    def proteins(self):
        """
        Returns a list of unique_IDs that correspond to the proteins in this Proteome.
        NOTE this returns a list of the IDs, not the actual Protein objects. To get the
        corresponding protein object one must use the ``.protein(<unique_ID>)`` notation.

        Returns
        --------

        ``list`` of ``str``
            Returns a list of unique_IDs

        """

        return list(self._records.keys())



    ## ------------------------------------------------------------------------
    ##
    def protein(self, unique_ID, safe=True):
        """
        Returns the ``Protein`` object associated with the passed unique_ID. If there is no
        Protein associated with the provided unique_ID then if ``safe=True`` (default) an 
        exception is raised, while if ``safe=False`` then ``None`` is returned.

        Parameters
        -----------
        unique_id : string
            String corresponding to a unique_ID associated with some protein

        safe : bool 
            [**Default = ``True``]** If set to True then a missing unique_ID will raise an exception. If ``False``
            then a missing unique_ID will simply return None

        Returns
        --------
        
        Protein Object or None
            Depending on if the passed unique_ID is found in the ``Proteome``, a ``Protein`` 
            object or None will be returned

        """

        try:
            return self._records[unique_ID]
        except KeyError:
            if safe:
                raise ProteomeException('unique_ID not found in proteome')
            else:
                return None



    ## ------------------------------------------------------------------------
    ##
    def add_protein(self, sequence, name, unique_ID, attribute_dictionary, safe=True):
        """
        Function that allows the user to add a new protein to a Proteomes in an 
        ad-hoc fashion. In general most of the time it will make sense to add
        proteins all at once from some input source, but the ability to add
        proteins one at a time is also useful
        
        Parameters
        -----------
        sequence : string
            Amino acid sequence of the protein. Note  - no sanity check of the
            sequence is performed.

        name : string
            String reflecting the protein name. Again this can be anything
        
        unique_id : string
            String corresponding to a unique_ID associated with some protein

        attribute_dictionary : dict
            The attribute_dictionary provides a key-value pairing for arbitrary information.
            This could include gene names, different types of identifies, protein copy number,
            a set of protein partners, or anything else one might wish to associated with the
            protein as a whole.

        safe : boolean (default = True)
            If set to True then a duplicate unique_ID will raise an exception. If false
            then a duplicate unique_ID will simply return None

        Returns
        --------
        
        Protein object or None
            Depending on if the passed unique_ID is found in the Proteomes, a Protein 
            object or None will be returned

        """

        if unique_ID in self._records:
            if safe:
                raise ProteomeException('Non-unique unique_ID passed [%s]' % (unique_ID))
            else:
                return

        self._records[unique_ID] = Protein(sequence, name, self, unique_ID, attribute_dictionary)
        

     
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
             ``<Proteome>.attributes()`` (which returns a list of the valid names)

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
                raise ProteomeException('Requesting attribute [%s] from Proteome [%s] but this attribute has not been assigned' % (name, str(self))) 

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
                raise ProteomeException("Trying to add attribute [%s=%s] to Proteome [%s] but this attribute is already set.\nPossible options are: %s" %(name,val, str(self), str(self._attributes.keys())))
                
        self._attributes[name] = val        
        

    ## ------------------------------------------------------------------------
    ##
    @property
    def domains(self):
        """
        Function that returns a list of all domain objects associated with the Proteome. Unlike
        when .domains is called on a Protein object (where a list of domain_IDs are returned)
        the complete list of Domain objects is returned in a list. 

        This function is useful if you wish to indiscriminately ask questions of domains without
        considering the proteins they come from. However, each Domain has a .protein object 
        associated with it, so one can always map a Domain back to a Protein.

        Returns
        --------------
        list of Domains
             A list of all the Domains from every protein in the Proteome
        
        """

        all_domains = []

        for prot in self:
            for domain_name in prot.domains:
                all_domains.append(prot.domain(domain_name))
        return all_domains



    ## ------------------------------------------------------------------------
    ##                
    @property
    def sites(self):
        """
        Function that returns a list of all Site objects associated with the Proteome. Unlike
        when .sites is called on a Protein object (where a list of domain_IDs are returned)
        the complete list of Domain objects is returned in a list. 

        This function is useful if you wish to indiscriminately ask questions of sites without
        considering the proteins they come from. However, each Site has a .protein object 
        associated with it, so one can always map a Site back to a Protein.        

        Returns
        --------------
        list of Domains
             A list of all the Domains from every protein in the Proteome
        
        """

        all_sites = []

        for prot in self:
            for site_name in prot.sites:
                sitelist = prot.site(site_name)
                for s in sitelist:
                    all_sites.append(s)
        return all_sites
                

        
    ## ------------------------------------------------------------------------
    ##                
    @property
    def unique_domain_types(self):
        """
        Returns the list of unique Domain types associated with this Proteome.

        Return
        -------

        list of strings
            Each element in the list is a string that corresponds to a Domain type

        """

        # Some description of what's going on here is in order. Every time a new domain
        # is added, the Domain constructor calls the function _Domain__update_domain_types
        # which checks if the domain_type of the domain being added is already in the
        # __unique_domain_types list. If yes, fine, if no, it gets added. This means
        # _unique_domain_types keeps track of a count of the complete number of unique
        # domains in the Proteome. An analogous setup holds true for the sites.        
        
        return self._unique_domain_types
    


    ## ------------------------------------------------------------------------
    ##                
    @property
    def unique_site_types(self):
        """
        Returns the list of unique Site types associated with this Proteome.

        Return
        -------

        list of strings
            Each element in the list is a string that corresponds to a Site type

        """

        # Some description of what's going on here is in order. Every time a new site
        # is added, the Site constructor calls the function _Site__update_site_types
        # which checks if the domain_type of the domain being added is already in the
        # __unique_site_types list. If yes, fine, if no, it gets added. This means
        # _unique_site_types keeps track of a count of the complete number of unique
        # sites in the Proteome. An analogous setup holds true for the domains.        

        return self._unique_site_types
    


    ## ------------------------------------------------------------------------
    ##                
    def check_unique_ID(self, unique_id):
        """
        Function that checks if a given unique ID is found. Note that this function is not needed
        for testing if a unique_ID is present if the goal is to request Protein Objects (or not).
        Instead, one can use the .protein(<unique_ID>, safe=False). By setting safe=False if the
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


    ####################################
    ##                                ##
    ##       INTERNAL FUNCTIONS       ##
    ##                                ##
    ####################################

    ## ------------------------------------------------------------------------
    ##                
    def __len__(self):
        """
        The length of the Proteome is defined as the number of proteins in it
        
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
        Allows a Proteome object to act as a generator that yields actual proteins,
        so the syntax

        >>> for protein in ProteomeObject:
        >>>    print(protein.sequence)

        is be valid and would iterate through the proteins in the Proteome. 
        
        This makes performing some analysis over all proteins quite easy        

        """

        for i in self._records:
            yield self._records[i]

    ## ------------------------------------------------------------------------
    ##                
    def __contains__(self, m):
        if m in self._records.keys():
            return True
        else:
            return False

    ## ------------------------------------------------------------------------
    ##                        
    def __getitem__(self, key):
        """
        Allows slicing index into Proteome to retrieve subsets of protein
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

        
        Note - we this function is named as __Domain_... so it can be specifically
        and uniquely be called from a Domain object. This function is ONLY called
        last thing in the Domain constructor where it allows the Proteome object
        to keep track of the total number of unique domain types in the Proteome.

        The function is (by default) called by the Domain constructor

        Parameters
        ----------------
        domain_type : string
            String that defines a domain type

        Returns
        ---------------

        No return value, but will appropriately update the Proteome object

        """
        if domain_type not in self.unique_domain_types:
            self._unique_domain_types.append(domain_type)



    ## ------------------------------------------------------------------------
    ##                
    def _Site__update_site_types(self, site_type):
        """
        INTERNAL FUNCTION (not for public API use)


        Note - we this function is named as __SITE_... so it can be specifically
        and uniquely be called from a Site object. This function is ONLY called
        last thing in the Site constructor where it allows the Proteome object
        to keep track of the total number of unique domain types in the Proteome

        The function is (by default) called by the Site constructor

        Parameters
        ----------------
        site_type : string
            String that defines a site type

        Returns
        ---------------

        No return value, but will appropriately update the Proteome object

        """
        if site_type not in self.unique_site_types:
            self._unique_site_types.append(site_type)



    ## ------------------------------------------------------------------------
    ##                
    def __build_diagnosis_string_proteome_construction(self, entry):
        """
        INTERNAL FUNCTION (not for public API use)
        
        Function that builds a string to help in diagnosing what might go wrong during
        Proteome construction. Called by the Proteome constructor.        
        
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
            s = str(entry['attribute_dictionary'])
        except Exception:
            s = 'FAILED'

        ds = ds +"attribute_dictionary: %s\n" %(s)

        return ds
                                        


    
        
    
