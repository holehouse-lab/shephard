"""
proteome.py

From the SHEPHARD package
Sequence-based Hierachical and Extendable Platform for High-throughput Analysis of Region of Disorder
Ginell & Holehouse, 2020

Handles the primary functions
"""

from .exceptions import ProteomeException
from .protein import Protein
    

class Proteome:

    def __init__(self, input_dictionary):
        """
        To create 

        """

        # initiallize book-keeping instruments
        self._records = {}
        self._unique_domain_types = []
        self._unique_site_types = []

        for i in input_dictionary:

            sequence   = input_dictionary[i][0]
            name       = input_dictionary[i][1]
            unique_ID  = input_dictionary[i][2]
            attribute_dict = input_dictionary[i][3]

            if unique_ID in self._records:
                raise ProteomeException('Non-unique unique_ID passed [%s]' % (unique_ID))

            # add in a new protein
            self._records[unique_ID] = Protein(sequence, name, self, unique_ID, attribute_dict)

    @property
    def proteins(self):
        return list(self._records.keys())

    def protein(self, unique_ID):
        return self._records[unique_ID]

    @property
    def unique_domain_types(self):
        return self._unique_domain_types
    
    @property
    def unique_site_types(self):
        return self._unique_site_types
    

    def __len__(self):
        return len(self._records)

    def __repr__(self):
        return "[Proteome]: Sequence dataset with %i protein records" %(len(self))

    def __iter__(self):
        for i in self._records:
            yield self._records[i]

    def _Domain__update_domain_types(self, domain_type):
        """
        Note - we this function is named as __Domain_... so it can be specifically
        and uniquely be called from a Domain object. This function is ONLY called
        last thing in the Domain constructor where it allows the proteome object
        to keep track of the total number of unique domain types in the proteome

        """
        if domain_type not in self.unique_domain_types:
            self._unique_domain_types.append(domain_type)

    def _Site__update_site_types(self, site_type):
        """
        Note - we this function is named as __SITE_... so it can be specifically
        and uniquely be called from a Site object. This function is ONLY called
        last thing in the Site constructor where it allows the proteome object
        to keep track of the total number of unique domain types in the proteome

        """
        if site_type not in self.unique_site_types:
            self._unique_site_types.append(site_type)

            
    def get_number_of_proteins(self):
        return len(self._records)


    def unique_ID_present(self, unique_id):
        """
        Checks if a given unique ID is inside 

        """
        if unique_id in self._records.keys():
            return True
        else:
            return False


    
        
    
