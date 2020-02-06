"""
quickstart.py

From the SHEPHARD package
Sequence-based Hierachical and Extendable Platform for High-throughput Analysis of Region of Disorder
Ginell & Holehouse, 2020

Handles the primary functions
"""


from pyfaidx import Fasta
from .proteome import Proteome
from .exceptions import UtilitiesException


## ------------------------------------------------------------------------
##
def quickstart(filename, extract_unique_ID=None):
    """
    Quickstart provides an easy way to create an initial Proteome object by passin in a FASTA
    file. 

    """

    # read in the fasta file. In the future we will probably need to write out own more
    # robust FASTA parser, but for now this will do
    # IN = Fasta(filename, duplicate_action='longest', read_long_names=True)
    IN = Fasta(filename, read_long_names=True)

    # extract the keys (FASTA headers) and initialize the record_index (internal
    # numbering used for construction. Also initialize the proteom_dict, which is
    # a dictionary of protein entries we passe to Proteome.
    KEYZ = IN.keys()        
    record_index  = 0
    proteome_list = []

    # for each entry
    for k in KEYZ:

        # create a key-value pair where 
        #   key = the unique record_index (this is only used for internal structure
        #         within this function to assure we never overwright in this dictionary
        #
        #  value = a four-position list where the positions reflect the following
        #        [0] = amino acid sequence
        #        [1] = name (this can be anything)
        #        [2] = unique_ID - this should be a unique identifier that can be used
        #              to cross-reference this entry to other data. If extrat_unique_ID
        #              is passed we try to use this 
        #        [3] = attribute dictionary (we set this to None)

        if extract_unique_ID:
            unique_ID = extract_unique_ID(k)
        else:
            unique_ID = record_index
            
        newdict = {}
        newdict['sequence'] = str(IN[k][:])
        newdict['name'] = k
        newdict['unique_ID'] = unique_ID
        newdict['attribute_dictionary'] = None

        proteome_list.append(newdict)

        record_index = record_index + 1

    S = Proteome(proteome_list)

    return S



