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
def quickstart(filename, extract_unique_id=None):
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
    proteome_dict = {}

    # for each entry
    for k in KEYZ:

        # create a key-value pair where 
        #   key = the unique record_index (this is only used for internal structure
        #         within this function to assure we never overwright in this dictionary
        #
        #  value = a four-position list where the positions reflect the following
        #        [0] = amino acid sequence
        #        [1] = name (this can be anything)
        #        [2] = unique_id - this should be a unique identifier that can be used
        #              to cross-reference this entry to other data. If extrat_unique_id
        #              is passed we try to use this 
        #        [3] = attribute dictionary (we set this to None)

        if extract_unique_id:
            unique_id = extract_unique_id(k)
        else:
            unique_id = None


        proteome_dict[record_index] = [str(IN[k][:]), k, unique_id, None]
        record_index = record_index + 1

    S = Proteome(proteome_dict)

    return S



