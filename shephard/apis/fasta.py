"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

import protfasta
from shephard.proteome import Proteome

## ------------------------------------------------------------------------
##
def fasta_to_proteome(filename, build_unique_ID=None, invalid_sequence_action='fail'):
    """
    Stand alone function that allows the user to build a Proteome from a standard
    FASTA file. 

    The input filename must be a FASTA file without duplicate headers. If the file
    has duplicate headers and these have to be further processed we suggest using
    the protfasta (https://protfasta.readthedocs.io/) package to parse through the
    FASTA file first creating a santizied input FASTA.
    
    Each protein in a Proteome must have a unique_ID associated with it. There
    are two ways a FASTA file can be used to generate a unique ID:

        1. By parsing the FASTA header, which could be as much as simply reading 
           the header or couple involve some more complex logic.

        2. By incrementing an automatically unique ID.

    IF the argument ``build_unique_ID`` is not provided, the ``fasta_to_proteome`` function
    will automatically generate a unique numerical ID for each protein.

    However, if the ``build_unique_ID`` argument *is* provided, this function is used to
    convert the header into a unique key.
    
    Parameters
    ------------

    filename : string
        Name of the FASTA file we're going to parse in. Note the protein name will be
        defined as the full FASTA header for each entry **unless** a ``header_parser``
        function is provided

    build_unique_ID : function
        [**Default = None**] ``build_unique_ID`` allows a user-defined function that is 
        used to convert the FASTA header to a (hopefully) unique string. This can be 
        useful if the FASTA header is well structured and includes a specific, useful
        unique string that can be used as the unique_ID.

    invalid_sequence_action : ``'ignore'``, ``'fail'``, ``'remove'``, ``'convert'``, ``'convert-ignore'``
        [**Default = 'fail'**] Selector that determines how to deal with invalid sequences. If ``convert``
        or ``convert-ignore`` are chosen, then conversion is completed with either the standard conversion 
        table (shown under the ``correction_dictionary`` documentation) or with a custom conversion dictionary 
        passed to ``correction_dictionary``. 
        Options are as follows: 
            * ``ignore``  - invalid sequences are completely ignored

            * ``fail``    - invalid sequence cause parsing to fail and throw an exception
  
            * ``remove`` - invalid sequences are removed

            * ``convert`` - invalid sequences are convert

            * ``convert-ignore`` - invalid sequences are converted to valid sequences and any remaining invalid residues are ignored
        
    Returns 
    --------
    Proteome
        Returns an initialized Proteome object 
    
    """

    # read in the fasta file
    fasta_dictionary = protfasta.read_fasta(filename, invalid_sequence_action=invalid_sequence_action)

    # extract the keys (FASTA headers) and initialize the record_index (internal
    # numbering used for construction. Also initialize the proteom_dict, which is
    # a dictionary of protein entries we passed to Proteome.
    record_index  = 0
    proteome_list = []

    # for each entry
    for k in fasta_dictionary:

        # create a key-value pair where 
        #   key = the unique record_index (this is only used for internal structure
        #         within this function to assure we never overwrite in this dictionary
        #
        #  value = a four-position list where the positions reflect the following
        #        [0] = amino acid sequence
        #        [1] = name (this can be anything)
        #        [2] = unique_ID - this should be a unique identifier that can be used
        #              to cross-reference this entry to other data. If extrat_unique_ID
        #              is passed we try to use this 
        #        [3] = attribute dictionary (we set this to None)

        if build_unique_ID:
            unique_ID = build_unique_ID(k)
        else:
            unique_ID = record_index
            
        newdict = {}
        newdict['sequence'] = str(fasta_dictionary[k])
        newdict['name'] = k
        newdict['unique_ID'] = unique_ID
        newdict['attribute_dictionary'] = None

        proteome_list.append(newdict)

        record_index = record_index + 1

    # finally actually build the proteome and return it
    return Proteome(proteome_list)


## ------------------------------------------------------------------------
##
def proteome_to_fasta(filename, proteome, include_unique_ID_in_header=False):
    """
    Stand alone function that allows the user to write a FASTA file from a
    Proteome file. 


    Parameters
    ------------
    filename : string
        Name of the FASTA file we're going to write to. We will automatically overwrite 
        a file if it's there, so be careful!

    proteome : Proteome object
        The proteome object we wish to write to disk.

    include_unique_ID_in_header : bool (default False)
        Sometimes it may be desirable to 

    """

    # build output list with or without the unique_ID 
    outlist = []
    for protein in proteome:

        # this is where we define the FASTA header...
        if include_unique_ID_in_header:
            fasta_header = "%s | UID=%s" %(protein.name, protein.unique_ID)
        else:
            fasta_header = protein.name
            
        outlist.append([fasta_header, protein.sequence])
        
    # use the protfasta library to write the file to disk
    protfasta.write_fasta_file(outlist, filename, linelength=80)
                


