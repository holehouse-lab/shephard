"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

import protfasta
from shephard.proteome import Proteome
from shephard.exceptions import APIException

SHEPHARD_ATTRIBUTE_SPLITTER='SHPRD_ATTRIBUTES='

## ------------------------------------------------------------------------
##
def fasta_to_proteome(filename, 
                      proteome=None, 
                      build_unique_ID=None, 
                      build_attributes=None, 
                      use_header_as_unique_ID=False, 
                      force_overwrite=False,
                      invalid_sequence_action='fail'):
    """
    Stand alone function that allows the user to build a Proteome from a standard
    FASTA file, or add sequences in a FASTA file to an existing Proteome.

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
        function is provided.

    proteome : Proteome
        If a Proteome object is provided the FASTA file will be read and added to the existing
        proteome, whereas if set to None a new Proteome will be generated.

    build_unique_ID : function
        [**Default = None**] ``build_unique_ID`` allows a user-defined function that is 
        used to convert the FASTA header to a (hopefully) unique string. This can be 
        useful if the FASTA header is well structured and includes a specific, useful
        unique string that can be used as the unique_ID.
        
    build_attributes : function
        [**Default = None**] ``build_attributes`` allows a user-defined function that allows meta-information
        from the FASTA header to be converted into protein attributes. Specifically, build_attributes 
        should be a function which takes in the FASTA header as a string and returns a dictionary where
        key:value pairs are assigned as protein attributes. This can be useful if the FASTA header is well
        structured and includes a specific, useful information relivent to protein of interest. 
    
    use_header_as_unique_ID : bool
        [**Default = False**] ``user_header_as_unique_ID`` is a boolean flag which, if set to true
        means the unique_ID is set to the FASTA file header. NOTE that the combination of this parameter being
        set to true and `build_unique_ID` function not being set to None will trigger an exception as this means
        there are two conflicting definitions of how the unique_ID should be defined. Note that if non-unique
        headers are found this will trigger an exception.

    force_overwrite : bool
        [**Default = False**] Flag that if set to true and we encounter a unique_ID that is already in the proteome
        the newer value overwrites the older one without predudice. This is mostly useful if you are adding in a file
        with known duplicate entries OR combining multiple FASTA files where you know there's some duplications. Note
        that if build_unique_ID = None and user_header_as_unique_ID = None then fasta_to_proteome guarentees that every
        FASTA entry will be given a unique_ID (meaning force_overwrite is irrelevant in this case).

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

    # parameter sanity checking
    if use_header_as_unique_ID is True and build_unique_ID is not None:
        raise APIException('Cannot simultaneously set use_header_as_unique_ID = True and build_unique_ID to not None')
        
    # read in the fasta file using protfasta
    fasta_dictionary = protfasta.read_fasta(filename, invalid_sequence_action=invalid_sequence_action)

    # extract the keys (FASTA headers) and initialize the record_index (internal
    # numbering used for construction. Also initialize the proteom_dict, which is
    # a dictionary of protein entries we passed to Proteome.
    record_index  = 0

    # IF we're adding to a new proteome this bit of code sets the record_index to the largest new integer
    # such that we can add multiple proteomes in succession and we'll get a proteome where there are numerically
    # contigous unique_IDs.  Note we only do this if we'll be using the record_index
    if proteome is not None and (build_unique_ID is None or use_header_as_unique_ID is None):
        numeric_record_ids = []
        for uid in proteome.proteins:
            try:
                numeric_record_ids.append(int(uid))
            except ValueError:
                pass
        if len(numeric_record_ids) > 0:
            record_index = max(numeric_record_ids)+1

    # initialize the empty list
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
        
        
        # get unique_ID 
        if build_unique_ID:
            unique_ID = build_unique_ID(k)
        elif use_header_as_unique_ID is True:
            unique_ID = k
        else:
            unique_ID = record_index
        
        # build an attributes dictionary using the user-provided custom function
        if build_attributes:
            attributes = build_attributes(k)
        else:
            attributes = {}
            
        # now create an input dictionary orbject
        newdict = {}
        newdict['sequence'] = str(fasta_dictionary[k])
        newdict['name'] = k
        newdict['unique_ID'] = unique_ID
        newdict['attributes'] = attributes

        proteome_list.append(newdict)

        record_index = record_index + 1
        
    # finally if a proteome was provided then 
    if proteome is not None:
        proteome.add_proteins(proteome_list, force_overwrite=force_overwrite)
        return proteome
    else:
        # no proteome provided so build a new proteome and return it    
        return Proteome(proteome_list, force_overwrite=force_overwrite)



## ------------------------------------------------------------------------
##
def shephard_fasta_to_proteome(filename, 
                              proteome = None,
                              force_overwrite=False,
                              invalid_sequence_action='fail'):
                              
    """
    Stand alone function that allows the user to build a proteome from a FASTA
    file generated by SHEPHARD (using the proteome_to_fasta() function. When
    SHEPHARD generates a FASTA file it uses a general convention for encoding the
    unique ID, protein name,

    This function assumes the uniprot-standard format for the header
    file has been maintained - i.e.

    >xx|ACCESSION|xxxx

    Where ACCESSION is the uniprot accession and will be used as the unique_ID
    
    Parameters
    ------------

    filename : string
        Name of the FASTA file we're going to parse in. Note the protein name will be
        defined as the full FASTA header for each entry.

    proteome : Proteome
        If a Proteome object is provided the FASTA file will be read and added to the existing
        proteome, whereas if set to None a new Proteome will be generated.


    force_overwrite : bool
        [**Default = False**] Flag that if set to true and we encounter a unique_ID that is already in the proteome
        the newer value overwrites the older one without predudice. This is mostly useful if you are adding in a file
        with known duplicate entries OR combining multiple FASTA files where you know there's some duplications. Note
        that if build_unique_ID = None and user_header_as_unique_ID = None then fasta_to_proteome guarentees that every
        FASTA entry will be given a unique_ID (meaning force_overwrite is irrelevant in this case).

    invalid_sequence_action : ``'ignore'``, ``'fail'``, ``'remove'``, ``'convert'``, ``'convert-ignore'``
        [**Default = 'fail'**] Selector that determines how to deal with invalid sequences that contain invalid/non-standard 
        amino acids. If ``convert`` or ``convert-ignore`` are chosen, then conversion is completed with either the standard 
        conversion table (shown under the ``correction_dictionary`` documentation) or with a custom conversion dictionary 
        passed to ``correction_dictionary``. 

        Options are as follows: 
            * ``ignore``  - invalid sequences are completely ignored
            * ``fail``    - invalid sequence cause parsing to fail and throw an exception
            * ``remove`` - invalid sequences are removed
            * ``convert`` - invalid sequences are convert
            * ``convert-ignore`` - invalid sequences are converted to valid sequences and any remaining invalid residues are ignored
    
    Returns 
    --------
    Proteome Object
        Returns an initialized Proteome object 
    
    """
        
    # read in the fasta file using protfasta
    fasta_dictionary = protfasta.read_fasta(filename, invalid_sequence_action=invalid_sequence_action)

    # initialize the empty list
    proteome_list = []

    # for each entry    
    for k in fasta_dictionary:
        
        # because we know what the header format will be we can be definitive about extracting the relevant information
        fasta_split = k.split('|')

        # ENSURE EVERY single line is a valid 
        if fasta_split[0] != "SHPRD":
            raise APIException('Trying to parse a FASTA file that is expected to be SHEPHARD generated but formatting does not comply [on entry %s in file %s]' % (k, filename))
        
        # extract out 
        try:
            # get the unique ID
            unique_ID = fasta_split[1]

            # then take everything after the unique_ID
            tmp = "|".join(fasta_split[2:])
            attributes_string = tmp.split(SHEPHARD_ATTRIBUTE_SPLITTER)
            name = attributes_string[0]
        except IndexError:
            raise APIException('Trying to parse a FASTA file that is expected to be SHEPHARD generated but formatting does not comply [on entry %s in file %s]' % (k, filename))

        attributes_dict = {}

        if len(attributes_string) > 1:            
            attributes_string_s = attributes_string[1].split('\t')

            for a in attributes_string_s:
                local_k = a.strip().split('=')[0].strip()
                local_v = a.strip().split('=')[0].strip()
                attributes_dict[local_k] = local_v
                
                        
        # now create an protein dictionary object and populate!
        newdict = {}
        newdict['sequence'] = str(fasta_dictionary[k])
        newdict['name'] = name
        newdict['unique_ID'] = unique_ID
        newdict['attributes'] = attributes_dict

        proteome_list.append(newdict)
        
    # finally if a proteome was provided then 
    if proteome is not None:
        proteome.add_proteins(proteome_list, force_overwrite=force_overwrite)
        return proteome
    else:
        # no proteome provided so build a new proteome and return it    
        return Proteome(proteome_list, force_overwrite=force_overwrite)


## ------------------------------------------------------------------------
##
def proteome_to_fasta(filename, proteome, include_attributes_in_header=False):
    """
    Stand alone function that allows the user to write a FASTA file from a
    Proteome object.


    Parameters
    ------------
    filename : string
        Name of the FASTA file we're going to write to. We will automatically overwrite 
        a file if it's there, so be careful! Note that no extension is added in part because
        FASTA files can be .f/.fa/.fasta. Recommended a .fasta file extension.

    proteome : Proteome object
        The proteome object that will be written to discu
        
    include_attributes_in_header : bool (default False)
        Sometimes it may be desirable to annotate all attributes in header as key=value pair.
        This is useful because it means we can read-back protein attributes using the 
        ``build_attributes`` parameter in the ``fasta_to_proteome`` function 

    Returns
    -----------
    None
        No return object but a new file will be written

    """

    # build output list with or without the unique_ID 
    outlist = []
    for protein in proteome:

        # this is where we define the FASTA header...
        fasta_header = "SHPRD|%s|%s" % (protein.unique_ID, protein.name)
            
        # this is where we append the FASTA header with attributes    
        if include_attributes_in_header:

            fasta_header = fasta_header + " " + SHEPHARD_ATTRIBUTE_SPLITTER
            
            for k in protein.attributes:
                
                # these lines ensure there are no tabs in the attribute names of
                # values before we append them to the header file, ensuring that
                # IF we want to read the fasta header info back in to attributes
                # we can be confident that hidden tabs in the variables won't
                # mess things up!
                k_fixed = k.replace('\t', ' ' )
                i = protein.attribute(k)
                i_fixed = i.replace('\t', ' ')
                
                fasta_header = fasta_header + '\t' + "%s=%s" %(k_fixed, i_fixed)
            
        outlist.append([fasta_header, protein.sequence])
        
    # use the protfasta library to write the file to disk
    protfasta.write_fasta(outlist, filename, linelength=80)
                


