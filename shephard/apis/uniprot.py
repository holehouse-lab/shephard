"""
uniprot.py

From the SHEPHARD package
Sequence-based Hierachical and Extendable Platform for High-throughput Analysis of Region of Disorder
Ginell & Holehouse, 2020

Handles all I/O associated with uniprot-derived files.

"""


from . import fasta

## ------------------------------------------------------------------------
##
def uniprot_accession_from_line(line):
    """
    Function that converts a header from a uniprot fasta file
    to extract the uniprot ID. This an example of the type of function
    that can be passed to quickstart using the extract_unique_id argument.

    This function assumes the uniprot-standard format for the header
    file has been maintained - i.e.

    >xx|ACCESSION|xxxx

    where ACCESSION is the uniprot accession. 

    Parameters
    -----------

    line : string
        String where we expect the uniprot ID to be contained within two 'pipe' 
        characters ('|'). 

    Returns
    -----------
    string
        Returns the uniprot ID, although this is not formally validated. However,
        assuming the string follows standard uniprot fasta header conventions this
        should be true!
    

    """
    try:
        return line.split('|')[1].strip()
    except:
        raise UtilitiesException('Unable to parse string [%s] to identify uniprot ID' %(keystring))


## ------------------------------------------------------------------------
##
def uniprot_attributes_from_header(line):
    """
    Function that converts the header information from a uniprot fasta file
    to key value pairs in an attribute in dictionary. This an example of the type of function
    that can be passed to quickstart using the build_attributes argument.

    This function assumes the uniprot-standard format for the header
    file has been maintained - i.e.

    >xx|ACCESSION|EntryName Isoform/ProtienName key=value ...

    where ACCESSION is the uniprot accession. 

    Parameters
    -----------

    line : string
        String where we expect the EntryName to be contained within the second 'pipe' 
        characters ('|') and a space (' '). 

    Returns
    -----------
    dictionary
        Returns the attribute dictionary, although this is not formally validated. However,
        assuming the string follows standard uniprot fasta header conventions this
        should be true!
    

    """
    try:
        return line.split('|')[1].strip()
    except:
        raise UtilitiesException('Unable to parse string [%s] to identify uniprot ID' %(keystring))


        
        
## ------------------------------------------------------------------------
##
def uniprot_fasta_to_proteome(filename, invalid_sequence_action='fail', build_attributes=False):
    """
    Stand alone function that allows the user to build a proteome from a standard
    FASTA file downloaded from UniProt

    This function assumes the uniprot-standard format for the header
    file has been maintained - i.e.

    >xx|ACCESSION|xxxx

    Where ACCESSION is the uniprot accession and will be used as the unique_ID
    
    Parameters
    ------------

    filename : string
        Name of the FASTA file we're going to parse in. Note the protein name will be
        defined as the full FASTA header for each entry.


    invalid_sequence_action : ``'ignore'``, ``'fail'``, ``'remove'``, ``'convert'``, ``'convert-ignore'``
        [**Default = 'fail'**] Selector that determines how to deal with invalid sequences that contain invalid/non-standard amino acids. If ``convert`` or ``convert-ignore`` are chosen, then conversion is completed with either the standard conversion table (shown under the ``correction_dictionary`` documentation) or with a custom conversion dictionary passed to ``correction_dictionary``. 

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
    
    return fasta.fasta_to_proteome(filename, build_unique_ID=uniprot_accession_from_line, invalid_sequence_action=invalid_sequence_action)
