"""
uniprot.py

From the SHEPHARD package
Sequence-based Hierachical and Extendable Platform for High-throughput Analysis of Region of Disorder
Ginell & Holehouse, 2020-2022

Handles all I/O associated with uniprot-derived files.

"""

from shephard.exceptions import UtilitiesException
import protfasta
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

    >>> >xx|ACCESSION|xxxx

    where ACCESSION is the uniprot accession. 

    Parameters
    -----------

    line : string
        String where we expect the uniprot ID to be contained within two 'pipe' 
        characters ('|'). 

    Returns
    -----------
    string
        Returns the uniprot ID, although this is not formally validated. 
        However, assuming the string follows standard uniprot fasta header 
        conventions this should be true!

    """
    try:
        return line.split('|')[1].strip()
    except:
        raise UtilitiesException('Unable to parse string [%s] to identify uniprot ID' %(line))

        
        
## ------------------------------------------------------------------------
##
def uniprot_fasta_to_proteome(filename, 
                              proteome = None,
                              force_overwrite=False,
                              invalid_sequence_action='fail'):
                              
    """
    Stand alone function that allows the user to build a proteome from a 
    standard FASTA file downloaded from UniProt

    This function assumes the uniprot-standard format for the header
    file has been maintained - i.e.

    >>> >xx|ACCESSION|xxxx

    Where ACCESSION is the uniprot accession and will be used as the 
    unique_ID
    
    Parameters
    ------------

    filename : string
        Name of the FASTA file we're going to parse in. Note the protein 
        name will be defined as the full FASTA header for each entry.
        
    proteome : Proteome
        If a Proteome object is provided the FASTA file will be read and 
        added to the existing proteome, whereas if set to None a new 
        Proteome will be generated.

    force_overwrite : bool (default  = False)
        If this flag is set to true  and we encounter a unique_ID that is 
        already in the proteome the newer value overwrites the older one. 
        This is mostly useful if you are adding in a file with known 
        duplicate entries OR combining multiple FASTA files where you know 
        there's some duplications. Important - if we're building unique IDs
        based on numerical record indices then EVERY FASTA entry will be given 
        a unique_ID (meaning force_overwrite is irrelevant in this case).

    invalid_sequence_action : str (default = 'fail')
        Selector which defines the behaviour if a sequence with a non-
        standard amino acid is encountered. Valid options and their meaning
        are listed below:

            * ``ignore``  - invalid sequences are completely ignored

            * ``fail``    - invalid sequence cause parsing to fail and throw an exception
  
            * ``remove`` -  invalid sequences are removed

            * ``convert`` - invalid residues are converted to valid residues                            

            * ``convert-ignore`` - invalid sequences are converted to valid sequences and any remaining invalid residues are ignored.
    
    Returns 
    --------
    Proteome
        Returns an initialized Proteome object 
    
    """
    
    return fasta.fasta_to_proteome(filename, proteome=proteome, build_unique_ID=uniprot_accession_from_line, force_overwrite=force_overwrite, invalid_sequence_action=invalid_sequence_action)


## ------------------------------------------------------------------------
##
def uniprot_proteome_to_fasta(filename, proteome):                              
    """
    Stand alone function that allows the user to write a FASTA file from
    a Proteome under the assumption that the Proteome was built from a 
    uniprot FASTA.

    Practically, this just means that the Protein.name variable is used
    for the FASTA header, although the function will fail if duplicate
    headers are found.

    
    Parameters
    ------------

    filename : string
        Name of the FASTA file we're going to write sequences to

    proteome : Proteome
        The Proteome object from which FASTA file will be generated

    
    Returns 
    --------
    None
        No return variable but wll write to file 
    """

    out_dict = {}
    for p in proteome:

        header = p.name

        if header in out_dict:
            raise UtilitiesException(f'Duplicate name entries found in Proteome ({header}). Should not happen for UniProt headers')
        
        out_dict[header] = p.sequence

    protfasta.write_fasta(out_dict, filename, linelength=80)
