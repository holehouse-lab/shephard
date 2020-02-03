"""
Functions associated with interfacing with uniprot data

"""


## ------------------------------------------------------------------------
##
def extract_unique_id_uniprot(keystring):
    """
    Function that converts a header from a uniprot fasta file
    to extract the uniprot ID. This an example of the type of function
    that can be passed to quickstart using the extract_unique_id argument.

    Parameters
    -----------

    keystring : string
        String where we expect the uniprot ID to be contained within two 'pipe' 
        characters ('|'). This function extracts out the uniprot ID a

    Returns
    -----------
    string
        Returns the uniprot ID, although this is not formally validated. However,
        assuming the string follows standard uniprot fasta header conventions this
        should be true!
    

    """
    try:
        return keystring.split('|')[1]
    except:
        raise UtilitiesException('Unable to parse string [%s] to identify uniprot ID' %(keystring))


