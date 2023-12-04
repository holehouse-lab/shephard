"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (alex.holehouse@wustl.edu, g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

from shephard.exceptions import InterfaceException
from . import interface_tools 
from shephard.exceptions import ProteinException, ProteomeException
import shephard.exceptions as shephard_exceptions

MAX_BAD_COUNT  = 10

class _ProteinsInterface:

    """
    Class whose sole purpose is to encapsulate and then store
    parsed Proteins files. This is a hidden class and is not 
    accessible outside of this file.
    
    """

    def __init__(self, filename, delimiter='\t', skip_bad=True):
        r"""
        Expect files of the following format:

        Unique_ID, name, sequence, key1:value1, key2:value2, ..., keyn:valuen

        NOTE that each unique_ID can ONLY appear once!

        Parameters
        ----------------
                
        filename : str
            Name of the shephard proteins file to read

        delimiter : str (default = '\\t')
            String used as a delimiter on the input file. 

        skip_bad : bool (default = True)
            Flag that means if bad lines (lines that trigger an exception) 
            are encountered the code will just skip them. By default this is 
            true, which adds a certain robustness to file parsing, but could 
            also hide errors. Note that if lines are skipped a warning will be 
            printed (regardless of verbose flag). 


        """

        bad_count = 0

        if delimiter == ':':
            raise InterfaceException('When parsing protein file cannot use ":" as a delimeter because this is used to delimit key/value pairs (if provided)')

        with open(filename,'r') as fh:
            content = fh.readlines()
            

        ID2protein = {}
        linecount = 0

        for line in content:

            linecount = linecount + 1

            # skip comment lines
            if interface_tools.is_comment_line(line):
                continue

            sline = line.strip().split(delimiter)

            # try
            try:
                unique_ID = sline[0].strip()
                
                name = sline[1].strip()
                sequence = sline[2].strip()
                attributes = {}                
            except Exception as e:
                msg = f'Failed parsing file [{filename}] on line [{linecount}].\n\nException raised: {str(e)}\n\nline printed below:\n{line}'

                # should update this to also display the actual error...
                if skip_bad and bad_count < MAX_BAD_COUNT:
                    bad_count = bad_count + 1
                    shephard_exceptions.print_warning(msg + f"\nSkipping this line (count {bad_count} of {MAX_BAD_COUNT} ...)")                    
                    continue
                else:
                    raise InterfaceException(msg)

            # if some key/value pairs were included then parse these out one at a time
            if len(sline) > 3:
                attributes = interface_tools.parse_key_value_pairs(sline[3:], filename, linecount, line)
            else:
                # skip over empty entries
                pass
  
            if unique_ID in ID2protein:
                raise InterfaceException("Duplicate protein found in the file %s (offending UID=%s). This cannot be skipped" % (filename, UID))            
            else:
                ID2protein[unique_ID] = {'name':name, 'sequence':sequence, 'attributes':attributes}


        self.data = ID2protein



##############################################
##                                          ##
##     PUBLIC FACING FUNCTIONS BELOW        ##
##                                          ##
##############################################


## ------------------------------------------------------------------------
##
def add_proteins_from_file(proteome, filename, delimiter='\t', return_dictionary = False, safe=True, skip_bad=True, verbose=True):
    r"""
    Function that takes a correctly formatted 'protein' file and reads 
    every protein into the passed proteome.

    The function expects protein files to have the following format:

    >>> Unique_ID name sequence key_1:value_1 key_2:value_2 ... key_n:value_n

    One protein defined per line (with NO duplicates allowed - duplicate 
    entries on the file will trigger an un-rescuable error) where key:values 
    are optional and can be between 0 and n.
        
    **A couple of key points here**:

    * The default delimiter is tabs ('\\t') but this can be changed with the delimiter argument
    * Key value must be separated by a ':', as a result any delimiter (other than ':') can be used, but ':' is reserved for this role.
    * If a protein with the UID from the file exists in the passed proteome then this will throw an exception unless safe=False 
          
    Parameters
    ----------
    proteome : Proteome
        Proteome object to which attributes will be added

    filename : str
        Name of the shephard protein attributes file to read

    Other Parameters
    ----------------

    delimiter : str (default = '\\t')
        String used as a delimiter on the input file. 

    return_dictionary : bool (default = False)
        If set to true, this function will return the protein dictionary 
        and will NOT add that dictionary to the proteome - i.e. the function 
        basically becomes a parser for SHEPHARD-compliant protein files. 
        Default = False

    safe : bool (default = True)
        If set to True then any exceptions raised during the protein-adding 
        process are acted on. Specifically this becomes relevant if we wish 
        to overwrite duplicates (or throw an exception on duplicates).

    skip_bad : bool (default = True)
        Flag that means if bad lines (lines that trigger an exception) are 
        encountered the code will just skip them. By default this is true, 
        which adds a certain robustness to file parsing, but could also hide 
        errors. Note that if lines are skipped a warning will be printed 
        (regardless of verbose flag). skip_bad exclusively influences the 
        file-reading part of the process.
    
    verbose : bool (default = True)
        Flag that defines how 'loud' output is. Will warn about errors on 
        adding attributes.


    Returns
    -----------
    None or dict
        If return_dictionary is set to False (default) then this function 
        has no return value, but the proteins are added to the Proteome 
        object passed as the first argument. If return_dictionary is set 
        to True the function returns the parsed proteins dictionary without
        adding the newly-read proteins to the proteome.
        
    """        
    # check first argument is a proteome
    interface_tools.check_proteome(proteome, 'add_proteins_from_file (si_protein)')

    # next read in the file
    proteins_interface = _ProteinsInterface(filename, 
                                            delimiter=delimiter,
                                            skip_bad=skip_bad)

    if return_dictionary:
        return proteins_interface.data

    # finally add the proteins from the dictionary generated by the ProteinsInterface parser
    add_proteins_from_dictionary(proteome, 
                                 proteins_interface.data, 
                                 safe=safe, 
                                 verbose=verbose)



## ------------------------------------------------------------------------
##
def add_proteins_from_dictionary(proteome, protein_dictionary, safe=True, verbose=True):
    """
    Function that takes a correctly formatted protein dictionary and will 
    add those proteins to the Proteome.

    protein dictionaries are key-value pairs, where the key is a unique 
    ID and the value is itself a dictionary which has the following keys:
    
    * **name** - Protein name (uncontrolled vocabulary, but should be a string)
    * **sequence** - Amino acid sequence for the protein (note that no sanity checking is done here)
    * **attributes** - Dictionary of arbitrary key:value pairings (optional)

    Parameters
    ----------
    proteome : Proteome
        Proteome object to which attributes will be added

    protein_dictionary : dict
        Dictionary that defines proteins. The keys for this dictionary is 
        a unique protein IDs and the values is a list of dictionaries. Each 
        of THOSE sub dictionaries contains key-value pairs are described 
        above.

    safe : bool (default = True)
        If set to True then any exceptions raised during the protein-adding 
        process are acted on. If set to False, exceptions simply mean the 
        protein_attribute in question is skipped. Note if set to False, 
        pre-existing protein_attributes with the same name would be silently 
        overwritten (although this is not consider an error), while overwriting 
        will trigger an exception. 
               
        The only reason protein attribute addition could fail is if the 
        attribute already exists, so this is effectively a flag to define 
        if pre-existing attributes should be overwritten (False) or not (True).
        Default = True.
    
    verbose : bool (default = True)
        Flag that defines how 'loud' output is. Will warn about errors on 
        adding attributes.


    Returns
    -----------
    None
        No return value, but attributes are added to proteins in the Proteome 
        object passed as the first argument.
            
    """

    # check first argument is a Proteome
    interface_tools.check_proteome(proteome, 'add_protein_from_dictionary (si_proteins)')

    if safe is False:
        force_overwrite = True
    else:
        force_overwrite = False

    # for each entry in the overall dictionary
    for UID in protein_dictionary:
            
        # if attributes are included read these out. Note we expect
        # ats to be a dictionary
        try:
            ats = protein_dictionary[UID]['attributes']
        except:
            ats = None
                
        s = protein_dictionary[UID]['sequence']

        # note we use the clean_string to remove tab characters from 
        # the name should they exist
        n = interface_tools.clean_string(protein_dictionary[UID]['name'])
            
        try:
            proteome.add_protein(s, n, UID, attributes=ats, force_overwrite=force_overwrite)
        except (ProteinException, ProteomeException) as e:
            msg='- skipping protein %s (name = %s, len=%i' %(UID, n, len(s))
            if safe:
                shephard_exceptions.print_and_raise_error(msg, e)
            else:
                if verbose:
                    shephard_exceptions.print_warning(msg)
                    continue




                
## ------------------------------------------------------------------------
##
def write_proteins(proteome, filename, delimiter='\t'):
    """
    Function that writes out proteins to file in a standardized format.
    Note that attributes are converted to a string, which for simple 
    attributes is reasonable but is not really a viable stratergy for 
    complex objects, although this will not yeild and error.

    Writes out files with the format:

    >>> Unique_ID name sequence key_1:value_1 key_2:value_2 ... key_n:value_n

    
    Parameters
    -----------

    proteome : Proteome 
        Proteome object from which the proteins will be extracted from

    filename : str
        Filename that will be used to write the new proteins file


    Other Parameters
    ----------------

    delimiter : str (default = '\\t')
        Character (or characters) used to separate between fields. 
        Default is '\\t', which is recommended to maintain compliance 
        with default `add_protein_attributes_from_file()` function


    Returns
    --------
    None
        No return type, but generates a new file with the complete set 
        of protein attributes from this proteome written to disk.        

    """

    with open(filename, 'w') as fh:
        for protein in proteome:

            uid = protein.unique_ID
            n = interface_tools.clean_string(protein.name, delimiter=delimiter)
            s = protein.sequence            

            line = uid
            line = line + delimiter + n
            line = line + delimiter + s
            
            if len(protein.attributes) > 0:

                for k in protein.attributes:

                    atrbt = interface_tools.full_clean_string(protein.attribute(k))

                    line = line + delimiter +  "%s:%s" %(k, atrbt)

                    

            line = line + "\n"

            fh.write(line)
