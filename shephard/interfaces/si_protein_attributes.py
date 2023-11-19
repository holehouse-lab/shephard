"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (alex.holehouse@wustl.edu, g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

from shephard.exceptions import InterfaceException
from . import interface_tools 
from shephard.exceptions import ProteinException
import shephard.exceptions as shephard_exceptions

MAX_BAD_COUNT  = 10

class _ProteinAttributesInterface:

    """
    Class whose sole purpose is to encapsulate and then store
    parsed Protein Attribute files. This is a hidden class and is not 
    accessible outside of this file
    
    """

    def __init__(self, filename, delimiter='\t', skip_bad=True, preauthorized_uids=None):
        r"""
        Expect files of the followin format:

        Unique_ID, key1:value1, key2:value2, ..., keyn:valuen


        Parameters
        ----------------
        
        filename : str
            Name of the shephard domains file to read


        Other Parameters
        ----------------

        delimiter : str (default = '\\t')
            String used as a delimiter on the input file. 

        skip_bad : bool (default = True)
            Flag that means if bad lines (lines that trigger an exception) 
            are encountered the code will just skip them. By default this is 
            true, which adds a certain robustness to file parsing, but could 
            also hide errors. Note that if lines are skipped a warning will be 
            printed (regardless of verbose flag). 
            Default = True

        preauthorized_ids : list of str (default = None)
            List of unique_IDs that are expected to have relevant protein attributes
            If None then all protein attributes are parsed. Avoids parsing
            lines that are not needed into the interface objects.

        """

        bad_count = 0

        if delimiter == ':':
            raise InterfaceException('When parsing domain file cannot use ":" as a delimeter because this is used to delimit key/value pairs (if provided)')

        with open(filename,'r') as fh:
            content = fh.readlines()

        # convert the preauthorized uids to a set for faster lookup
        if preauthorized_uids is not None:
            preauthorized_uids = set(preauthorized_uids)
            
        ID2ADs={}

        linecount=0
        for line in content:

            linecount = linecount + 1

            # skip comment lines
            if interface_tools.is_comment_line(line):
                continue

            sline = line.strip().split(delimiter)

            # try
            try:
                unique_ID = sline[0].strip()
                
                # check if UID associated with this line is found in the
                # preauthorized list. If  not then skip this line
                if preauthorized_uids is not None and unique_ID not in preauthorized_uids:
                    continue
                
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
            if len(sline) > 1:
                attributes = interface_tools.parse_key_value_pairs(sline[1:], filename, linecount, line)
            else:
                # skip over empty entries
                continue
  
            if unique_ID in ID2ADs:
                ID2ADs[unique_ID].append(attributes)
            else:
                ID2ADs[unique_ID] = [attributes]

        self.data = ID2ADs



##############################################
##                                          ##
##     PUBLIC FACING FUNCTIONS BELOW        ##
##                                          ##
##############################################

## ------------------------------------------------------------------------
##
def add_protein_attributes_from_file(proteome, 
                                     filename, 
                                     delimiter='\t', 
                                     return_dictionary=False, 
                                     safe=True, 
                                     skip_bad=True, 
                                     verbose=True):
    r"""
    Function that takes a correctly formatted 'protein attributes' file and 
    reads all attributes into the proteins in the passed proteome.

    The function expects protein attribute files to have the following 
    format:

    One protein defined per line (although the same protein can appear 
    multiple times)

    >>> Unique_ID, key1:value1, key2:value2, ..., keyn:valuen

    A couple of key points here:

    - The default delimiter is tabs ('\\t') but this can be changed with 
      the delimiter argument

    - Key value must be separated by a ':', as a result any delimiter 
      (other than ':') can be used, but ':' is reserved for this role
      
    Parameters
    ----------
    proteome : Proteome Object
        Proteome object to which attributes will be added.

    filename : str
        Name of the shephard protein attributes file to read.

    delimiter : str (default = '\\t')
        String used as a delimiter on the input file. 

    return_dictionary : bool (default = False)
        If set to True, this function will return the protein_attributes 
        dictionary and will NOT add that dictionary to the proteome - 
        i.e. the function basically becomes a parser for SHEPHARD-compliant        
        protein_attributes files. 

    safe : bool (default = True)
        If set to True then any exceptions raised during the 
        protein_attribute-adding process are acted on. If set to False, 
        exceptions simply mean the protein_attribute in question is skipped.         
        Note if set to False, pre-existing protein_attributes with the same 
        name would be silently overwritten (although this is not consider an 
        error), while overwriting will trigger an exception if safe=True.
        
        The only reason protein attribute addition could fail is if the 
        attribute already exists, so this is effectively a flag to define 
        if pre-existing attributes should be overwritten (False) or not 
        (True).

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
        has no return value, but the protein_attributes are added to the 
        Proteome object passed as the first argument. If return_dictionary
        is set to True the function returns the parsed domains_dictionary 
        without adding the newly-read protein_attributes to the proteome.

    """

    # check first argument is a proteome
    interface_tools.check_proteome(proteome, 'add_attributes_from_file (si_protein_attributes)')

    # next read in the file
    protein_attribute_interface = _ProteinAttributesInterface(filename, 
                                                              delimiter=delimiter,
                                                              skip_bad=skip_bad,
                                                              preauthorized_uids=proteome.proteins)

    if return_dictionary:
        return protein_attribute_interface.data


    # finally add the domains from the dictionary generated by the ProteinAttributesInterface parser
    add_protein_attributes_from_dictionary(proteome, 
                                           protein_attribute_interface.data, 
                                           safe=safe, 
                                           verbose=verbose)



## ------------------------------------------------------------------------
##
def add_protein_attributes_from_dictionary(proteome, protein_attribute_dictionary, safe=True, verbose=True):
    r"""
    Function that takes a correctly formatted protein_atttribute dictionary
    and will add those attributes to the proteins in the Proteome.
    
    protein attribute dictionaries are key-value pairs, where the key is a 
    unique ID and the value is a list of dictionaries. For each sub-dictionary, 
    the key-value pair reflects the attribute key-value pairing.

    Parameters
    ----------
    proteome : Proteome Object
        Proteome object to which attributes will be added

    protein_attribute_dictionary : dict
        Dictionary that defines protein attributes. This is slightly 
        confusing, but the keys for this dictionary is a unique 
        protein IDs and the values is a list of dictionaries. Each of 
        THOSE sub-dictionaries has one (or more) key:value pairs that 
        define key:value pairs that will be associated with the protein 
        of interest.

    safe : boolean (default = True)
        If set to True then any exceptions raised during the process of 
        adding a protein_attribute are further raised. If set to False, 
        exceptions simply mean the protein_attribute in question is skipped.         
        Note if set to False, pre-existing protein_attributes with the same 
        name would be silently overwritten (although this is not consider an 
        error), while overwriting will trigger an exception if safe=True.
        Default = True
        
        The only reason protein attribute addition could fail is if the 
        attribute already exists, so this is effectively a flag to define 
        if pre-existing attributes should be overwritten (False) or not 
        (True).
    
    verbose : bool (default = True)
        Flag that defines how 'loud' output is. Will warn about errors on 
        adding attributes.

    Returns
    -----------
    None
        No return value, but attributes are added to proteins in the Proteome 
        object passed as the first argument
    """

    # check first argument is a Proteome
    interface_tools.check_proteome(proteome, 'add_protein_attributes (si_protein_attributes)')
    
    for protein in proteome:
        if protein.unique_ID in protein_attribute_dictionary:

            # note here each AD is its own dictionary
            for AD in protein_attribute_dictionary[protein.unique_ID]:

                # for each attribute-key
                for k in AD:            

                    # get the value
                    v = AD[k]

                    try:
                        protein.add_attribute(k, v, safe=safe)
                    except ProteinException as e:
                        msg='- skipping attribute entry on protein %s (key: %s) ' % (protein.unique_ID, k)
                        if safe:
                            shephard_exceptions.print_and_raise_error(msg, e)
                        else:
                            if verbose:
                                shephard_exceptions.print_warning(msg)
                                continue




## ------------------------------------------------------------------------
##
def write_protein_attributes(proteome, filename, delimiter='\t'):
    r"""
    Function that writes out protein attributes to file in a standardized 
    format. Note that attributes are converted to a string, which for simple 
    attributes is reasonable but is not really a viable stratergy for 
    complex objects, although this will not yeild and error.
    
    Parameters
    -----------
    proteome :  Proteome object
        Proteome object from which the domains will be extracted from

    filename : str
        Filename that will be used to write the new domains file

    delimiter : str (default = '\\t')
        Character (or characters) used to separate between fields. 
        Default is '\t', which is recommended to maintain compliance 
        with default `add_protein_attributes_from_file()` function.        
        
    Returns
    --------
    None
        No return type, but generates a new file with the complete set 
        of protein attributes from this proteome written to disk.
        
    """

    with open(filename, 'w') as fh:
        for protein in proteome:
            if len(protein.attributes) > 0:

                line = protein.unique_ID

                for k in protein.attributes:

                    atrbt = interface_tools.full_clean_string(protein.attribute(k))

                    line = line + delimiter +  "%s:%s" %(k, atrbt)

                line = line + "\n"

                fh.write(line)


## ------------------------------------------------------------------------
##
def write_protein_attributes_from_dictionary(protein_attribute_dictionary, filename, delimiter='\t'):
    r"""
    Function that writes out protein attributes to file in a standardized 
    format. Note that attributes are converted to a string, which for simple 
    attributes is reasonable but is not really a viable stratergy for 
    complex objects, although this will not yeild and error.
    
    Parameters
    -----------
    protein_attribute_dictionary :  dictionary
        protein_attribute_dictionary for which the protein IDs are keys 
        and the values are dictionaries with key:value pairs of attributes
        which are to be writen
    filename : str
        Filename that will be used to write the new domains file
    delimiter : str (default = '\\t')
        Character (or characters) used to separate between fields. 
        Default is '\t', which is recommended to maintain compliance 
        with default `add_protein_attributes_from_file()` function.        
        
    Returns
    --------
    None
        No return type, but generates a new file with the complete set 
        of protein attributes from this proteome written to disk.
    """

    with open(filename, 'w') as fh:

        for protein in protein_attribute_dictionary:

            local_attributes = protein_attribute_dictionary[protein]
            
            if len(local_attributes) > 0:

                line = protein

                for k,v in local_attributes.items():

                    atrbt = interface_tools.full_clean_string(v)

                    line = line + delimiter +  "%s:%s" %(k, atrbt)

                line = line + "\n"

                fh.write(line)
