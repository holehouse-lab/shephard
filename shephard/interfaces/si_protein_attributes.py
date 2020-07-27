"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

from shephard.exceptions import InterfaceException
from . import interface_tools 
from shephard.exceptions import ProteinException
import shephard.exceptions as shephard_exceptions

class _ProteinAttributesInterface:

    """
    Class whose sole purpose is to encapsulate and then store
    parsed Protein Attribute files. This is a hidden class and is not 
    accessible outside of this file
    
    """

    def __init__(self, filename, delimiter='\t', skip_bad=True):
        """
        Expect files of the followin format:

        Unique_ID, key1:value1, key2:value2, ..., keyn:valuen

        """

        if delimiter == ':':
            raise InterfaceException('When parsing domain file cannot use ":" as a delimeter because this is used to delimit key/value pairs (if provided)')

        with open(filename,'r') as fh:
            content = fh.readlines()

        ID2ADs={}

        linecount=0
        for line in content:

            linecount = linecount + 1

            sline = line.strip().split(delimiter)

            # try
            try:
                unique_ID = sline[0].strip()
                attributes = {}                
            except Exception:

                msg = 'Failed parsing file [%s] on line [%i]... line printed below:\n%s'%(filename, linecount, line)

                # should update this to also display the actual error...
                if skip_bad:
                    print('Warning: %s'%msg)
                    print('Skipping this line...')
                    continue
                else:
                    raise InterfaceException('Error: %s'%msg)

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
def add_protein_attributes_from_dictionary(proteome, protein_attribute_dictionary, safe=True, verbose=True):
    """
    Function that takes a correctly formatted protein_atttribute dictionary and will add those 
    attributes to the proteins in the Proteome.

    protein_attribute dictionaries are key-value pairs, where the key is a unique ID and the value
    is a list of dictionaries, where for each dictionary the key-value pair is a key-value pair of protein
    attributes.

    Parameters
    ----------
    proteome : Proteome Object
        Proteome object to which domains will be added

    protein_attribute_dictionary : dict
        Dictionary that defines protein attributes. This is slightly confusing, but the keys for this
        dictionary is a unique protein IDs and the values is a list of dictionaries. Each of THOSE sub
        dictionaries has one (or more) key:value pairs that define key:value pairs that will be associated
        with the protein of interest.

    safe : boolean 
        If set to True then any exceptions raised during the protein_attribute-adding process are acted
        on. If set to False, exceptions simply mean the protein_attribute in question is skipped. 
        Note if set to False, pre-existing protein_attributes with the same name would be silently 
        overwritten (although this is not consider an error), while overwriting will trigger an 
        exception in safe=True.
        
        The only reason protein attribute addition could fail is if the attribute already exists, so
        this is effectively a flag to define if pre-existing attributes should be overwritten (False) 
        or not (True).

        Default = True.
    
    verbose : boolean
        Flag that defines how 'loud' output is. Will warn about errors on adding domains.

    Returns
    -----------
    None
        No return value, but domains are added to the Proteome object passed as the first argument
    
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
                                shephard_exceptions.print_warning(msg,e)
                                continue


## ------------------------------------------------------------------------
##
def add_protein_attributes_from_file(proteome, filename, delimiter='\t', safe=True, skip_bad=True, verbose=True):
    """
    Function that takes a correctly formatted Proteome 'domains' file and reads 
    all domains into the passed Proteome.

    Expect Domain files to have the following format:

    One domain per line where:
    
    Unique_ID    start    stop    domain_type    key_1:value_1    key_2:value_2, ...,     key_n:value_n

    A couple of key points here:
    - The default delimiter is tabs ('\t') but this can be changed with the delimiter argument
    - The first four arguments are required, while all of the key:value pairs are optional
    - Key value must be separated by a ':', as a result any delimiter (other than ':') can be 
      used, but ':' is reserved for this role
    

    Parameters
    ----------
    proteome : Proteome Object
        Proteome object to which domains will be added

    filename : str
        Name of the shephard domains file to read

    delimiter : str 
        String used as a delimiter on the input file. Default = '\t'

    safe : boolean 
        If set to True then any exceptions raised during the protein_attribute-adding process are acted
        on. If set to False, exceptions simply mean the protein_attribute in question is skipped. 
        Note if set to False, pre-existing protein_attributes with the same name would be silently 
        overwritten (although this is not consider an error), while overwriting will trigger an 
        exception in safe=True.
        
        The only reason protein attribute addition could fail is if the attribute already exists, so
        this is effectively a flag to define if pre-existing attributes should be overwritten (False) 
        or not (True).

        Default = True.

    safe : boolean 
        If set to True over-writing attributes will raise an exception. If False, overwriting
        an attribute will silently over-write. Default = True.

    skip_bad : boolean
        Flag that means if bad lines (lines that trigger an exception) are encountered the code 
        will just skip them. By default this is true, which adds a certain robustness to file 
        parsing, but could also hide errors. Note that if lines are skipped a warning will be 
        printed (regardless of verbose flag). Default = True
    
    verbose : boolean
        Flag that defines how 'loud' output is. Will warn about errors on adding domains.

    Returns
    -----------
    None
        No return value, but domains are added to the Proteome object passed as the first argument
            

        
    """        
    # check first argument is a proteome
    interface_tools.check_proteome(proteome, 'add_domains_from_file (si_domains)')

    # next read in the file
    protein_attribute_interface = _ProteinAttributesInterface(filename, 
                                                              delimiter, 
                                                              skip_bad=skip_bad)

    # finally add the domains from the dictionary generated by the DomainsInterface parser
    add_protein_attributes_from_dictionary(proteome, 
                                           protein_attribute_interface.data, 
                                           safe=safe, 
                                           verbose=verbose)



                

## ------------------------------------------------------------------------
##
def write_protein_attributes(proteome, filename, delimiter='\t'):
    """
    Function that writes out protein attributes to file in a standardized format.

    
    Parameters
    -----------

    proteome :  Proteome object
        Proteome object from which the domains will be extracted from

    filename : str
        Filename that will be used to write the new domains file

    delimiter : str
        Character (or characters) used to separate between fields. Default is '\t'
        Which is recommended to maintain compliance with default `add_protein_attributes_from
        file()` function

    Returns
    --------
    None
        No return type, but generates a new file with the complete set of protein attributes
        from this proteome written to disk.

    """

    with open(filename, 'w') as fh:
        for protein in proteome:
            if len(protein.attributes) > 0:

                line = protein.unique_ID

                for k in protein.attributes:
                    line = line + delimiter
                    line = line + "%s:%s" %(k, protein.attribute(k))

                line = line + "\n"

                fh.write(line)
