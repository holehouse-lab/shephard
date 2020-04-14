from .interface_exceptions import InterfaceException
from . import interface_tools 
from shephard.exceptions import ProteinException

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
                attribute_dictionary = {}                
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
                attribute_dictionary = interface_tools.parse_key_value_pairs(sline[1:], filename, linecount, line)
            else:
                # skip over empty entries
                continue
  
            if unique_ID in ID2ADs:
                ID2ADs[unique_ID].append(attribute_dictionary)
            else:
                ID2ADs[unique_ID] = [attribute_dictionary]

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
    attributes to the proteins in the proteome.

    protein_attribute dictionaries are key-value pairs, where the key is a unique ID and the value
    is a list of dictionaries, where for each dictionary the key-value pair is a key-value pair of protein
    attributes.

    
    
    """

    # check first argument is a proteome
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
                        protein.add_attribute(k, v, safe=False)
                    except ProteinException as e:
                        msg='- skipping attribute entry at %i-$i on %s' %(start, end, protein)
                        if safe:
                            # if safe=True fail on any errors
                            print('Error %s' %(msg)) 
                            raise e
                        else:
                            if verbose:
                                print('Warning %s' %(msg))
                            continue


## ------------------------------------------------------------------------
##
def add_protein_attributes_from_file(proteome, filename, delimiter='\t', safe=True, skip_bad=True, verbose=True):
    """
    Function that takes a correctly formatted shephard 'domains' file and reads 
    all domains into the passed proteome.

    Expect Domain files to have the followin format:

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
        Proteome object 

    filename : string

    <TO DO>

        
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
    Function that writes out domains to file in a standardized format.
    
    <TO DO>

    """

    with open(filename, 'w') as fh:
        for protein in proteome:
            if len(protein.attributes) > 0:

                fh.write('%s%s'%(protein.unique_ID, delimiter))

                for k in protein.attributes:
                    fh.write('%s:%s%s' % (k, protein.attribute(k), delimiter))
                fh.write('\n')
