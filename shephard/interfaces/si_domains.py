from .interface_exceptions import InterfaceException
from . import interface_tools 

class DomainsInterface:

    def __init__(self, filename, delimiter='\t', safe=True):
        """
        Expect files of the followin format:

        Unique_ID, start, stop, domain_type, key1:value1, key2:value2, ..., keyn:valuen

        Note that the first four arguments are required, while all of the key:value pairs 
        are optional. Key value must be separated by a ':', but any delimiter (other than ':') 
        is allowed

        """

        if delimiter == ':':
            raise InterfaceException('When parsing domain file cannot use ":" as a delimeter because this is used to delimit key/value pairs (if provided)')

        with open(filename,'r') as fh:
            content = fh.readlines()

        ID2domain={}

        linecount=0
        for line in content:

            linecount=linecount+1

            sline = line.split(delimiter)
            
            try:
                unique_ID = sline[0].strip()
                start = int(sline[1].strip())
                end = int(sline[2].strip())
                domain_type = sline[3].strip()
                attribute_dictionary = {}
            except Exception:

                # should update this to also display the actual error...
                raise InterfaceException('Failed parsing file [%s] on line [%i]... line printed below:\n%s'%(filename, linecount, line))

            # if some key/value pairs were included then parse these out one at a time
            if len(sline) > 4:
                attribute_dictionary = interface_tools.parse_key_value_pairs(sline[4:], filename, linecount, line)
                              
  
            if unique_ID in ID2domain:
                ID2domain[unique_ID].append([start, end, domain_type, attribute_dictionary])
            else:
                ID2domain[unique_ID] =[[start, end, domain_type, attribute_dictionary]]

        self.data = ID2domain


def read_in_domains(proteome, filename, delimiter='\t', safe=True, autoname=False):
    """
    Function that takes a correctly formatted shephard 'domains' file and reads 
    all domains into the passed proteome
    
    Parameters
    ----------
    proteome : 

        
    """

        
    # check first argument is a proteome
    interface_tools.check_proteome(proteome, 'add_domains (si_domains)')

    domains_interface = DomainsInterface(filename, delimiter)

    for protein in proteome:
        if protein.unique_ID in domains_interface.data:
            for domain in domains_interface.data[protein.unique_ID]:

                start = domain[0]
                end = domain[1]
                domain_type = domain[2]
                ad = domain[3]
                

                try:
                    protein.add_domain(start, end, domain_type, ad, safe, autoname)
                except InterfaceException as e:
                    if skip_bad:
                        print('Warning - skippiing domain at %i-$i on %s' %(start, end, protein))
                        continue
                    else:
                        raise e

                


def write_domains(proteome, filename, delimiter='\t', skip_bad=False):

    with open(filename, 'w') as fh:
        for protein in proteome:
            for d in protein.domains:
                start = protein.domain(d).start
                end = protein.domain(d).end
                domain_type = protein.domain(d).domain_type
                fh.write('%s%s%i%s%i%s%s\n' % (protein.unique_ID, 
                                               delimiter, 
                                               start,
                                               delimiter,
                                               end,
                                               delimiter,
                                               domain_type))
    
