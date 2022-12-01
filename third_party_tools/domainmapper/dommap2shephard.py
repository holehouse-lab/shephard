#!/usr/bin/env python
"""
dommap2shephard - command line tool for converting output files from DomainMapper into SHEPHARD-compliant domain files.

Note this module can also be used and integrated into code - the main function of interes

"""

import argparse
import time
from os.path import exists


_VERSION='1.0'
_DATE='November 2022'


class Dommap2shephardException(Exception):
    # KEEPING THIS FOR NOW but it's not actually ysed
    pass


def welcome():
    """
    Argument-less function that just prints a welcome message, because
    everyone likes being welcomed!

    """
    print('')
    print('.....................................................')
    print('Running dommap2shephard')
    print(f'Version: {_VERSION} [{_DATE}]')
    print('Started %s '% time.strftime("%Y-%m-%d %H:%M"))
    print('.....................................................\n')




class DomainRecord:

    # -------------------------------------------------------------
    #
    def __init__(self, uid, shprd_domain_type, domain_type, start, end, non_contiguous, insertional, circularly_permuted,  evalue, ecod_fid, ecod_tgroup, ecod_xgroup, ecod_arch, domain_index, domain_max_index, shprd_domain_count):
        """
        DomainRecord is a class for storing information about protein domains. It provides
        a coherent data structure with a couple of methods that enabling simple formatting
        of data for generating SHEPHARD compliant domain domains from each DomainRecord
        object.

        Note: DomainRecords correspond to a single contigous set of residues. For domains
        that are 'simple' domains, this means each domain has a single record. However,
        for domains with non-contigous regions there will be multiple records for a 
        single protein.

        This is kept track of using three distinct indexing approaches.

        Firstly, proteins are indexed based on their unique ID (i.e. UniProt ID). This
        defines the protein from which a DomainRecord emerges. Secondly, all domains 
        (continous or non-continous) posses an index for the number of records they
        are associated with. For continous domains this index is by definition '1 of 1'
        but for non-continous domains it could be '1 of y', where 'y' is the total 
        number of DomainRecords in the domain. FINALLY, in addition to these two,
        each new domain has a unique integer code that increments from 1 and is
        just a unique identifier within a given protien. We do this because it 
        is possible 2 non-contiguous domains overlap perfectly for a subset of 
        their DomainRecords, so this index lets us distiguish between this being
        a desirved behavior vs. a bug.

        To be completed with a second monitor

        Parameters
        ---------------
        uid : str
            Unique protein ID

        shprd_domain_type : str
            Must be of format DomainMapper_{domain_count} where domain_count is
            the index of domain found in this protein. Note that here a domain 
            with 3 non-continous domais would all have the same domain_count.
            
        domain_type : str

        start : str (or int)

        end : str (or int)
        
        non_contiguous :  bool
        
        insertional : bool

        circularly_permuted : bool
        
        evalue : str (or float)
        
        ecod_fid : str

        ecod_tgroup : str 

        ecod_xgroup : str

        ecod_arch : str

        domain_index : str (or int)

        domain_max_index : str (or int) 

        shprd_domain_count : str (or 

        

        """

        # note - we don't cast any of these 
        
        
        self.uid = uid
        self.domain_type = domain_type
        self.shprd_domain_type = shprd_domain_type
        self.start = str(start)
        self.end = str(end)
        self.non_contiguous = non_contiguous
        self.insertional = insertional
        self.circularly_permuted = circularly_permuted
        self.evalue = str(evalue)
        self.ecod_fid = ecod_fid
        self.ecod_tgroup = ecod_tgroup
        self.ecod_xgroup = ecod_xgroup
        self.ecod_arch = ecod_arch
        self.domain_index = domain_index
        self.domain_max_index = domain_max_index
        self.shprd_domain_count = shprd_domain_count


        
        # validation to catch inappropriate assignments. Note this is FAR from fullproof but is better
        # than nothing.
        try:
            tmp = int(self.start)
        except ValueError:
            print('On {self.unique_name}')
            print('Start value could not be parsed - must be a value, and was passed as {self.start}')
            exit(1)
        
        try:
            tmp = int(self.end)
        except ValueError:
            print('On {self.unique_name}')
            print('End value could not be parsed - must be a value, and was passed as {self.end}')
            exit(1)

        try:
            tmp = int(self.domain_index)
        except ValueError:
            print('On {self.unique_name}')
            print('Domain index value could not be parsed - must be a value, and was passed as {self.domain_index}')
            exit(1)

        try:
            tmp = int(self.domain_max_index)
        except ValueError:
            print('On {self.unique_name}')
            print('Domain max index value could not be parsed - must be a value, and was passed as {self.domain_max_index}')
            exit(1)
                        
        try:
            tmp = float(self.evalue)
        except ValueError:
            print(self.evalue)
            print('E-value could not be parsed - must be a value, and was passed as {self.evalue}')
            exit(1)

        if self.insertional not in [True,False]:
            print("Insertional argument not passed as 'True' or 'False'")
            exit(1)
            
        if self.non_contiguous not in [True,False]:
            print("Insertional argument not passed as 'True' or 'False'")
            exit(1)
        
        if self.circularly_permuted not in [True,False]:
            print("Insertional argument not passed as 'True' or 'False'")
            exit(1)

    # -------------------------------------------------------------
    #
    @property
    def unique_name(self):
        """
        Returns a unique string that defines the domain with the following format
        
        <uid>_domain_<shprd domain count>_<domain type>_<start>_<end>_<index>_of_<total components for this domain>

        """
        return f"{self.uid}_domain_{self.shprd_domain_count}_{self.domain_type}_idx_{self.start}_{self.end}_{self.index_string}"

    
    # -------------------------------------------------------------
    #
    @property
    def index_string(self):
        """
        Returns a string that defines the domains relative position in the set of domains
        protein positions associated with this domain. For most 'normal' domains this
        simply returns 1_of_1, but for non-contigous domains this will return x_of_y, where
        x is the index of this particular record and y is the total number of records
        that make up this protein.

        Returns
        -------------
        str
            String as described above

        """
        return f"{self.domain_index}_of_{self.domain_max_index}"

    
    # -------------------------------------------------------------
    #
    def generate_shephard_domain_line(self):

        
        
        return f"{self.uid}\t{self.start}\t{self.end}\t{self.shprd_domain_type}\tdomain_type:{self.domain_type}\tnon_contiguous:{self.non_contiguous}\tinsertional:{self.insertional}\tcircularly_permuted:{self.circularly_permuted}\te_value:{self.evalue}\tECOD_F_ID:{self.ecod_fid}\tECOD_T_group:{self.ecod_tgroup}\tECOD_X_group:{self.ecod_xgroup}\tECOD_architecture:{self.ecod_arch}\tdomain_index:{self.index_string}\tdomain_count:{self.shprd_domain_count}\tdomain_unique_string:{self.unique_name}\tfrom_domain_mapper:True\n"


    ## END OF CLASS
    ## 
    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 


# -------------------------------------------------------------
#
def build_domain_record(uid, sline, domain_count):
    """
    Function that takes a unique ID (UID) and a splitline (sline)
    list and returns a list of DomainRecord objects which can then
    be used to write a SHEPHARD domains file 

    Parameters
    -------------------
    uid : str
        String defining the protein ID

    sline : list
        The split line for the DomainMapper record that will be used
        to construct the one (or more) DomainRecord objects

    domain_count = 

    Returns
    ------------
    list
        The function returns a list of one or more DomainRecord objects,
        where DomainRecords 

    """

    shprd_domain_type = f'DomainMapper_{domain_count}'

    # determine domain types based on presence of
    # flags in the domain_type column
    if sline[3].find('NC') > -1:
        non_contiguous = True
    else:
        non_contiguous = False


    if sline[3].find('IS') > -1:
        insertional = True
    else:
        insertional = False
        
    if sline[3].find('CP') > -1:
        circularly_permuted = True
    else:
        circularly_permuted = False

    # extract out some other info
    evalue      = sline[1]
    ecod_fid    = sline[8] # ECOD F-ID
    ecod_tgroup = sline[6] # ECOD T-group
    ecod_xgroup = sline[5] # ECOD X-group
    ecod_arch   = sline[4] # ECOD Architecture
    domain_type = sline[7] # ECOD F group
    
        
    # having assigned these things we finally deal
    # with group start and end which will determine
    # the number of records associated with a given
    # domain

    start_end = sline[2].split(',')

    # if this group has non NC regions
    if len(start_end) == 1:
        tmp = start_end[0].split('-')
        start = tmp[0]
        end = tmp[1]
        return [DomainRecord(uid, shprd_domain_type, domain_type, start, end, non_contiguous, insertional, circularly_permuted, evalue, ecod_fid, ecod_tgroup, ecod_xgroup, ecod_arch, 1, 1, domain_count)]

    # else we have some NC regions
    else:

        total = len(start_end)
        idx = 1
        return_list = []
        for start_end in start_end:
            tmp = start_end.split('-')
            start = tmp[0]
            end = tmp[1]
            return_list.append(DomainRecord(uid, shprd_domain_type, domain_type, start, end, non_contiguous, insertional, circularly_permuted, evalue, ecod_fid, ecod_tgroup, ecod_xgroup, ecod_arch, idx, total, domain_count))
            idx = idx + 1

        return return_list
    

def read_domainmapper_file(fn, uniprot=False):
    """
    Function that takes in a DomainMapper output files and parses it out into
    a dictionary where keys are domain IDs and values is a list of elements
    that corresponds to the same 9 elements as defined in the DomainMapper
    output file.

    Note that for 'common' user errors an Exception isn't raised, but the 
    function prints an error message and exist with status code (1). These
    scenarios include:

    1. Input file being missing
    2. Error in parsing number of lines in the DomainMapper file (suggests
       this may not be a valid DomainMapper file). 

    Parameters
    ---------------
    fn : str
        Filename for input file. All sanity checking for file parsing is 
        done inside this function.

    uniprot : bool
        Flag which, if set to true, means the assumption is the protein
        IDs are of the format xx|UNIPROT_ID|YYY..., such that the UniProt
        ID is extracted out. Default = False

    Returns
    --------------
    dict
        Returns a dictionary where key is the protein ID and the value
        is a list of DomainRecord objects associated with that protein.

    """

    # check file exists and then read it in
    if not exists(fn):
        print(f'File [{fn}] could not be found. Exiting...')
        exit(1)

    with open(fn,'r') as fh:
        content = fh.readlines()


    # set some variables
    data = {}
    domain_count = {}

    # cycle over each line of the read file
    for idx, line in enumerate(content):

        # skip empty line
        if len(line.strip()) == 0:
            continue
        
        # skip comment line
        if line.strip()[0] == '#':
            continue

        sline = line.strip().split('\t')

        if len(sline) != 9:
            print(f'\nERROR: Could not parse line {idx} of {fn}; implies malformatted DomainMapper file.\n')
            print(f'Full linew written below for clarity:\n{line}\n')
            print('Exiting...')
            exit(1)
            

    
        if uniprot is True:
            try:
                uid = sline[0].split('|')[1]
            except Exception:
                print(f'Error when parsing protein ID on line {idx} using the --uniprot flag. Expecting IDs of the format xx|UniProt Accession|yy')
                print(f'Full linew written below for clarity:\n{line}')
                print('Exiting...')
                exit(1)                

        else:
            uid = sline[0]

        # initialize is first time we're seeing this protein
        if uid not in data:
            data[uid] = []
            domain_count[uid] = 0

        # increment counter 
        domain_count[uid] = domain_count[uid] + 1

        # add record
        data[uid].extend(build_domain_record(uid, sline, domain_count[uid]))

    return data
        

## <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##
## <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


if __name__=="__main__":

    welcome()
    parser = argparse.ArgumentParser()




    parser = argparse.ArgumentParser()
    parser.add_argument("--info", help='Print a help message explaining how to use this dommap2shephard',action='store_true')
    parser.add_argument("-f", help='Inputfile generated by DomainMapper')
    parser.add_argument("--uniprot", help='If provided, assumes the input FASTA file fed into DomainMapper originated from UniProt, so extracts out the UniProt ID, otherwise the full protein ID is used', action='store_true')
    parser.add_argument("-o", help='Defines filename to be written. If not provided, default is "shprd_domain_DomainMapper.tsv"')
    args = parser.parse_args()

    if args.info:

        print('.....................................................')
        print('                      INFO                           ')
        print('.....................................................')
        print('This script converts the output from the dommap executable')
        print('into a SHEPHARD-compliant domains file.\n')
        print('Specifically, if the dommap executable is run as:\n')
        print('   dommap -f input.hmm.out -o mapped_input.hmm.out\n')
        print('Then this output file (mapped_input.hmm.out) can be converted')
        print('to a SHEPHARD-complaint domain file by running:\n')
        print('   python -f dommap2shephard.py mapped_input.hmm.out -o shprd_domain_mapped_domains.tsv\n')
        print('The resulting shephard domains file maps each domain with the following information\n')
        print('ID               : The ID of the protein associated with the domain')
        print('start            : The position where the domain starts (Where first position in a protein = 1)')
        print('end           : The position where the domain ends (inclusive)')
        print('domain_type  : SHEPHARD enables different types of domains, so this reflects the fact this')
        print('               domain came from DomainMapper (as opposed the type of domain it is in terms')
        print('               of structural annotation). Note that All domains will be of type DomainMapper_X')
        print('               where X is an ID that is unique to the domain in that protein. This enables multiple')
        print('               domains of the same type to not clash')
        print('')
        print('Beyond these required components of a domain file, we then also include the following named attributes')
        print('in the SHEPHARD domains file:\n')
        print('domain_type:         : ECOD Family (F-) group - a controled vocabulary of Family (F-group) domains.')
        print('non_contiguous       : Defines if the domain contains additional non-congitous segments. Will be True or False.')            
        print('insertional          : Defines if the domain has been inserted inside another domain. Will be True or False.')           
        print('circularly_permuted  : Defines if the domain is circularly permutaed.Will be True or False.')           
        print('e_value              : Confidence on domain assigment. Smaller is better. Will be a float')
        print('ECOD_F_ID            : ECOD F group ID - e.g. 304.9.1.30')
        print('ECOD_T_group         : ECOD T group name - string defining a protein group based on topology')
        print('ECOD_X_group         : ECOD X group name - string defining a protein group based on similar structure but lack')
        print('                       convincing evidence for bona fide homology')
        print('ECOD_architecture    : ECOD architecture string')
        print('domain_index         : For non_contigous domains, this string is of format X_of_Y where X is the current index and')
        print("                       Y is the total number of non_contigous domains. For 'normal' domains this is just 1_of_1") 
        print('domain_count         : The Y described in the domain index above (i.e. how parts are associated with the underlying')
        print('                     : biological domain that this SHEPHARD domain is associated with. Will be 1 most of the time.')
        print('domain_unique_string : A unique string defining this domain for easy extranl referencing of format:')
        print('                       <UID>_domain_<domain_count>_<domain_type>_idx_<start>_<end>_<domain_index>.')
        print('                       e.g. P04147_domain_1_RRM_1_idx_38_115_1_of_1')
        print('from_domain_mapper   : Flag which is always set to True defining this domain as from DomainMapper.')
        print('\n\n')
        print('If you have any issues with this script please raise an issue on the SHEPHARD GitHub page!')
        print('Last updated 2022-11-30')
        print('.....................................................')
        exit(0)



    # --------------------------------------------------------------------

    if args.f is None:
        print('ERROR: Please provide an DomainMapper input file')
        print('Exiting...')
        exit(1)

    if args.uniprot is True:
        uniprot = True
    else:
        uniprot = False

    if args.o is not None:
        outfile = args.o
    else:
        outfile = 'shprd_domain_DomainMapper.tsv'

    print('Reading input file... ', end='')
    data = read_domainmapper_file(args.f, uniprot)
    print('Done')

    print('Writing output file... ', end='')
    with open(outfile, 'w') as fh:
        
        for protein in data:
            for record in data[protein]:
                fh.write(record.generate_shephard_domain_line())
    print('Done')
    




    
