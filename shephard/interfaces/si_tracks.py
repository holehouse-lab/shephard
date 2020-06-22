"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""


from .interface_exceptions import InterfaceException
from . import interface_tools 
from shephard import general_utilities
from shephard.exceptions import ProteinException

class _TracksInterface:

    def __init__(self, filename, delimiter='\t', mode='values', skip_bad=True):
        """
        
        Class for reading in correctly formatted tracks files for parsing into a
        Proteome object.

        Tracks files must adhere to the following specification

            unique_ID, track_name, val_1, val_2, ...., val_n 

        where n = length of protein.

        This class allows a tracksfile to be read in and defined as either a values
        track file, or a symbols track file, returning a tracks dictionary. 

        


        """
        
        with open(filename,'r') as fh:
            content = fh.readlines()

        ID2track={}

        linecount=0
        for line in content:

            linecount=linecount+1

            sline = line.strip().split(delimiter)
                        
            data_vector=[]
            try:

                # extract track name and unique_id
                unique_ID = sline[0].strip()
                track_name = sline[1].strip()
                
                # parse track values or symbols
                if mode == 'value':

                    # for each element in sline strip whitespace and convert to a float
                    data_vector = [float(i.strip()) for i in sline[2:]]

                else:
                    # for each element in sline strip whitespace 
                    data_vector = [i.strip() for i in sline[2:]]
                        
                if unique_ID not in ID2track:
                    ID2track[unique_ID] = [[track_name, data_vector]]
                else:
                    ID2track[unique_ID].append([track_name, data_vector])

            except Exception:

                msg = 'Failed parsing file [%s] on line [%i]... line printed below:\n%s'%(filename, linecount, line)

                # should update this to also display the actual error...
                if skip_bad:
                    print('Warning: %s'%msg)
                    print('Skipping this line...')
                    continue
                else:
                    raise InterfaceException('Error: %s'%msg)

        self.data = ID2track


##############################################
##                                          ##
##     PUBLIC FACING FUNCTIONS BELOW        ##
##                                          ##
##############################################

## ------------------------------------------------------------------------
##
def add_tracks_from_file(proteome, filename, mode, delimiter='\t', safe=True, skip_bad=True, verbose=True):
    """
    Function that takes a correctly formatted shephard 'tracks' file and reads 
    all Tracks into the passed proteome.

    Expect Track files to have the following format:

    One protein per line, where each protein has the following information:
    
    Unique_ID    track_name    res1    res2    res3 .... resn

    Where res1, res2, resn are symbol or values to be mapped to the 1st, 2nd, or nth residue.
    There should be the same number of res1,2,...n entries are there are residues in the 
    associated protein
    
    A couple of key points here:
    - The default delimiter is tabs ('\t') but this can be changed with the delimiter argument
    - Each track must assign a value or a symbol to EVERY residue in the protein
    
    Parameters
    ----------
    proteome : Proteome Object
        Proteome object 

    filename : str
        Name of the shephard domains file to read

    mode : string {'symbol','value'}
       A selector that defines the type of track file to be read. Must be either 'symbol' or 
       'value'

    delimiter : str 
        String used as a delimiter on the input file. Default = '\t'

    safe : boolean 
        If set to True over-writing tracks will raise an exception. If False, overwriting
        a track will silently over-write. Default = True.

    skip_bad : boolean
        Flag that means if bad lines (lines that trigger an exception) are encountered the code 
        will just skip them. By default this is true, which adds a certain robustness to file 
        parsing, but could also hide errors. Note that if lines are skipped a warning will be 
        printed (regardless of verbose flag). Default = True
    
    verbose : boolean
        Flag that defines how 'loud' output is. Will warn about errors on adding tracks.

    Returns
    -----------
    None
        No return value, but tracks are added to the Proteome object passed as the first argument
        
    """        

    # check first argument is a Proteome
    interface_tools.check_proteome(proteome, 'add_tracks_from_file (si_tracks)')

    # check mode is valid
    general_utilities.valid_keyword('mode', mode, ['symbol','value'])
    
    # next read in the file
    tracks_interface = _TracksInterface(filename, delimiter=delimiter, mode=mode, skip_bad=skip_bad)

    # finally add the domains from the dictionary generated by the DomainsInterface parser
    add_tracks_from_dictionary(proteome, tracks_interface.data, mode, safe=True, verbose=verbose)



## ------------------------------------------------------------------------
##
def add_tracks_from_dictionary(proteome, tracks_dictionary, mode, safe=True, verbose=True):
    """
    Function that takes a 
    Parameters
    ----------

    proteome : Proteome Object
        Proteome object which tracks will be added to

    tracks_dictionary : dict
        Dictionary in which keys are unique IDs for proteins and the value is a list of lists,
        where each sublist where element 0 is the track name and element 1 is itself a list that
        corresponds to the set of positions to be assigned to the track.
    

    mode : string {'symbol','value'}
       A selector that defines the type of track file to be read. Must be either 'symbol' or 
       'value'

    safe : bool (default = True)
        Flag which if true with throw an exception of a track with the same name already exists

    verbose : boolean
        Flag that defines how 'loud' output is. Will warn about errors on adding tracks.
        
    Returns
    -----------
    None
        No return value, but tracks are added to the Proteome object passed as the first argument

        
    """
        
    # check first argument is a proteome
    interface_tools.check_proteome(proteome, 'add_tracks_from_dictionary (si_tracks)')

    # check mode is valid
    general_utilities.valid_keyword('mode', mode, ['symbol','value'])
    
    # cycle through each protein in the proteome...
    for protein in proteome:
        if protein.unique_ID in tracks_dictionary:
            for track in track_dictionary[protein.unique_ID]:

                # get the track name
                track_name = track[0]

                # add the track as either values or symbols depending 
                # on what was provided
                try:
                    if mode == 'values':
                        protein.add_track(track_name, values=track[1], safe=False)
                    else:
                        protein.add_track(track_name, symbols=track[1], safe=False)
                except ProteinException as e:                    
                    msg='- not adding %s track [%s] to [%s] because it already exists ' % (mode, track_name, protein)
                    if safe:
                        # if safe=True fail on any errors
                        print('Error %s' %(msg)) 
                        raise e
                    else:
                        if verbose:
                            print('Warning %s' %(msg))
                        continue



                    
                        
