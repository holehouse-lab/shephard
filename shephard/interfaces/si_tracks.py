"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (alex.holehouse@wustl.edu, g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""


from shephard.exceptions import InterfaceException, ProteinException, TrackException
import shephard.exceptions as shephard_exceptions
from . import interface_tools 
from shephard import general_utilities
import os



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

        ID2track = {}

        linecount = 0

        # cycle over every line in the file
        for line in content:

            linecount = linecount + 1

            # extract chop off lagging whitespace and divide up using the delimiter
            sline = line.strip().split(delimiter)                        
            track_data = []
            
            # for this list 
            try:

                # extract track name and unique_id
                unique_ID = sline[0].strip()
                track_name = sline[1].strip()
                
                # parse track values or symbols
                if mode == 'values':

                    # for each element in sline strip whitespace and convert to a float
                    track_data = [float(i.strip()) for i in sline[2:]]

                elif mode == 'symbols':
                    # for each element in sline strip whitespace 
                    track_data = [i.strip() for i in sline[2:]]
                else:
                    raise InterfaceException('Error: %s'% "mode passed = %s, yet this does not match 'symbols' or 'values'")

                        
                if unique_ID in ID2track:
                    ID2track[unique_ID].append({'track_name':track_name, 'track_data':track_data})                    
                else:
                    ID2track[unique_ID] = [{'track_name':track_name, 'track_data':track_data}]
                    

            except Exception as e:

                msg = 'Failed parsing file [%s] on line [%i].\n\nException raised: %s\n\nline printed below:\n%s'%(filename, linecount, str(e), line)

                # should update this to also display the actual error...
                if skip_bad:
                    shephard_exceptions.print_warning(msg + "\nSkipping this line...")
                    continue
                else:
                    raise InterfaceException(msg)

        self.data = ID2track


##############################################
##                                          ##
##     PUBLIC FACING FUNCTIONS BELOW        ##
##                                          ##
##############################################

## ------------------------------------------------------------------------
##
def add_tracks_from_file(proteome, filename, mode, delimiter='\t', return_dictionary=False, safe=True, skip_bad=True, verbose=True):
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

    mode : string {'symbols','values'}
       A selector that defines the type of track file to be read. Must be either 'symbols' or 
       'values'

    delimiter : str 
        String used as a delimiter on the input file. Default = '\t'

    return_dictionary : bool
        If set to true, this function will return the tracks dictionary and will NOT add that
        dictionary to the proteome - i.e. the function basically becomes a parser for SHEPHARD-compliant
        tracks files. Default = False

    safe : boolean 
        If set to True then any exceptions raised during the Track-adding process (i.e. after file
        parsing) are acted on. If set to False, exceptions simply mean the site in question is skipped. 
        Note if set to False  pre-existing tracks with the same name would be silently overwritten (although 
        this is not consider an error), while overwriting will trigger an exception in safe=True.
        There are various reasons site addition could fail (e.g. track does not match length of protein)
        so if verbose=True then the cause of an exception is also printed to screen. It is highly 
        recommend that if you choose to use safe=False you also set verbose=True. Default = True.

    skip_bad : boolean
        Flag that means if bad lines (lines that trigger an exception) are encountered the code 
        will just skip them. By default this is true, which adds a certain robustness to file 
        parsing, but could also hide errors. Note that if lines are skipped a warning will be 
        printed (regardless of verbose flag). Default = True
    
    verbose : boolean
        Flag that defines how 'loud' output is. Will warn about errors on adding tracks.

    Returns
    -----------
    None or dict
        If return_dictionary is set to False (default) then this function has no return
        value, but the tracks are added to the Proteome object passed as the first argument. If
        return_dictionary is set to True the function returns the parsed tracks dictionary without
        adding the newly-read tracks to the proteome.
        
    """        

    # check first argument is a Proteome
    interface_tools.check_proteome(proteome, 'add_tracks_from_file (si_tracks)')

    # check mode is valid
    general_utilities.valid_keyword('mode', mode, ['symbols','values'])
    
    # next read in the file
    tracks_interface = _TracksInterface(filename, delimiter=delimiter, mode=mode, skip_bad=skip_bad)

    if return_dictionary:
        return tracks_interface.data

    # finally add the domains from the dictionary generated by the DomainsInterface parser
    add_tracks_from_dictionary(proteome, tracks_interface.data, mode, safe=safe, verbose=verbose)



## ------------------------------------------------------------------------
##
def add_tracks_from_dictionary(proteome, tracks_dictionary, mode, safe=True, verbose=True):
    """

    Function that takes a correctly formatted tracks dictionary and will add those tracks to 
    the proteins in the Proteome.
    

    track dictionaries are key-value pairs, where the key is a unique ID and the value
    is a list of dictionaries. For each sub-dictionary, there are two key-value pairs that
    reflect:

        'track_name'  : name of the track (str)
        'track_data' : parsed list of floats (if expecting values) or strings (if expecting symbols)
                        that should equal the length of the associated protein.


    Parameters
    ----------

    proteome : Proteome Object
        Proteome object which tracks will be added to

    tracks_dictionary : dict
        Dictionary in which keys are unique IDs for proteins and the value is a list of dictionaries,
        where each subdictionary has the two key-value pairs:

        'track_name'  : name of the track (str)
        'track_data' : parsed list of floats (if expecting values) or strings (if expecting symbols)
                        that should equal the length of the associated protein.
    
    mode : string {'symbols','values'}
       A selector that defines the type of track file to be read. Must be either 'symbols' or 
       'values'

    safe : bool (default = True)
        If set to True then any exceptions raised during the track-adding process are acted
        on. If set to False, exceptions simply mean the Track in question is skipped. 
        Note if set to False, pre-existing Tracks with the same name would be silently overwritten (although 
        this is not consider an error), while overwriting will trigger an exception in safe=True
        There are various reasons Track addition could fail (length does not match the protein etc) 
        and so if verbose=True then the cause of an exception is also printed to 
        screen. It is highly recommend that if you choose to use safe=False you also set verbose=True. 
        Default = True.

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
    general_utilities.valid_keyword('mode', mode, ['symbols','values'])
    
    # cycle through each protein in the proteome...
    for protein in proteome:
        if protein.unique_ID in tracks_dictionary:
            for track in tracks_dictionary[protein.unique_ID]:

                # get the track name and vector info
                track_name = track['track_name']
                track_data = track['track_data']

                # add the track as either values or symbols depending 
                # on what was provided
                try:
                    if mode == 'values':
                        protein.add_track(track_name, values=track_data, safe=safe)
                    else:
                        protein.add_track(track_name, symbols=track_data, safe=safe)

                # if an ProteinException was raised when trying to add a track some
                # anticipated error occurred
                except (ProteinException, TrackException) as e:      

                    msg='- skipping track at %s on %s' %(track_name, protein)
                    if safe:
                        shephard_exceptions.print_and_raise_error(msg, e)
                    else:
                        if verbose:
                            shephard_exceptions.print_warning(msg)
                        continue




## ------------------------------------------------------------------------
##
def write_all_tracks(proteome, outdirectory='.', value_fmt = "%.3f", delimiter='\t'):
    """
    Function that writes all tracks associated with a proteome out, whereby the output
    filenames are defined as:
    
    shephard_track_<trackname>.tsv
    
    and are written to the outdirectory.

    Because track files MUST be written as one per track_name, this function is equivalent
    to cycling through each unique track name and writing it out sequentially. 
    
    Parameters
    -----------

    proteome :  Proteome object
        Proteome object from which the domains will be extracted from

    outdirectory : str
        String that defines the output directory. By default sets to the present
        working directory ('.').

    value_fmt : str
        Format string that will be used for values. Default = "%.3f"
        
    delimiter : str
        Character (or characters) used to separate between fields. Default is '\t'
        Which is recommended to maintain compliance with default `add_tracks_from_files()`
        function.
    
    Returns
    --------
    None
        No return type, but generates a new file with the complete set of domains
        from this proteome written to disk.


    """

    for t_name in protein.unique_track_names:
        outname = os.path.join(outdirectory, "shephard_track_%s.tsv" %( t_name))
        write_tracks(proteome, outname, t_name, value_fmt, delimiter)
    

            



## ------------------------------------------------------------------------
##
def write_track(proteome, filename, track_name, value_fmt = "%.3f", delimiter='\t'):
    """
    Function that writes out a specific track to file in a standardized format. Note that
    because track files are inevitably quite big default behaviour is to only write out a
    single track file at a time (i.e. unlike write_domains or write_sites where ALL domains
    or all sites are - by default - written out, here ONLY a single type of track, defined
    by track_name, can be written.

    To write ALL the tracks from a file, see si_tracks.write_all_tracks().
    
    Parameters
    -----------
    proteome :  Proteome object
        Proteome object from which the domains will be extracted from

    filename : str
        Filename that will be used to write the new domains file

    track_name : str
        Name of the track to be written out.

    value_fmt : str
        Format string that will be used for values. Default = "%.3f". Note that this is not
        a smart value so if the actual value used means that %.3f looses all meaning this will
        not trigger a warning, so, be careful!
        
    delimiter : str
        Character (or characters) used to separate between fields. Default is '\t'
        Which is recommended to maintain compliance with default `add_tracks_from_files()`
        function.
    
    Returns
    --------
    None
        No return type, but generates a new file with the complete set of domains
        from this proteome written to disk.

    """

    # test the passed value_fmt string works. This is not fullproof but at least validates that
    # the string can parse a float (which is a necessary requirement for tracks values to be read
    # back in again by shephard
    try:
        a = value_fmt %( 1.5 )

        if float(a) != 1.5:
            raise InterfaceException('Invalid value_fmt passed [%s]'%(str(value_fmt)))
    except TypeError:
        raise InterfaceException('Invalid value_fmt passed [%s]'%(str(value_fmt)))

            
    with open(filename, 'w') as fh:
        
        for protein in proteome:

            
            # try and extract out the track in question
            t = protein.track(track_name, safe=False)
            if t is not None:
                unique_ID = protein.unique_ID

                # build the initial string
                out_string = "%s%s%s%s" % (unique_ID, delimiter, t.name, delimiter)

                if t.values is not None:
                    for v in t.values:
                        out_string = out_string + "%s%s" % (value_fmt %(v), delimiter)
                else:
                    for v in t.symbols:
                        out_string = out_string + "%s%s" % (v, delimiter)


                fh.write('%s\n'%(out_string))
                    
                    


    
