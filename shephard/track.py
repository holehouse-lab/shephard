"""
Track class File - ShanEnt Suite 

version: 3

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

from . exceptions import TrackException
import numpy as np

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Class that defines a Site in sequence
#
class Track:
    """
    Tracks define information that maps along a protein sequence.

    A track is, fundementally, a vector which is the length of the sequence. This could be an way to re-code
    the amino acid sequence, or reflect some kind of sliding window analysis.

    Tracks can can either define a set of symbols that convert residues to symbols (i.e. discrete classifications)
    or values (i.e. floating number values associated with each position).


    """
    
    def __init__(self, name, protein, values=None, symbols=None):
        """

        Parameters
        ------------

        name : string
           Defines the name of the track. This can be any value, but should be something that makes sense. The name can be 
           used by analysis routines.

        protein : 


        """

        # if values was provided for the track...
        if values is not None:

            # check the values provided is the same length as the number of residues - if not raise an exception
            if len(protein.sequence) != len(values):
                raise TrackException('Track length of %i does not match protein length of %i (values track)\b Track = %s\nProtein=%s' %(len(values), len(protein.sequence), name, str(protein)))

            try:                
                # values should be a numpy array of float64 type
                values = np.array(values, dtype='float64')
            except ValueError:
                raise TrackException('Unable to convert values passed into float64 numpy array [Track=%s, Protein=%s' %( name, str(protein)))
            

        # if the symbols were provided
        if symbols is not None:

            # check lengths match up
            if len(protein.sequence) != len(symbols):
                raise TrackException('Track length (symbols) does not match protein length [Track=%s, Protein=%s' %( name, str(protein)))

            # if we passed a list (which we now know is the right length) we're good!
            if isinstance(symbols, list):
                pass

            # if we passed a string convert to a list
            elif(symbols, str):
                symbols = list(symbols)
            
            else:
                raise TrackExceptio('Unable to convert passed symbols track to a list of symbols [Track=%s, Protein=%s' %( name, str(protein)))


        # if NEITHER symbols nor track were provided through an exception
        if symbols is None and values is None:
            raise TrackException('Empty track provided [Track=%s, Protein=%s' %(name, str(protein)))


        self._values  = np.array(values)
        self._symbols = symbols
        self._name = name
        self._protein = protein
        
        

    ## ------------------------------------------------------------------------
    ##
    @property
    def name(self):
        return self._name

    ## ------------------------------------------------------------------------
    ##
    @property
    def values(self):
        return self._values

    ## ------------------------------------------------------------------------
    ##
    @property
    def symbols(self):
        return self._symbols

    ## ------------------------------------------------------------------------
    ##
    @property
    def protein(self):
        return self._protein
      
    ## ------------------------------------------------------------------------
    ##      
    def __repr__(self):             
        return "Track [name: %s] associated with protein %s" % (self.name, self.protein)

    ## ------------------------------------------------------------------------
    ##      
    def __len__(self):
        return len(self._protein)

    
    
