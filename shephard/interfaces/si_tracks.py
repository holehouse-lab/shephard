from .interface_exceptions import InterfaceException
from . import interface_tools 

class TracksInterface:

    def __init__(self, filename, delimiter='\t', mode='values', safe=True):
        """
        Expect files of the followin format:

        track_name, unique_ID, val1, val2, ...., valn where n = length of protein


        """

        if mode not in ['values','symbols']:
            raise InterfaceException("When parsing a track file mode must be either 'value' or 'symbol'")
            

        with open(filename,'r') as fh:
            content = fh.readlines()

        ID2track={}

        linecount=0
        for line in content:

            linecount=linecount+1

            sline = line.split(delimiter)
                        
            data_vector=[]
            try:

                # extract track name and unique_id
                track_name = sline[0].strip()
                unique_ID = sline[1].strip()

                # parse track values or symbols
                for i in sline[2:]:
                    if mode == 'value':
                        data_vector.append(float(i.strip()))
                    else:
                        data_vector.append(i.strip())
                        
                if unique_ID not in ID2track:
                    ID2track[unique_ID] = [[track_name, data_vector]]
                else:
                    ID2track[unique_ID].append([track_name, data_vector])

            except Exception:

                # should update this to also display the actual error...
                raise InterfaceException('Failed parsing file [%s] on line [%i]... line printed below:\n%s'%(filename, linecount, line))

        self.data = ID2track



def read_in_tracks(proteome, filename, delimiter='\t', mode='values', safe=True):
    """
    Function that takes a 
    Parameters
    ----------
    proteome : anything
        
    


    safe : bool (default = True)
        Flag which if true with throw an exception of a track with the same name already exists
        
    """

        
    # check first argument is a proteome
    interface_tools.check_proteome(proteome, 'read_in_tracks (si_tracks)')

    track_interface = TracksInterface(filename, delimiter, mode)

    for protein in proteome:
        if protein.unique_ID in track_interface.data:
            for track in track_interface.data[protein.unique_ID]:

                name = track[0]
                if mode == 'values':
                    protein.add_track(name, values=track[1], safe=safe)
                else:
                    protein.add_track(name, symbols=track[1], safe=safe)
                    
                    
                
                    
