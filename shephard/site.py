"""
Site Class File - ShanEnt Suite 

version: 3

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

from . import general_utilities


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Class that defines a Site in sequence
#
class Site:
    
    def __init__(self, position, site_type, protein, symbol=None, value=None, attributes={}):

        # absolute position in protein associated with the site
        self._position  = int(position)
        
        # reference to the protein object from which the site was taken
        self._protein   = protein

        # the site type identifier
        self._site_type = str(site_type)

        # a numerical value associated with the site
        self._value     = general_utilities.cast_or_none(value, float)

        # a symbol associated with the site
        self._symbol    = general_utilities.cast_or_none(symbol, str)

        self._attributes = attributes

        # update the proteome if this is a novel type of site
        protein.proteome.__update_site_types(self._site_type)


    @property
    def residue(self):
        # note we have to add a +1 offset because sequence is still 
        # 
        return self._protein.sequence[self._position + 1]

    @property
    def position(self):
        return self._position

    @property
    def protein(self):
        return self._protein

    @property
    def site_type(self):
        return self._site_type

    @property
    def symbol(self):
        return self._symbol

    @property
    def value(self):
        return self._value

    @property
    def attributes(self):
        return self._attributes

    def get_local_sequence_context(self, offset = 5):
        self._protein.get_local_sequence_context(self.position, offset)
        
        # compute start/end of context according to the offset
        p1 = max(1, self.position - offset)
        p2 = min(self._protein._len, self.position + offset)

        # note we subtract -1 here from the p1 site to shift from realworld indexing
        # to i0 indexing. However we do not offset p2 because this means the start
        # and end become inclusive
        return self._protein._sequence[p1-1:p2]
    

    def __repr__(self):             
        return "|Site: %i-%s| in %s" % (self.position, self._site_type, self.protein)
    
    
