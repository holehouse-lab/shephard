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

        # update the proteome if this is a novel type of site
        protein.proteome.__update_site_types(self._site_type)


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
    

    def __repr__(self):             
        return "|Site: %i-%s| in %s" % (self.position, self._site_type, self.protein)
    
    
