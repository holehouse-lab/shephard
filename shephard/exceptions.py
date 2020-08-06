"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis

"""

def print_warning(msg, e):
    print('Warning: %s' %(msg))
    print('%s\n'%(e))

def print_and_raise_error(msg, e):
    print('Error: %s' %(msg))
    raise e

    

class ShephardException(Exception):
    """
    General exception
    """
    pass


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
class SiteException(Exception):
    """
    Exception for the Domain class
    """
    pass


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
class TrackException(Exception):
    """
    Exception for the Domain class
    """
    pass


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
class DomainException(Exception):
    """
    Exception for the Domain class
    """
    pass


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
class ProteinException(Exception):
    """
    Exception for the Proteins class
    """
    pass


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
class ProteomeException(Exception):
    """
    Exception for the Proteome class
    """
    pass


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
class UtilitiesException(Exception):
    """
    Exception for general utility exceptions
    """
    pass


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
class InterfaceException(Exception):
    """
    Exception for general utility exceptions
    """
    pass
