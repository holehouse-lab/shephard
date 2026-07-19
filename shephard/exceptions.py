"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis

"""

def print_warning(msg):
    """
    Function that prints a warning message.
    """
    
    print(f'WARNING: {msg}')


def print_and_raise_error(msg, e):    
    print(f'ERROR: {msg}')
    raise e

    

class ShephardException(Exception):
    """
    Base exception for all SHEPHARD exceptions. Every other exception
    defined here inherits from this, so

        except ShephardException:

    will catch any error raised by SHEPHARD.
    """
    pass


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
class SiteException(ShephardException):
    """
    Exception for the Site class
    """
    pass

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
class TrackException(ShephardException):
    """
    Exception for the Track class
    """
    pass


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
class DomainException(ShephardException):
    """
    Exception for the Domain class
    """
    pass


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
class ProteinException(ShephardException):
    """
    Exception for the Proteins class
    """
    pass


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
class ProteomeException(ShephardException):
    """
    Exception for the Proteome class
    """
    pass


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
class UtilitiesException(ShephardException):
    """
    Exception for the general_utilities and sequence_utilities modules
    """
    pass


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
class InterfaceException(ShephardException):
    """
    Exception for the interfaces (si_*) modules
    """
    pass


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
class APIException(ShephardException):
    """
    Exception for the apis modules
    """
    pass
