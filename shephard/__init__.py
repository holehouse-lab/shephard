"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""
import os
from shephard.proteome import Proteome


# Generate _version.py if missing and in the Read the Docs environment
if os.getenv("READTHEDOCS") == "True" and not os.path.isfile('../shephard/_version.py'):   
    import versioningit            
    __version__ = versioningit.get_version('../')
else:
    from ._version import __version__

def get_version():
    return __version__

_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_data(path):
    """
    This function 
    """
    return os.path.join(_ROOT, 'data', path)


