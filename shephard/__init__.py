"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""
import os


# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions


def get_version():
    return __version__

_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_data(path):
    """
    This function 
    """
    return os.path.join(_ROOT, 'data', path)

from .proteome import Proteome
