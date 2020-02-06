"""
shephard
Sequence-based Hierachical and Extendable Platform for High-throughput Analysis of Region of Disorder
"""
import os

# Add imports here
#from .shephard import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions


_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_data(path):
    return os.path.join(_ROOT, 'data', path)
