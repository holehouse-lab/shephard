"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

# Import package, test suite, and other packages as needed
import shephard
import pytest
import sys

def test_shephard_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "shephard" in sys.modules
