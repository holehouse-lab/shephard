"""
Unit and regression test for the shephard package.
"""

# Import package, test suite, and other packages as needed
import shephard
import pytest
import sys

def test_shephard_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "shephard" in sys.modules
