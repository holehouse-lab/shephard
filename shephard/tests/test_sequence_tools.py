"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

import pytest
import numpy as np

from shephard.tools import sequence_tools

def test_build_mega_string(TS1):
    """
    This function tests the global concatination of sequences 
    from either protien or domain objects
    """
    
    # test passing of protien objects 
    assert sequence_tools.build_mega_string(TS1) == ''.join([i.sequence for i in TS1])
    
    # test passing of domain objects 
    assert sequence_tools.build_mega_string(TS1.domains) == ''.join([i.sequence for i in TS1.domains])

def test_find_string_positions():
    """
    This function tests the global concatination of sequences 
    from either protiens or domains
    """
    s0 = 'AACABBBBACABA'
    s1 = "A"
    s2 = 'CAB'
    s3 = 'Z'
    
    # test multible positions 
    assert sequence_tools.find_string_positions(s1, s0) == [1,2,4,9,11,13]
    
    # test multible characters 
    assert sequence_tools.find_string_positions(s2, s0) == [3,10]
    
    # test nothing 
    assert sequence_tools.find_string_positions(s3, s0) == []