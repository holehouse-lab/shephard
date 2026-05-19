"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

# must import pytest
import pytest

# import additional modules needed for this code
import shephard
from shephard.apis import uniprot  
from shephard import interfaces

# import the domain_tools module which we're going to test
from shephard.tools import domain_tools

