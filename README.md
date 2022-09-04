SHEPHARD
==============================
#### Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder


### Current major version: 0.1.16 (September 2022)

[//]: # (Badges)
[![Travis Build Status](https://travis-ci.com/REPLACE_WITH_OWNER_ACCOUNT/shephard.svg?branch=master)](https://travis-ci.com/REPLACE_WITH_OWNER_ACCOUNT/shephard)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/shephard/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/shephard/branch/master)


## About
SHEPHARD is a Python toolkit for integrative proteome-wide analysis. It was written by Garrett Ginell and Alex Holehouse.

SHEPHARD enables you to read in protein sequence data and annotate it with different types of sequence annotations (Sites, Domains, and Tracks). As an example


## Installation
Copy and paste into your terminal:

	pip install shephard

This installs the current stable release candidate from PyPi.

#### Installation from GitHub

Copy and paste into your terminal:

	pip install shephard@git+git://github.com/holehouse-lab/shephard.git

This installs the current bleeding-edge version directly from GitHub.


## Documentation
Online documentation for SHEPHARD can be found here:

[https://shephard.readthedocs.io/en/latest/](https://shephard.readthedocs.io/en/latest/)

## Tutorial Examples
Examples and Google Colab tutorials can be found here: 

[https://github.com/holehouse-lab/shephard-colab](https://github.com/holehouse-lab/shephard-colab)

## Status
SHEPHARD is fully released, and the SHEPHARD preprint is forthcoming. 

## Roadmap
SHEPHARD is the base code for a large body of sequence-based bioinformatic tools developed by the Holehouse lab. These include:

* [metapredict](https://github.com/idptools/metapredict) - high-performance disorder predictor
* [parrot](https://github.com/idptools/parrot) - a general tool for deep learning of sequence features
* [pipit](https://github.com/idptools/PIPIT) - A simple tool for sequential sequence shuffling, as implemented in Langstein *et al.* [preprint here](https://www.biorxiv.org/content/10.1101/2022.02.10.480018v1).
* [sparrow](https://github.com/idptools/sparrow) - a high-throughput tool for sequence analysis (*in development*)
* [goose](https://github.com/idptools/goose) - a general purpose tool for the rational design of disordered protein sequences (*in development*)


These tools together form the backbone of our informatics infrastructure, and SHEPHARD will contain direct or indirect API access to each of them (and various other tools).

## Change log
As we approach final release and versions of SHEPHARD are available for distribution, a change log is updated and changes that break backwards compatibility or introduce new features are tagged as minor/major increments. Bug fixes/docs/tests are simply tagged by their git hash.


#### Version 0.1.16 (September 2022)
* Update for [PyPI update](https://pypi.org/project/shephard/)
* Improved documentation ahead of final release (including tools docs).
* Added ability to return sites as lists for all site acquisition functions in proteins and domains.
* Added much more detailed tests for site acquisition functions 


#### Version 0.1.15 (September 2022)
* Update for [PyPI update](https://pypi.org/project/shephard/)

#### Version 0.1.10 (September 2022)
* Major update 
* Lots of new tests 
* Enable sites to read/write if values = None without throwing an exception
* Fixed bug in writing sites from list
* **BREAKING CHANGE**: Changed `shephard.protein.get_residue()` to `shephard.protein.residue()`, inkeeping with style for other getter functions


#### Version 0.1.9 (September 2022)
* Major update
* Lots of new tests
* Added ability to write lists of sites and tracks (as we can with domains)
* Refactoring of interface writing code
* Added explicitly checks for domain, site, and track types when writing from lists of these objects
* Added `Track.symbol()` and `Track.value()` functions to extract a single symbol or value at a specific position.
* Updated documentation to include these new functions
* Updated tests to encompass new features
* Fixed bugs in exception handling
* **BREAKING CHANGE**: Changed `shephard.interfaces.si_tracks.write_track()` to `shephard.interfaces.si_tracks.write_tracks()` (i.e. plural) to match names from other functions


#### Version 0.1.8 (August 2022)
* Bug fix in `domain_tools.py` for identifying overlap between two domains
* Fixed inconsistencies in writing domains that led to trailing whitespace
* Fixed bugs in exception throwing code
* More tests

#### Version 0.1.7 (April 2022)
* Improved documentation
* Added domain_to_track() function in tools.track_tools


#### Version 0.1.5 (April 2022)
* First version released to PyPI

#### Version 0.1.4 (Feb 2022)
* Added ability to remove Tracks, Sites and Domains from a Protein objects
* Track number of unique domains, sites, and tracks rather than just their presence/absence
* Updated Track writing
* Added Tracks MUST be either symbolic or values-based but cannot be both
* 

#### Version 0.1.3.1 (May 2021)
* Various bug fixes
* Improved performance 
* Updated interfaces for reading/writing different types of files
* Major updates to internal docs
* This release should be considered largely stable, although docs are lacking
* Expanded the test suite


### Version 0.1.2.1 (August 2020)
**WARNING**: This version breaks backwards compatibility with prior versions!

* `protein.get_domains_by_type()` now returns a list of domains instead of a dictionary. This helps bring consistency to how domains are retrieved and moves us away from dictionary returning.
* Various internal updates 

### Copyright

Copyright (c) 2019-2022, Garrett M. Ginell and Alex S. Holehouse  - [Holehouse lab](http://holehouse.wustl.edu/)

