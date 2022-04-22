SHEPHARD
==============================
####Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder


### Current major version: 0.1.5 (April 2022)

[//]: # (Badges)
[![Travis Build Status](https://travis-ci.com/REPLACE_WITH_OWNER_ACCOUNT/shephard.svg?branch=master)](https://travis-ci.com/REPLACE_WITH_OWNER_ACCOUNT/shephard)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/shephard/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/shephard/branch/master)


## About
SHEPHARD is a Python toolkit for integrative proteome-wide analysis. It was written by Garrett Ginell and Alex Holehouse.

## Installation
Copy and paste into your terminal:

	pip install shephard

This installs the current stable release candidate from PyPi.

#### Installation from GitHub

Copy and paste into your terminal:

	pip install shephard@git+git://github.com/holehouse-lab/shephard.git

This installs the current bleeding-edge version directly from GitHub.


## Documentation
Documentation: [https://shephard.readthedocs.io/en/latest/](https://shephard.readthedocs.io/en/latest/)

## Status
As of April 2022 the full release of SHEPHARD is available on PyPI. We are still adding some tweaks and tests but this will be the final version in the SHEPHARD paper.

If you plan to use SHEPHARD in your science please let [Alex](http://holehouse.wustl.edu/) know, only so he can warn you if bugs are found or backwards compatibility is broken.

## Roadmap
SHEPHARD is the base code for a large body of sequence-based bioinformatic tools developed by the Holehouse lab. These include:

* [metapredict](https://github.com/idptools/metapredict) - high-performance disorder predictor
* [parrot](https://github.com/idptools/parrot) - a general tool for deep learning of sequence features
* [pipit](https://github.com/idptools/PIPIT) - A simple tool for sequential sequence shuffling, as implemented in Langstein *et al.* [preprint here](https://www.biorxiv.org/content/10.1101/2022.02.10.480018v1).
* **sparrow** - a high-throughput tool for sequence analysis (*in development*)
* **goose** - a general purpose tool for the rational design of disordered protein sequences (*in development*)


These tools together form the backbone of our informatics infrastructure, and SHEPHARD will contain direct or indirect API access to each of them (and various other tools).

## Change log
As we approach final release and versions of SHEPHARD are available for distribution, a change log is updated and changes that break backwards compatibility or introduce new features are tagged as minor/major increments. Bug fixes/docs/tests are simply tagged by their git hash.

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

