SHEPHARD
==============================
####Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder


[//]: # (Badges)
[![Travis Build Status](https://travis-ci.com/REPLACE_WITH_OWNER_ACCOUNT/shephard.svg?branch=master)](https://travis-ci.com/REPLACE_WITH_OWNER_ACCOUNT/shephard)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/shephard/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/shephard/branch/master)


## About
SHEPHARD is a Python toolkit for integrative proteome-wide analysis. It was written by Garrett Ginell and Alex Holehouse.

## Installation
Copy and paste into your terminal:

	pip install shephard@git+git://github.com/holehouse-lab/shephard.git

This installs the current bleeding-edge version directly from GitHub (i.e. this repo!).

## Documentation
Incomplete, but available here: [https://shephard.readthedocs.io/en/latest/](https://shephard.readthedocs.io/en/latest/)

## Status
SHEPHARD is currently under an open beta release and we hope to reach public beta in summer 2021.

We are not yet at a release were backwards compatibility is guaranteed - i.e. we might make changes that break compatibility with prior versions if necessary. Be warned! 

If you plan to use SHEPHARD in your science please let [Alex](http://holehouse.wustl.edu/) know, only so he can warn you if bugs are found or backwards compatibility is broken.

## Roadmap
SHEPHARD is the base code for a large body of sequence-based bioinformatic tools developed by the holehouse lab. These include:

* [metapredict](https://github.com/idptools/metapredict) - high-performance disorder predictor
* [parrot](https://github.com/idptools/parrot) - a general tool for deep learning of sequence features
* **sparrow** - a high-throughput tool for sequence analysis (*in development*)
* **goose** - a general purpose tool for the rational design of disordered protein sequences (*in development*)
* [pipit](https://github.com/idptools/PIPIT) - A simple tool for sequential sequence shuffling, as implemented in Langstein *et al.* (to be published).

These tools together form the backbone of our informatics infrastructure, and SHEPHARD will contain direct or indirect API access to each of them (and various other tools).

## Change log
As we approach final release and versions of SHEPHARD are available for distribution, a change log is updated and changes that break backwards compatibility or introduce new features are tagged as minor/major increments. Bug fixes/docs/tests are simply tagged by their git hash.

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

Copyright (c) 2019-2021, Garrett M. Ginell and Alex S. Holehouse  - [Holehouse lab](http://holehouse.wustl.edu/)

