SHEPHARD
==============================
#### Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder


### Current major version: 0.1.21 (January 2024)
[//]: # (Badges)
[![Documentation Status](https://readthedocs.org/projects/shephard/badge/?version=latest)](https://readthedocs.org/projects/shephard/badge/?version=latest&style=for-the-badge)


## About
SHEPHARD is a Python toolkit for integrative proteome-wide analysis. It was written by Garrett Ginell and Alex Holehouse.

SHEPHARD enables you to read in protein sequence data and annotate it with different types of sequence annotations (Sites, Domains, and Tracks).
  

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
SHEPHARD is fully released, and the [SHEPHARD paper](http://dx.doi.org/10.1101/2022.09.18.508433) is out in Bioinformatics. Please cite SHEPHARD as:

Ginell, G. M., Flynn, A. J. & Holehouse, A. S. SHEPHARD: a modular and extensible software architecture for analyzing and annotating large protein datasets. Bioinformatics 39, (2023).

## Roadmap
SHEPHARD is the base code for a large body of sequence-based bioinformatic tools developed by the Holehouse lab. These include:

* [metapredict](https://github.com/idptools/metapredict) - high-performance disorder predictor [paper v1](http://dx.doi.org/10.1016/j.bpj.2021.08.039), [paper v2](http://dx.doi.org/10.1101/2022.06.06.494887), [paper v2-ff](http://dx.doi.org/10.1038/s41592-023-02159-5).
* [parrot](https://github.com/idptools/parrot) - a general tool for deep learning of sequence features [paper](http://dx.doi.org/10.7554/eLife.70576)
* [sparrow](https://github.com/idptools/sparrow) - a high-throughput tool for sequence analysis, including the [ALBATROSS networks](http://dx.doi.org/10.1038/s41592-023-02159-5) (*in development*)
* [goose](https://github.com/idptools/goose) - a general-purpose tool for the rational design of disordered protein sequences [paper](http://dx.doi.org/10.1101/2023.10.29.564547).

These tools together form the backbone of our informatics infrastructure, and SHEPHARD will contain direct or indirect API access to each of them (and various other tools).

## Change log
The Changelog below reports on changes as we updated SHEPHARD. Specific types of changes include **BUG FIXES**, **PERFORMANCE UPGRADES**, and **NEW FEATURES**, and these will be tagged as such.


#### Version 0.1.21-patch (May 2024)
* We added the `albatross_api` module to apis, which lets you pass in a Proteome and annotate at either the protein level or domain level all sequence predicted Rg and Re values. Right now this does both but better granularity and tests will be added before the bump to 0.1.22


#### Version 0.1.21 (January 2024)
* BREAKING CHANGE: We renamed shephard.apis.metapredict to shephard.apis.metapredict_api to avoid namespace clashing with metapredict the package. This is of course avoidable by aliasing one/both, but this was poor design. Going forward, we will append _api to the end of api modules.
* Including import of metapredict_api from apis such that `from shephard.apis import metapredict_api` syntax works
* Removed batch_mode as a variable to consider in the metapredict_api functions; size-collect is the only mode supported in metapredict; if this changes we'll revisit things but for now no need to add additional confusion.

#### Version 0.1.20 (December 2023)
* Fixed a minor but where the `shephard.interfaces.si_proteins` interface required proteins to ALREADY be in the proteome which proteins were being added to, which makes no sense, so we removed this constraint. 

#### Version 0.1.19 (November 2023)
* Added version requirement (3.7 to 3.11 inclusive)
* **PERFORMANCE UPGRADE**: Improved how large annotation files are parsed so we ONLY parse lines with unique IDs matching unique IDs in the associated Proteome we're annotating - massive improvement in performance when working with large (10,000 - 100,0000) annotation datasets. This should change nothing on the frontend or any of the behavior other than making SHEPHARD much faster for large datasets
* **PERFORMANCE UPGRADES** Changed some of the error message construction to avoid major overhead when many (1000s of sites) are added (specifically, we previously by default generated an error message that listed out all the sites in a protein when testing for a dictionary type in a Site construction line; this has been removed). 
* Better error handling for interface classes (print only the first 10 errors if many lines are read incorrectly - avoids a situation where the wrong file causes GBs of out text)
* Added explicit tests for all internal Interface classes.
* Added documentation for Protein interface files (as missing previously)


#### Version 0.1.18 (February 2023)
* Added defensive programming for writing sites and domains where if a `domain_type` or `site_type` variable is passed, we check explicitly that it's a list.
* Added ability to write_protein_attributes_from_dictionary (new function in `si_protein_attributes.py`.


#### Version 0.1.17 (September 2022)
* **BUG FIX** Fixed bug in writing domains from list.
* Added import from apis module such that `from shephard import apis` now enables `apis.<module>` to work

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

Copyright (c) 2019-2023, Garrett M. Ginell and Alex S. Holehouse  - [Holehouse lab](http://holehouse.wustl.edu/)

