SHEPHARD
==============================
#### Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder


### Current major version: 0.2.3 (May 2026)
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

Together, these tools form the backbone of our informatics infrastructure, and SHEPHARD provides direct or indirect API access to each of them (and various other tools).

## Change log
The Changelog below reports on changes as we updated SHEPHARD. Specific types of changes include **BUG FIXES**, **PERFORMANCE UPGRADES**, and **NEW FEATURES**, and these will be tagged as such.

#### Version 0.2.3 (May 2026)
* **BUG FIXES**:
  * Fixed a typo in `proteome.py` (`s.position. s.site_type`) that raised an `AttributeError` when copying any `Protein` containing sites into a new `Proteome` (e.g. `Proteome(list_of_Protein_objects)` or `add_proteins()` with `Protein` objects).
  * Fixed a broken type check in `track.py` (`elif(symbols, str):`, an always-true tuple) so symbols tracks are validated correctly and the invalid-input error path is reachable.
  * Fixed an undefined-name error in `track.py` where a non-dict `attribute_dictionary` raised `NameError` instead of a `TrackException`.
  * Fixed an undefined-name error in `Domain.site()` (`domain.start`/`domain.end`) that raised `NameError` instead of a `DomainException` for out-of-range positions.
  * Fixed `interface_tools.is_comment_line()` raising an `IndexError` on blank lines; blank lines are now skipped, so all `si_*` file parsers tolerate blank lines.
  * `tools.sequence_tools.build_mega_string()` now honors the `return_as_list` argument (previously ignored).
  * `tools.site_tools.build_site_density_vector()` now honors the `append_leading_lagging` argument (previously ignored) and raises a clear `ShephardException` when `window_size` exceeds the protein length (previously an `IndexError`).
  * `Proteome.add_proteins()` no longer raises an `IndexError` when passed an empty list.
  * Copying proteins between Proteomes now keys off the track type rather than list truthiness, so an all-zero/empty values track is no longer miscopied as a symbols track.
  * `Domain.site()` now raises a `DomainException` (instead of a raw `KeyError`) when no site exists at an in-range position.
  * Cleaned up the `safe=False` skip path in the `si_sites`, `si_domains`, `si_proteins`, and `si_protein_attributes` interfaces so the `continue` is no longer nested under `if verbose:`.
  * Fixed `Protein.convert_to_valid()`, which was effectively unusable: the internal indexing sentinel was counted in the length check so it *always* raised with the default `safe=True`, and with `copy=False` the cleaned sequence was never written back (silent no-op). It now operates on the user sequence and correctly mutates in place.
  * Fixed `apis.fasta.shephard_fasta_to_proteome()` attribute parsing, which set each round-tripped attribute's value to a copy of its key and produced a spurious empty-string attribute. Attribute name/value round trips now work (and values containing `=` are preserved).
  * Added the missing `random` and `site_tools` imports in `tools.experimental`; both public functions previously raised `NameError` on any call.
  * `Site.remove_attribute()` now raises a `SiteException` instead of a `ProteinException`, consistent with the other annotation objects.
  * `Site.get_track_value()` and `Site.get_track_symbol()` now return `None` (as documented) when the track is missing and `safe=False`, instead of raising a `TypeError` (`None[0]`).
  * Fixed the error path in `Protein.remove_track()` which referenced a non-existent `self.protein` attribute, raising an `AttributeError` instead of the intended `ProteinException` when removing a Track that belongs to a different protein.
  * `Track.attribute()` error message now correctly refers to the `Track` rather than a "protein" (copy-paste fix).
  * `Proteome.__contains__()` now returns `False` (rather than implicitly `None`) when tested against an object that is neither a `str` nor a `Protein`, so `x in proteome` always evaluates to a bool.
  * `Domain.get_track_values()` / `Domain.get_track_symbols()` now use `is None` instead of `== None` for their internal identity checks.
  * Fixed the coverage configuration in `setup.cfg`, which omitted `metapredict/_version.py` (a copy-paste from another project) instead of `shephard/_version.py`.
* **NEW FEATURES**:
  * Added `Proteome.tracks` (property) and `Proteome.get_tracks_by_name()`, mirroring the existing `Proteome.domains` / `Proteome.sites` and `get_domains_by_type()` / `get_sites_by_type()` accessors (previously there was no way to retrieve `Track` objects at the Proteome level).
* **PERFORMANCE UPGRADES**:
  * `Protein.add_domain()` now performs an O(1) uniqueness check against the internal domain dictionary instead of rebuilding and sorting the full domain list on every call (and on every iteration of the `autoname` loop). Loading *D* domains into a protein drops from ~O(D² log D) to ~O(D).
  * `si_domains.add_domains_from_dictionary()`, `si_sites.add_sites_from_dictionary()` and `si_tracks.add_tracks_from_dictionary()` now iterate the (typically small) annotation dictionary and use O(1) protein look-ups, rather than scanning every protein in the Proteome. Annotating *K* proteins in a *P*-protein Proteome drops from O(P) to O(K) — a large speed-up when annotating a subset of a big Proteome (this also makes `Protein.add_domains()`/`build_domain()` no longer O(P)).
  * `Proteome[i]` / slicing no longer materialises a list of every Protein on each access; indexing is now lazy (`P[0]` is ~O(1) rather than O(P)).
* **PACKAGING / INSTALL & CI**:
  * Added the `shephard/py.typed` PEP 561 marker. It was already declared in `package-data` but the file did not exist; SHEPHARD now correctly advertises itself as a typed package.
  * Fixed wheel packaging: the bundled data files (`shephard/data/`, including `test_data/` and `look_and_say.dat`) are now included via `package-data`. Previously they were absent from built wheels, so `shephard.get_data()` pointed at missing files for wheel installs.
  * `pyproject.toml`: removed the unnecessary `numpy` build-system requirement (SHEPHARD is pure Python; `numpy` remains a runtime dependency); bumped `requires-python` to `>=3.8` (3.7 is end-of-life); added trove classifiers, keywords, `[project.urls]`, a second author, a `pytest-cov` test extra, and `[tool.pytest.ini_options] testpaths` so `pytest` discovers the suite from any directory.
  * `setup.cfg`: removed the deprecated `[aliases] test = pytest` (removed in setuptools ≥ 72).
  * Removed dead CI configuration for shut-down services (`.travis.yml`, `.lgtm.yml`, `devtools/travis-ci/`) and added a modern GitHub Actions workflow (`.github/workflows/ci.yml`) running the test suite on Python 3.8–3.12 (Ubuntu + macOS) with coverage upload.
  * Modernized `.readthedocs.yaml` (`ubuntu-22.04`, Python 3.11) and made it install the package itself so autodoc can import `shephard`; clarified `docs/requirements.txt`.
* **TESTS**:
  * Added a comprehensive test suite (`test_comprehensive_*` and `test_bugfix_regressions`) covering the user-facing Proteome, Protein, Domain, Site, Track, interface, tools and FASTA/UniProt APIs, plus an explicit regression test for every bug fixed in this release (~285 new tests).
  * Added a session-scoped `conftest.py` fixture that pins the working directory, making the entire test suite runnable from any directory (~13 legacy `si_*` write tests previously only passed when run from inside the tests directory).
  * Relocated the stray `shephard/test_domains_sites.py` into `shephard/tests/`.
* **REFACTOR**:
  * Converted all `%`-style string formatting in the package to f-strings.

#### Version 0.2.2 (November 2024)
* Updated and fixed `metapredict_api` and `albatross_api` including adding tests
* Defaulted to use metapredict V3
* Restructured organization to use `pyproject.toml`

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

