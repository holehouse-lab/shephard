# Changelog

All notable changes to SHEPHARD are recorded here. Changes are tagged by type -
**BUG FIXES**, **NEW FEATURES**, **PERFORMANCE UPGRADES**, **PACKAGING / INSTALL & CI**,
**DOCUMENTATION**, **TESTS**, and **REFACTOR** - so it is easy to see at a glance
whether a release is likely to affect you.

SHEPHARD uses [semantic versioning](https://semver.org/) loosely: patch releases are
bug fixes, minor releases add functionality, and any change that breaks backwards
compatibility is called out explicitly at the top of the relevant entry.

#### Version 0.2.4 (July 2024)
* **BUG FIXES**:
  * `Proteome[i]` could not be indexed from the end: a negative index raised a `KeyError` and a slice with any negative bound (e.g. `P[-2:]`) raised a `ValueError` from `islice`. Indexing now follows the standard Python conventions, so negative indices count back from the last protein, and an out-of-range integer index raises an `IndexError` rather than a `KeyError`.
  * Over-writing an existing Track (`add_track(..., safe=False)`, `build_track_values_from_sequence()`, `build_track_symbols_from_sequence()`, `build_track()`) left the replaced Track counted in the Proteome's book-keeping, so a track name could linger in `Proteome.unique_track_names` (and `track_names_to_track_type`) after every Track with that name had been removed. The same was true of `Protein.add_domain()` and `Proteome.unique_domain_types`. Note that the replacement annotation is now built *before* the old one is discarded, so a failed overwrite leaves the pre-existing annotation intact.
  * `Protein.build_domain()` accepted and documented `safe` and `autoname` but forwarded neither, so `autoname=True` still raised on a duplicate domain and `safe=False` could not be used to over-write one.
  * Track attributes were silently dropped when Proteins were copied into a new Proteome (`Proteome(list_of_Protein_objects)` or `add_proteins()`), because `Protein.add_track()` had no way to set them. `add_track()` now takes an `attributes` dictionary (matching `add_domain()` and `add_site()`), and Track attributes are deep-copied alongside Domain and Site attributes.
  * `Protein.sites` is documented as returning Sites sorted N- to C-terminal but returned them in the order they happened to be added, so `protein.sites` (and therefore `Proteome.sites`, and the order in which Sites were written to file) was unordered whenever Sites were added out of sequence.
  * Three internal book-keeping exceptions in `proteome.py` were missing their `f` prefix, so a genuine book-keeping failure reported a literal `{track_name}`/`{domain_type}`/`{site_type}` placeholder rather than the offending name.
  * `Site.get_track_values()` on a symbols track (and `Site.get_track_symbols()` on a values track) raised a raw `TypeError` (`None` is not subscriptable) rather than an informative `SiteException`. The same was true of `Track.values_region()` and `Track.symbols_region()`, which now raise a `TrackException`. The equivalent `Domain` functions already reported this properly and continue to raise a `DomainException`.
  * A Site with no symbol was written out as the literal string `None` (as an unset value always has been) but read back in as the *string* `'None'` rather than `None`, so Sites did not round-trip cleanly through a Sites file. Symbols are now handled exactly as values are.
* **PERFORMANCE UPGRADES**:
  * `si_protein_attributes.add_protein_attributes_from_dictionary()` now iterates the (typically small) attributes dictionary and uses O(1) protein look-ups rather than scanning every protein in the Proteome, matching the equivalent change made to the domains, sites and tracks interfaces in 0.2.3.
* **PACKAGING / INSTALL & CI**:
  * Added a [tox](https://tox.wiki/) configuration so the full test suite can be run against every supported Python (3.10 - 3.14) with a single `tox` command, mirroring the CI matrix. Individual versions can be run with `tox -e py312`, pytest arguments can be passed through after a `--`, and `tox -e coverage` reproduces the coverage run CI does. The configuration lives in `pyproject.toml` (tox's native TOML format, so it requires tox 4.21 or later) rather than in a separate `tox.ini`, keeping project configuration in one place.
* **DOCUMENTATION**:
  * Corrected the stated Python requirement in the installation docs (3.8 -> 3.10, matching `pyproject.toml`), the release date in the README (May -> July 2026), and the copyright range in the Sphinx config.
  * Documented Proteome indexing/slicing semantics, Track attributes (including the fact that they are in-memory only and are not written to Tracks files), and the `None` symbol/value round-trip convention in Sites files.
  * `si_protein_attributes.write_protein_attributes_from_dictionary()` is a public function but was missing from the interfaces documentation; also fixed the `si_proteins` section, which named the wrong module.
  * Fixed the parameter name in the `metapredict_api.annotate_proteome_with_disorder_track()` docstring (`track_name` -> `name`) and removed a stale claim from `Site.update_site_symbol()` that it updated the Proteome's site book-keeping.
* **TESTS**:
  * Added `test_bugfix_regressions_0_2_4.py`, an explicit regression test for each of the bugs listed above.
  * Documented how to run the suite across Python versions with tox in the code convention documentation.
* **REFACTOR**:
  * Removed unused imports (`numpy` in `protein.py`, `re` in `sequence_utilities.py`, `ProteinException` in `site.py`, and `random`/`site_tools` in `domain_tools.py`).

#### Version 0.2.3 (July 2026)
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
  * **Exception hierarchy**: every SHEPHARD exception (`SiteException`, `TrackException`, `DomainException`, `ProteinException`, `ProteomeException`, `UtilitiesException`, `InterfaceException`, `APIException`) subclassed bare `Exception` rather than `ShephardException`, so `except ShephardException:` caught essentially nothing raised by the library. All of them now inherit from `ShephardException` as intended. Also corrected `TrackException`'s docstring (it described the Domain class) and the three exceptions that shared a copy-pasted "general utility exceptions" description.
  * `Domain.get_sites_by_type()` accepted and documented a `return_list` argument but never forwarded it, so `return_list=True` silently returned a dictionary.
  * Fixed the `extend_ends` block in `tools.domain_tools.build_domains_from_track_values()`, which could not work as written: it used a float as a list index (`len(B)/2`, a `TypeError`) and compared a list *slice* against an integer (always `False`, so the C-terminal extension never happened). The two ends are now handled symmetrically, and a zero-width extension is a no-op rather than an error.
  * `tools.experimental.get_site_density_in_domain_normalized_by_protein()` returned the *minimum* enrichment value for a domain that was maximally enriched (sites in the domain, none in the sampled background) and the maximum for a depleted one — i.e. the two extremes were inverted. It also sampled background regions using 0-indexed positions (including position 0, which is not a valid protein position) while querying with SHEPHARD's 1-indexed inclusive convention.
  * `tools.experimental.get_site_density_percentile_normalized_by_protein()` searched the density distribution using exact float equality, so it fell off the end of the function and returned `None` whenever no exact match existed. Both functions in this module now have proper docstrings.
  * **Track file format**: track lines were written with a trailing delimiter. This is invisible with the default tab delimiter, but with any other delimiter it produced an empty final field that could not be read back (`float('')` in values mode, or a spurious extra symbol in symbols mode). Track files are now written without a trailing delimiter, matching the domain and site writers. *Note this changes the byte content of written track files by one byte per line; files written by older versions still read correctly.*
  * `si_tracks`: an invalid `mode` raised from inside the per-line `try`, so `skip_bad=True` swallowed it — passing a bad mode produced ten confusing warnings and an empty dictionary instead of an error. The mode is now validated before parsing begins.
  * `si_tracks.write_tracks_from_list()` raised an `IndexError` on an empty list (the domain and site equivalents handle this fine); it now writes an empty file. `write_tracks()` also leaked its file handle if writing raised, and would call `open(None)` if given neither a filename nor a file handle; both now use a context manager and the argument combination is validated.
  * `si_protein_attributes.write_protein_attributes_from_dictionary()` could not consume the dictionary the corresponding parser produces: the reader returns `{unique_ID: [dict, ...]}` but the writer assumed `{unique_ID: dict}` and raised an `AttributeError` on the list. It now accepts both shapes.
  * **`skip_bad` now covers attribute parsing**: in every `si_*` module the `key:value` parsing sat *outside* the guarded `try`, so a malformed attribute raised an `InterfaceException` even with `skip_bad=True`. The documented hard cap on skipped lines (`MAX_BAD_COUNT`, 10) is also now documented in every `skip_bad` docstring — previously `skip_bad=True` would still raise on the eleventh bad line with no indication that this was expected.
  * **Attribute keys are now sanitized on write**: all the `si_*` writers cleaned attribute *values* but interpolated *keys* raw, so a key containing a delimiter or a colon produced a file that could not be read back. Relatedly, `interface_tools.full_clean_string()` hardcoded the tab character regardless of the delimiter the writer was actually using; it now takes the delimiter.
  * `si_sites.add_sites_from_dictionary()` never validated that its first argument was a `Proteome` (every sibling `*_from_dictionary()` function does), and `SiteException` was imported but never caught, so it escaped even with `safe=False`. Several error messages in `si_sites` and `si_domains` also re-cast the offending value (e.g. `int(start)`) inside the exception handler, which could raise a second error while reporting the first.
  * `apis.metapredict_api`: `disorder_threshold` was accepted and documented by both domain-annotating functions but never passed to metapredict, so a non-default threshold was silently ignored. The default is now `None` (metapredict's version-appropriate default, which is what these functions effectively used before) and the value is passed through.
  * `apis.metapredict_api` ran an actual disorder prediction at import time as a version probe, so `import shephard.apis` loaded a model as a side effect. The version check is now lazy. A redundant `try/except` that retried an identical call was also removed.
  * **Missing optional dependencies now fail clearly**: both `metapredict_api` and `albatross_api` only *printed* on a failed import, leaving the imported name undefined so the first real call died with an opaque `NameError`. Both now raise an informative `APIException` at call time.
  * `apis.albatross_api.annotate_domains_with_dimensions()` keyed predictions on `unique_ID + '_' + domain_name`, which is ambiguous if a unique_ID itself contains an underscore (silently mis-assigning or dropping predictions). Keys are now unambiguous. The unused `batch_mode` argument is retained but documented as a deprecated no-op.
  * `apis.fasta.fasta_to_proteome()`: a guard testing `use_header_as_unique_ID is None` on a boolean was always `False`, and auto-generated unique_IDs were `int` while every other code path produces `str` (so downstream string look-ups missed).
  * `apis.fasta.proteome_to_fasta()` raised an `AttributeError` for any non-string attribute value (attributes are explicitly documented as arbitrary), and round-tripping a Proteome through a FASTA file with `include_attributes_in_header=True` appended a trailing space to every protein name each time.
  * `tools.attribute_tools.cast_attributes()` caught only `ValueError`, so `skip_failed=True` did not skip the most common failure (`TypeError` from e.g. `float(None)`).
* **NEW FEATURES**:
  * Added `Proteome.tracks` (property) and `Proteome.get_tracks_by_name()`, mirroring the existing `Proteome.domains` / `Proteome.sites` and `get_domains_by_type()` / `get_sites_by_type()` accessors (previously there was no way to retrieve `Track` objects at the Proteome level).
  * **Live UniProt annotation**: `apis.uniprot` can now pull annotations directly from the UniProt REST API and convert them into SHEPHARD annotations, via `annotate_protein_with_uniprot()` (a single protein) and `annotate_proteome_with_uniprot()` (a whole Proteome, batching accessions so a few thousand proteins cost a few tens of API calls rather than one call each). Features that span a region become Domains, single-residue features become Sites, and disulfide bonds/cross-links become a pair of Sites recording each other's position. Regions covered by an experimental structure (from the PDB cross-references) can be added as Domains carrying the PDB ID, method, and resolution, and secondary structure and structural coverage can be added as Tracks. Every annotation carries the UniProt feature type, description, evidence codes, and source accession as attributes. Which annotations are fetched is controlled per-class by keyword, and only the data required is downloaded. Because UniProt annotations are positional, the local sequence is verified against the UniProt sequence by default and proteins whose sequences have drifted are refused rather than silently mis-annotated. Requires network access, but adds no new dependencies. See the [API documentation](https://shephard.readthedocs.io/en/latest/apis.html) for details and worked examples.
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
  * **Raised the minimum supported Python to 3.10.** Python 3.8 and 3.9 are both end-of-life, and on those versions pip back-solves to a numpy old enough that we never test against it. Classifiers now cover Python 3.10-3.14, and the CI matrix was updated to match (3.10-3.14 on Ubuntu, plus 3.12 on macOS).
  * Adopted [PEP 639](https://peps.python.org/pep-0639/) license metadata (`license = "MIT"` plus `license-files`), replacing the deprecated license table and the `License :: OSI Approved` classifier. Built distributions now carry a `License-Expression` field and ship `LICENSE` in the wheel. This requires `setuptools>=77`, so the build requirement was bumped accordingly.
  * Added `Development Status :: 5 - Production/Stable` (SHEPHARD has been released and in use for several years), `Natural Language :: English`, `Topic :: Scientific/Engineering :: Chemistry`, and `Typing :: Typed` classifiers; expanded the keyword list.
  * Added `Changelog` and `Research paper` entries to `[project.urls]`, so both are linked directly from the PyPI sidebar.
  * `MANIFEST.in` now ships `CHANGELOG.md` and `CITATION.cff` in the source distribution.
  * Moved the coverage configuration from `setup.cfg` into `pyproject.toml` (`[tool.coverage.run]`), so pytest, coverage and packaging config all live in one place. `setup.cfg` now holds only the flake8 and yapf style settings, since flake8 has no `pyproject.toml` support.
* **DOCUMENTATION**:
  * Documented the new UniProt REST API in `docs/apis.rst`, covering how UniProt features map onto SHEPHARD annotations, the metadata attached to each annotation, a full table of the annotation groups and their defaults, sequence verification, the summary dictionary and `safe` semantics, and batching behaviour.
  * Added three worked examples to `docs/examples.rst`: annotating a proteome directly from UniProt, combining UniProt PTM Sites with metapredict IDR Domains, and using the structural coverage Track to find regions with no experimental structure. All three are runnable against the bundled test set.
  * Added `demo/uniprot/uniprot_api_examples.ipynb`, an executed notebook walking through the UniProt API end to end - annotation groups, single-protein and batch annotation, disulfide bonds, provenance and evidence codes, experimental structures and Tracks, error handling, saving annotations, and a worked enrichment question.
  * Moved the change log out of `README.md` into `CHANGELOG.md`, and added the usual status badges and a citation section to the README.
  * Added `CITATION.cff`, so GitHub renders a "Cite this repository" entry pointing at the SHEPHARD paper.
  * Fixed the install command in `README.md`, which still used the `git://` protocol GitHub disabled in 2022, and updated the paper citation to the published DOI (`10.1093/bioinformatics/btad488`) rather than the preprint.
* **TESTS**:
  * Added a comprehensive test suite (`test_comprehensive_*` and `test_bugfix_regressions`) covering the user-facing Proteome, Protein, Domain, Site, Track, interface, tools and FASTA/UniProt APIs, plus an explicit regression test for every bug fixed in this release (~285 new tests).
  * Added a session-scoped `conftest.py` fixture that pins the working directory, making the entire test suite runnable from any directory (~13 legacy `si_*` write tests previously only passed when run from inside the tests directory).
  * Relocated the stray `shephard/test_domains_sites.py` into `shephard/tests/`.
  * Added `test_bugfix_regressions_0_2_3.py`, an explicit regression test for each of the bugs listed above, and `test_api_uniprot_annotation.py`, a fully offline suite for the new UniProt REST API (parsing, annotation, batching, and every error path, plus the HTTP layer's retry/backoff behaviour). Live tests against the real UniProt API are in `test_api_uniprot_live.py` and are skipped unless `SHEPHARD_TEST_NETWORK=1`, so CI stays hermetic.
  * Fixed three test modules that defined two test functions with the same name, where the second silently shadowed the first so it never ran (`test_api_albatross.py`, `test_si_tracks.py`, `test_proteome.py`). One shadowed test was a strict duplicate and was removed; the rest are now named distinctly and run.
  * `test_api_albatross.py` now skips (rather than failing with a `NameError`) when `sparrow` is not installed, matching the metapredict test. This was a live CI failure, since CI installs neither predictor.
  * The `TS1_domains2_sites` fixture in `conftest.py` was both session-scoped and `autouse`, so a single mutable `Proteome` was built for every test in the session and shared between the tests that used it. It is now function-scoped like every other fixture.
  * Removed three dead files: `shephard/data/test_data/tst.py` (scratch code importing APIs that no longer exist — and, because the `data/` directory is shipped as package data, installed into every wheel), `shephard/tests/test_logic.py` (empty), and `shephard/tests/test_domains_sites.py` (imports only, no test functions).
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
