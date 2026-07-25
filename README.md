SHEPHARD
==============================
#### Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

### Current major version: 0.2.3 (July 2026)

[//]: # (Badges)
[![CI](https://github.com/holehouse-lab/shephard/actions/workflows/ci.yml/badge.svg)](https://github.com/holehouse-lab/shephard/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/shephard/badge/?version=latest)](https://shephard.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/holehouse-lab/shephard/branch/master/graph/badge.svg)](https://codecov.io/gh/holehouse-lab/shephard)
[![PyPI](https://img.shields.io/pypi/v/shephard.svg)](https://pypi.org/project/shephard/)
[![Python versions](https://img.shields.io/pypi/pyversions/shephard.svg)](https://pypi.org/project/shephard/)
[![Downloads](https://static.pepy.tech/badge/shephard)](https://pepy.tech/project/shephard)
[![License: MIT](https://img.shields.io/pypi/l/shephard.svg)](https://github.com/holehouse-lab/shephard/blob/master/LICENSE)
[![Last commit](https://img.shields.io/github/last-commit/holehouse-lab/shephard.svg)](https://github.com/holehouse-lab/shephard/commits/master)
[![DOI](https://img.shields.io/badge/DOI-10.1093%2Fbioinformatics%2Fbtad488-blue.svg)](https://doi.org/10.1093/bioinformatics/btad488)


## About
SHEPHARD is a Python toolkit for integrative proteome-wide analysis. It was written by Garrett Ginell and Alex Holehouse.

SHEPHARD enables you to read in protein sequence data and annotate it with different types of sequence annotations (Sites, Domains, and Tracks).
  

## Installation
Copy and paste into your terminal:

	pip install shephard

This installs the current stable release candidate from PyPi. SHEPHARD is pure Python and requires Python 3.10 or higher; the only runtime dependencies are [numpy](https://numpy.org/) and [protfasta](https://protfasta.readthedocs.io/), both installed automatically.

#### Installation from GitHub

Copy and paste into your terminal:

	pip install shephard@git+https://github.com/holehouse-lab/shephard.git

This installs the current bleeding-edge version directly from GitHub.


## Documentation
Online documentation for SHEPHARD can be found here:

[https://shephard.readthedocs.io/en/latest/](https://shephard.readthedocs.io/en/latest/)

## Tutorial Examples
Examples and Google Colab tutorials can be found here: 

[https://github.com/holehouse-lab/shephard-colab](https://github.com/holehouse-lab/shephard-colab)

Runnable demo notebooks also ship with the package in [`demo/`](demo/), including a [walkthrough of the UniProt API](demo/uniprot/uniprot_api_examples.ipynb).

## Status
SHEPHARD is fully released and in active use, and the SHEPHARD paper is out in *Bioinformatics*.

## Citation
If you use SHEPHARD in your work, please cite:

> Ginell, G. M., Flynn, A. J. & Holehouse, A. S. SHEPHARD: a modular and extensible software architecture for analyzing and annotating large protein datasets. *Bioinformatics* **39**, btad488 (2023). [https://doi.org/10.1093/bioinformatics/btad488](https://doi.org/10.1093/bioinformatics/btad488)

In BibTeX:

```bibtex
@article{ginell2023shephard,
  author  = {Ginell, Garrett M. and Flynn, Aidan J. and Holehouse, Alex S.},
  title   = {{SHEPHARD}: a modular and extensible software architecture for analyzing and annotating large protein datasets},
  journal = {Bioinformatics},
  volume  = {39},
  number  = {8},
  pages   = {btad488},
  year    = {2023},
  doi     = {10.1093/bioinformatics/btad488}
}
```

This repository also includes a [`CITATION.cff`](CITATION.cff) file, so GitHub's "Cite this repository" button will generate the citation for you in whichever format you need.

## Roadmap
SHEPHARD is the base code for a large body of sequence-based bioinformatic tools developed by the Holehouse lab. These include:

* [metapredict](https://github.com/idptools/metapredict) - high-performance disorder predictor [paper v1](http://dx.doi.org/10.1016/j.bpj.2021.08.039), [paper v2](http://dx.doi.org/10.1101/2022.06.06.494887), [paper v2-ff](http://dx.doi.org/10.1038/s41592-023-02159-5).
* [parrot](https://github.com/idptools/parrot) - a general tool for deep learning of sequence features [paper](http://dx.doi.org/10.7554/eLife.70576)
* [sparrow](https://github.com/idptools/sparrow) - a high-throughput tool for sequence analysis, including the [ALBATROSS networks](http://dx.doi.org/10.1038/s41592-023-02159-5) (*in development*)
* [goose](https://github.com/idptools/goose) - a general-purpose tool for the rational design of disordered protein sequences [paper](http://dx.doi.org/10.1101/2023.10.29.564547).

Together, these tools form the backbone of our informatics infrastructure, and SHEPHARD provides direct or indirect API access to each of them (and various other tools).

## Change log

The full change log lives in [CHANGELOG.md](CHANGELOG.md).

### Copyright

Copyright (c) 2019-2026, Garrett M. Ginell and Alex S. Holehouse  - [Holehouse lab](http://holehouse.wustl.edu/)
