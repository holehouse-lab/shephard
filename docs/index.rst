.. shephard documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SHEPHARD
=========================================================

**SHEPHARD** (**S**\ equence-based **H**\ ierarchical and **E**\ xtendable **P**\ latform for **H**\ igh-throughput **A**\ nalysis of **R**\ egions of **D**\ isorder) is a Python-based software framework for reading, annotating, and analyzing large protein datasets. Its objective is to make it easy to perform reproducible, error-free analyses of large protein datasets with arbitrarily complex sequence annotations.

The emergence of high-throughput experiments and high-resolution computational predictions has led to an explosion in the quality and volume of protein sequence annotations at proteomic scales. Unfortunately, sanity-checking, integrating, and analyzing complex sequence annotations remains logistically challenging and introduces a major barrier to entry for even superficial integrative bioinformatics. SHEPHARD addresses this technical burden by combining an object-oriented hierarchical data structure with database-like features, enabling programmatic annotation, integration, and analysis of complex datatypes at proteome-wide scales (millions of unique annotations).

The general use-case for SHEPHARD would include scenarios such as:

1. Wanting to ask large-scale statistical questions about big protein datasets.
2. Asking how different types of sequence or protein annotations relate to one another.
3. Easily associating experimentally-generated data with extant sequence information.
4. Sharing analyses and annotated data with the broader community in a reproducible, accessible way.

SHEPHARD is deliberately **not** an analysis library: it does not provide built-in analysis routines. Instead it functions as the underlying, internally-consistent backend for sequence- or structure-centric analysis pipelines, taking care of the data-cleaning, data-structure, and I/O complexity so that the user can focus on the science.


Installation
--------------

SHEPHARD is pure Python (no compiled extensions) and supports Python 3.8 and above. It is distributed via the Python Package Index (PyPI), so the current public release can be installed using::

	pip install shephard

Alternatively, you can install the current bleeding-edge version directly from GitHub using

.. code-block:: bash

	pip install shephard@git+https://github.com/holehouse-lab/shephard.git

The only runtime dependencies are `numpy <https://numpy.org/>`_ and `protfasta <https://protfasta.readthedocs.io/>`_, both of which are installed automatically.

To test if the installation has worked, open the Python interpreter or a Jupyter notebook and run

.. code-block:: python

	import shephard

If this works, you should be good to go!

A note on Python environments
.................................

We strongly recommend installing into an isolated environment (e.g. a :code:`conda` environment or a :code:`venv`/:code:`virtualenv`). Isolated environments let you pin a specific Python version and set of packages, keeping SHEPHARD insulated from your system Python. For more information there are many `examples of conda introduction tutorials online, like this one here <https://astrobiomike.github.io/unix/conda-intro>`_.

Demos and examples
--------------------
Once installed, SHEPHARD makes it very easy to work with large protein datasets. As an example, we have a collection of basic notebooks showing functionality on the `SHEPHARD colab repository <https://github.com/holehouse-lab/shephard-colab>`_. Precomputed annotations for many model proteomes (including proteome-wide annotations and third-party studies) are available at the `shephard-data repository <https://github.com/holehouse-lab/shephard-data>`_.


Colab notebook
--------------------
If you don't want to install SHEPHARD on your local computer, we also provide several Google Colab notebooks with SHEPHARD pre-installed and a collection of precomputed proteome-wide annotations. Check out our `colab notebook repository on GitHub <https://github.com/holehouse-lab/shephard-colab/>`_.


Citing SHEPHARD
----------------
If you use SHEPHARD in your work, please cite:

	Ginell, G. M., Flynn, A. J. & Holehouse, A. S. SHEPHARD: a modular and extensible software architecture for analyzing and annotating large protein datasets. *Bioinformatics* **39**, btad488 (2023). https://doi.org/10.1093/bioinformatics/btad488


About
---------
SHEPHARD was developed and written by Garrett Ginell and Alex Holehouse in the Holehouse lab. For issues, bugs, or feature requests `please raise an issue on GitHub <https://github.com/holehouse-lab/shephard/issues/>`_.


Documentation index
----------------------


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   overview
   getting_started
   proteome
   protein
   domain
   site
   track
   shephard_file_types
   interfaces
   tools
   apis
   examples
   code_convention


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
