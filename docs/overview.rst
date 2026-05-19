Overview
===============

SHEPHARD (**S**\ equence-based **H**\ ierarchical and **E**\ xtendable **P**\ latform for **H**\ igh-throughput **A**\ nalysis of **R**\ egions of **D**\ isorder) is a modular framework for performing programmatic analyses of large protein datasets. Specifically, SHEPHARD makes it easy to:

* Read in large protein datasets
* Annotate with pre-defined annotations, or compute and annotate with novel analysis routines
* Perform large-scale analysis of annotations associated with those protein datasets
* Share both analyses and annotated data with the broader community in a straightforward way

SHEPHARD was designed to be a simple, efficient, and easy-to-use data management package. It does not offer built-in analysis routines, but instead was built to function as the underlying backend for sequence- or structure-centric analysis pipelines.

It does not itself perform any kind of analysis, but instead provides all of the nitty-gritty tooling that is necessary but complicated for doing large-scale analysis of protein datasets. By hiding all of the logistical complexity, a user can quickly ask complex questions. The goal here is to make it trivial to ask what might otherwise be computationally complex things to address.

SHEPHARD addresses three main challenges in the context of data integration and analysis:

1. It enables a clear and syntactically simple programming interface for asking broad, integrative statistical questions across large datasets.
2. It makes it easy to read annotations from external sources and write annotations to files, making it easy to integrate into existing workflows.
3. It is provided as a locally-installable Python package *and* as Google Colab notebooks that ship with preloaded annotations for the human proteome, so it serves both seasoned bioinformaticians and novice users.


-------------------------
SHEPHARD architecture
-------------------------
SHEPHARD stores information in a bidirectional, object-oriented hierarchical format. This hierarchy is built from a few key building blocks: **Proteomes**, **Proteins**, **Domains**, **Sites**, **Tracks**, and **attributes**.

* A **Proteome** is the base container. It holds *n* **Proteins** and can itself carry arbitrary key-value **attributes**.
* A **Protein** is the major unit of data storage. Each Protein is associated with exactly one amino acid sequence and can be annotated with *m* **Domains**, *p* **Tracks**, *q* **Sites**, and *r* **attributes**.
* A **Domain** is a contiguous subregion of the protein (defined by a start and end position and a free-form ``domain_type``).
* A **Site** is a residue-specific annotation at a single position; multiple Sites can coexist at the same position, and each Site can carry a symbol and/or a numerical value.
* A **Track** is a numerical (*values*) or symbolic (*symbols*) vector with a one-to-one mapping to each residue in the protein sequence.
* **Attributes** are arbitrary key-value pairs that can be attached to a Proteome, Protein, Domain, Site, or Track.

The hierarchy is *bidirectional*: every annotation knows the Protein it belongs to (via its ``.protein`` accessor) and every Protein knows its Proteome (via ``.proteome``), so you can always traverse from an annotation back up to its parent and vice versa. This makes it trivial to, for example, take a Site, ask which Domains it falls inside, and then read out a Track value at that position.

SHEPHARD makes it easy to define your own Domains, Tracks, and Sites, and to read these into and out of files. It does this by defining specific, well-documented tab-separated file formats (described in the :doc:`shephard_file_types` and :doc:`interfaces` documentation pages). Sequences are read in via FASTA format, and SHEPHARD piggybacks on the robust FASTA parser `protfasta <https://protfasta.readthedocs.io/>`_ developed by the Holehouse lab.


-------------------------------------
Loose coupling and software longevity
-------------------------------------

Because SHEPHARD interacts with the 'outside world' (other types of data) via **interfaces** and **APIs** that require a specific file format (or programmatic input format), it maintains a loose coupling between the internal architecture and any other tool, pipeline, or dataset. The interfaces ensure a contract between SHEPHARD and the outside world, which means the user can write robust code using SHEPHARD knowing that they are insulated from changes to those external tools.

This design strategy is deliberate: SHEPHARD does not make assumptions about the outside world, avoiding the challenges introduced by changing dependencies or deprecation of external libraries. The practical consequence is **software longevity** — analysis pipelines written today can be re-run and re-evaluated years from now, against changing datasets, without modifying the pipeline itself. Customizable and reproducible analysis becomes extremely easy, and the same code can be used to perform analogous analyses for very different questions.

In addition to loose coupling, SHEPHARD performs behind-the-scenes consistency and sanity checking when annotations are added (for example, verifying a Track is the same length as the sequence, that Domain boundaries fall within the protein, or that a unique_ID is genuinely unique). This helps catch simple errors, malformatted data, or inconsistent annotations early.


---------------------------------
General workflow in SHEPHARD
---------------------------------

The general workflow we envisage for people using SHEPHARD is as follows:

1. Read in a set of sequences from a FASTA file and generate a new ``Proteome`` using one of the :code:`apis` functions (described elsewhere in the documentation).

2. Load in some annotations (Domains, Sites, Tracks, or protein attributes) via the :code:`interfaces` functions, or compute them on the fly using the :code:`tools` modules / your own functions.

3. Perform analysis in a way that queries the relationship between sequence and those different annotations.

4. (Optionally) Write the resulting annotations back out to SHEPHARD files so they can be shared as plain, human-readable supplementary data.
