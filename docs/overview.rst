Overview
===============

SHEPHARD (**S**\equence-based **H**\ierarchical and **E**\xtendable **P**\latform for **H**\igh-throughput **A**\nalysis of **R**\egion of **D**\isorder) is a modular framework for performing programmatic analyses of proteins. It does not itself perform any kind of analysis, but instead provides all of the nitty-gritty tools that are necessary but complicated for doing large-scale analysis of protein datasets. By hiding all of the logistical complexity, a user can quickly ask complex questions. The goal here is to make it trivial to ask what might otherwise be computationally complex things to address.



--------------
Installation
--------------

Installation for now should be performed using the current RC (release candidate). This will be file that looks something like :code:`shephard.....tar.gz.`. To install we strongly recommend having a :code:`conda` environment set up (the scope of :code:`conda` setup is beyond that of this documentation) and from within that :code:`conda` environment use :code:`pip` ::

    pip install shephard@git+git://github.com/holehouse-lab/shephard.git


This should install without issue, and once installed, SHEPHARD is available for import in any Python code you write when executed from within that :code:`conda` environment.

If installing a package via :code:`pip` is new to you, we recommend spending a moment getting used to the Python programming language and the concept of :code:`conda` environments.



--------------------
Demos and examples
--------------------
Once installed, SHEPHARD makes it very easy to work with large protein datasets. As an example, we have a collection of basic notebooks showing functionality located at on the `SHEPHARD GitHub supporting data page <https://github.com/holehouse-lab/supportingdata/tree/master/2022/ginell_2022/example_notebooks>`_. 




-------------------------
SHEPHARD architecture
-------------------------
SHEPHARD stores information in a bidirectional hierarchical format. This hierarchy is built from a few key building blocks: **Proteomes**, **Proteins**, **Domains**, **Sites** and **Tracks**.

A collection of proteins is called a Proteome. A Proteome contains 1 or more Proteins. Each Protein contains zero or more Domains, Sites and Tracks. A Domain is a contiguous subregion in the protein. A site is a residue-specific annotation, where multiple sites can coexist at the same position. A Track is a numerical or symbolic vector that maps some arbitrary value to each residue in the protein sequence.

SHEPHARD makes it easy to define your own domains, tracks, and sites, and read this into and out of files. It does this by defining specific file formats for these types of files, where those file formats are described in the interfaces documentation page. Sequences are read in via FASTA format, and SHEPHARD piggy backs on the robust FASTA parser `protfasta <https://protfasta.readthedocs.io/>`_  that was developed by the Holehouse lab and released in early 2020.

Because SHEPHARD interacts with the 'outside world' (other types of data) via **interfaces** that require a specific file format (or programmatic input format), it maintains a loose coupling between the internal architecture and any other tool, pipeline, or dataset. The interfaces ensure a contract between SHEPHARD and the outside world, which means the user can write robust code using SHEPHARD knowing that they are insulated from changes to those external tools. This removes direct dependencies on other tools while allowing data from those tools to be readily accessible.

In this way, customizable and reproducible analysis becomes extremely easy, and the same code can be used to perform analogous analyses for very different questions. This means you can spend more time analyzing the results and thinking about the questions, and less time on the software development side.


---------------------------------
General workflow in SHEPHARD
---------------------------------

The general work flow we envisage for people using SHEPHARD is as follows:

1. Read in a set of sequences from a FASTA file and generate a new Proteome using one of the :code:`api` functions (described elsewhere in the documentation)

2. Load in some annotations (domains, sites, tracks or protein_attributes) via the :code:`interface` functions.

3. Perform analysis in a way that queries the relationship between sequence and those different annotations.





