Overview
===============

**Version 0.1.0**

SHEPHARD (Sequence-based Hierachical and Extendable Platform for High-throughput Analysis of Region of Disorder) is a modular framework for building programatic analyses of proteins. It does not itself perform any kind of analysis, but instead provides all of the nitty-gritty tools that are necessary but complicated for doing largescale analysis of protein datasets. By hidding all of the logistical complexity, a user can quickly ask complex questions.

SHEPHARD stores information in a bidirectional hierachical format. A collection of proteins is called a Proteome. A Proteome contains 1 or more Proteins. Each Protein contains one or more Domains, Sites and Tracks. A Domain is a contigous subregion in the protein. A site is a residue-specific annotation, where multiple sites can coexist at the same position. A Track is a numerical or symbolic vector that maps some arbitrary value to each residue in the protein sequence.

SHEPHARD makes it easy to define your own domains, tracks, and sites, and read this into and out of files. In this way, customizable and reproducible analysis becomes extremley easy, and the same code can be used to perform analagous analyses for very different questions. This means you can spend more time analyzing the results and thinking about the questions, and less time on the software development side.

--------------
Installation
--------------

Installation for now should be performed using the release candidate you were sent. This will look something like :code:`shephard.....tar.gz.`. To install we strongly recommend having a conda environment set up (the scope of conda setup is beyond that of this documentation) and from within that conda environment use :code:`pip` ::

    pip install shephard...tar.gz

This should install without issue, and once installed, shephard is available for import in any Python code you write when executed from within that conda environment.
