apis
=================

the :code:`apis` package contains a collection of tools for working with non-SHEPHARD format files. In particular, SHEPHARD uses FASTA files for sequence storage, which is not a format SHEPHARD has control over, so interaction with FASAT files occurs via APIs.

We anticipate the number of :code:`api` modules to remain small, but for some tools that the lab has control over, direct interaction via an API makes sense. In general, we use interfaces for interacting with data to ensure weak coupling between SHEPHARD and other tools.

fasta
------------
A non-uniprot FASTA file can also be read in using the :code:`fasta` module

.. autofunction:: shephard.apis.fasta.fasta_to_proteome

uniprot
------------
The uniprot module provides tools for working with uniprot data. Right now, only an automatic uniprot FASTA file parser is in place, but over time we plan to add more generic file I/O for uniprot derived files, given the robustness and broad usership.

.. autofunction:: shephard.apis.uniprot.uniprot_fasta_to_proteome
.. autofunction:: shephard.apis.uniprot.uniprot_proteome_to_fasta







