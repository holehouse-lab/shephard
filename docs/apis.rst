apis
=================

the :code:`apis` package contains a collection of tools for working with non-SHEPHARD format files. In particular, SHEPHARD uses FASTA files for sequence storage, which is not a format SHEPHARD has control over, so interaction with FASTA files occurs via APIs.

We anticipate the number of :code:`api` modules to remain small, but for some tools that the lab has control over, direct interaction via an API makes sense. In general, we use interfaces for interacting with data to ensure weak coupling between SHEPHARD and other tools.

fasta
------------
A non-uniprot FASTA file can also be read in using the :code:`fasta` module

.. autofunction:: shephard.apis.fasta.fasta_to_proteome
.. autofunction:: shephard.apis.fasta.proteome_to_fasta

uniprot
------------
The uniprot module provides tools for working with uniprot data. Right now, only an automatic uniprot FASTA file parser is in place, but over time we plan to add more generic file I/O for uniprot derived files, given the robustness and broad usership.

.. autofunction:: shephard.apis.uniprot.uniprot_fasta_to_proteome
.. autofunction:: shephard.apis.uniprot.uniprot_proteome_to_fasta


metapredict_api
-----------------
The metapredict module provides tools for annotating proteome-scale information using metapredict. This depends on having metapredict V2-FF available, which is not defined as a hard requirement, but, if used, enables the annotation of entire proteomes in terms of IDRs and disorder scores with a single function.

.. autofunction:: shephard.apis.metapredict_api.annotate_proteome_with_disorder_track
.. autofunction:: shephard.apis.metapredict_api.annotate_proteome_with_disordered_domains
.. autofunction:: shephard.apis.metapredict_api.annotate_proteome_with_disorder_tracks_and_disordered_domains

albatross_api
-----------------
The ALBATROSS API interfaces with the radius of gyration (Rg) and end-to-end distance (Re) predictions provided by ALBATROSS (Lotthammer et al. Nat. Meth. 2024). This depends on `sparrow <https://github.com/idptools/sparrow>`_ being installed, but enables you to annotate at either a proteine level or a domain level sequences with predicted Rg and Re values. 

.. autofunction:: shephard.apis.albatross_api.annotate_proteome_with_dimensions
.. autofunction:: shephard.apis.albatross_api.annotate_domains_with_dimensions




