Proteome
=================

The ``Proteome`` is the base container in SHEPHARD. A Proteome holds one or more ``Protein`` objects (keyed by a unique ID), and can itself carry arbitrary key-value attributes (useful for dataset-level metadata such as the source database, build date, or organism).

Proteomes are usually constructed indirectly by reading a FASTA file with one of the :doc:`apis` functions (e.g. ``uniprot.uniprot_fasta_to_proteome`` or ``fasta.fasta_to_proteome``), but can also be built directly from a list of protein dictionaries or a list of existing ``Protein`` objects (in which case the proteins are deep-copied so the new Proteome is fully independent).

A Proteome also behaves like a standard Python container: ``len(P)`` gives the number of proteins, ``for protein in P`` iterates over ``Protein`` objects directly, ``'some_unique_ID' in P`` (or ``protein_object in P``) tests membership, and ``P[0]``, ``P[-1]``, or ``P[10:20]`` index or slice into it. Indexing follows the usual Python conventions, so negative indices count back from the last protein and an out-of-range index raises an ``IndexError``.

In addition to per-protein access, the Proteome provides convenient *aggregate* accessors that operate across every protein at once — ``domains``, ``sites``, and ``tracks`` return every annotation of that type in the whole Proteome, while ``get_domains_by_type``, ``get_sites_by_type``, and ``get_tracks_by_name`` return the subset matching a given identifier. Because the hierarchy is bidirectional, each returned annotation can be mapped back to its parent Protein via its ``.protein`` accessor.

.. autoclass:: shephard.proteome.Proteome


Proteome functions
....................

.. autofunction:: shephard.proteome.Proteome.check_unique_ID
.. autoproperty:: shephard.proteome.Proteome.attributes
.. autofunction:: shephard.proteome.Proteome.attribute
.. autofunction:: shephard.proteome.Proteome.add_attribute
.. autofunction:: shephard.proteome.Proteome.remove_attribute
.. autofunction:: shephard.proteome.Proteome.__iter__
.. autofunction:: shephard.proteome.Proteome.__contains__
.. autofunction:: shephard.proteome.Proteome.__getitem__
.. autofunction:: shephard.proteome.Proteome.__len__



Protein functions
....................

.. autoproperty:: shephard.proteome.Proteome.proteins
.. autofunction:: shephard.proteome.Proteome.protein
.. autofunction:: shephard.proteome.Proteome.add_protein
.. autofunction:: shephard.proteome.Proteome.add_proteins
.. autofunction:: shephard.proteome.Proteome.remove_protein
.. autofunction:: shephard.proteome.Proteome.remove_proteins


Domain properties
....................

.. autoproperty:: shephard.proteome.Proteome.domains
.. autoproperty:: shephard.proteome.Proteome.unique_domain_types
.. autofunction:: shephard.proteome.Proteome.get_domains_by_type


Site properties
....................

.. autoproperty:: shephard.proteome.Proteome.sites
.. autoproperty:: shephard.proteome.Proteome.unique_site_types
.. autofunction:: shephard.proteome.Proteome.get_sites_by_type


Track properties
....................

.. autoproperty:: shephard.proteome.Proteome.tracks
.. autofunction:: shephard.proteome.Proteome.get_tracks_by_name
.. autoproperty:: shephard.proteome.Proteome.unique_track_names
.. autoproperty:: shephard.proteome.Proteome.track_names_to_track_type
