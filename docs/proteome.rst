Proteome
=================

The ``Proteome`` is the base container in SHEPHARD. A Proteome holds one or more ``Protein`` objects (keyed by a unique ID), and can itself carry arbitrary key-value attributes (useful for dataset-level metadata such as the source database, build date, or organism).

Proteomes are usually constructed indirectly by reading a FASTA file with one of the :doc:`apis` functions (e.g. ``uniprot.uniprot_fasta_to_proteome`` or ``fasta.fasta_to_proteome``), but can also be built directly from a list of protein dictionaries or a list of existing ``Protein`` objects (in which case the proteins are deep-copied so the new Proteome is fully independent).

In addition to per-protein access, the Proteome provides convenient *aggregate* accessors that operate across every protein at once — ``domains``, ``sites``, and ``tracks`` return every annotation of that type in the whole Proteome, while ``get_domains_by_type``, ``get_sites_by_type``, and ``get_tracks_by_name`` return the subset matching a given identifier. Because the hierarchy is bidirectional, each returned annotation can be mapped back to its parent Protein via its ``.protein`` accessor.

.. autoclass:: shephard.proteome.Proteome


Proteome functions
....................

.. autofunction:: shephard.proteome.Proteome.check_unique_ID
.. autofunction:: shephard.proteome.Proteome.attributes
.. autofunction:: shephard.proteome.Proteome.attribute
.. autofunction:: shephard.proteome.Proteome.add_attribute
.. autofunction:: shephard.proteome.Proteome.remove_attribute
.. autofunction:: shephard.proteome.Proteome.__iter__
.. autofunction:: shephard.proteome.Proteome.__contains__
.. autofunction:: shephard.proteome.Proteome.__getitem__
.. autofunction:: shephard.proteome.Proteome.__len__



Protein functions
....................

.. autofunction:: shephard.proteome.Proteome.proteins
.. autofunction:: shephard.proteome.Proteome.protein
.. autofunction:: shephard.proteome.Proteome.add_protein
.. autofunction:: shephard.proteome.Proteome.add_proteins
.. autofunction:: shephard.proteome.Proteome.remove_protein
.. autofunction:: shephard.proteome.Proteome.remove_proteins


Domain properties
....................

.. autofunction:: shephard.proteome.Proteome.domains
.. autofunction:: shephard.proteome.Proteome.unique_domain_types
.. autofunction:: shephard.proteome.Proteome.get_domains_by_type


Site properties
....................

.. autofunction:: shephard.proteome.Proteome.sites
.. autofunction:: shephard.proteome.Proteome.unique_site_types
.. autofunction:: shephard.proteome.Proteome.get_sites_by_type


Track properties
....................

.. autofunction:: shephard.proteome.Proteome.tracks
.. autofunction:: shephard.proteome.Proteome.get_tracks_by_name
.. autofunction:: shephard.proteome.Proteome.unique_track_names
.. autofunction:: shephard.proteome.Proteome.track_names_to_track_type
