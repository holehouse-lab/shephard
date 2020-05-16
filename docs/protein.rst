Protein
=================
Each individual sequence entry is represented as its own Protein. Proteins contain four types of metadata

1. **Attributes**: Arbitrary key-value pairings where the key must be a string and the value can be any Python object (:code:`str`, :code:`int`, :code:`dict`, figure handle, simulation object)
2. **Tracks**: Vectors that are the same length of the protein, are either numerical or symbolic, and contain residue-specific metadata projected over the whole sequence
3. **Domains**: Continuous sub-regions within the protein
4. **Sites**: Individual site positions within the protein

Each of these can be accessed using associated functions.

**IMPORTANT NOTE**: In the field of biology and bioinformatics the indexing system used to describe regions and sites (1) starts at 1 and (2) is inclusive. For example, if we had a sequence of `MAEPQRDG` and wanted the region defined by 2-4 this would reflect `AEP`. In contrast the Python programming language (1) indexes from 0 and (2) is exclusive for ranges, so (using slice notation, for example) `MAEPQRDG[2:4]` would yield EP. Within SHEPHARD, all user-facing tools operate using biology-style indexing. This means you can read in data directly from native databases without worrying about conversion, and use the same number scheme always. Because of this, to sub-select regions of the protein sequence we STRONGLY recommend using the functions :code:`.get_sequence_region(start,end)` or :code:`.get_sequence_context(pos)`, rather than directly slicing the :code:`.sequence` property. 


.. autoclass:: shephard.protein.Protein

Protein properties
...........

.. autofunction:: shephard.proteome.Protein.name
.. autofunction:: shephard.proteome.Protein.proteome
.. autofunction:: shephard.proteome.Protein.sequence
.. autofunction:: shephard.proteome.Protein.unique_ID
.. autofunction:: shephard.proteome.Protein.attributes

.. autofunction:: shephard.proteome.Protein.get_sequence_region
.. autofunction:: shephard.proteome.Protein.get_sequence_context
.. autofunction:: shephard.proteome.Protein.check_sequence_is_valid
.. autofunction:: shephard.proteome.Protein.convert_to_valid


Attribute functions
.......................
.. autofunction:: shephard.proteome.Protein.attributes
.. autofunction:: shephard.proteome.Protein.attribute
.. autofunction:: shephard.proteome.Protein.add_attribute

Domain functions
..................

.. autofunction:: shephard.proteome.Protein.domains
.. autofunction:: shephard.proteome.Protein.domain
.. autofunction:: shephard.proteome.Protein.domain_types
.. autofunction:: shephard.proteome.Protein.get_domains_by_type
.. autofunction:: shephard.proteome.Protein.add_domain
.. autofunction:: shephard.proteome.Protein.add_domains
.. autofunction:: shephard.proteome.Protein.build_domains


Track functions
..................
.. autofunction:: shephard.proteome.Protein.tracks
.. autofunction:: shephard.proteome.Protein.track
.. autofunction:: shephard.proteome.Protein.get_track_values
.. autofunction:: shephard.proteome.Protein.get_track_symbols
.. autofunction:: shephard.proteome.Protein.add_track
.. autofunction:: shephard.proteome.Protein.build_track
.. autofunction:: shephard.proteome.Protein.buld_track_values_from_sequence
.. autofunction:: shephard.proteome.Protein.buld_track_symbols_from_sequence



Site functions
..................
.. autofunction:: shephard.proteome.Protein.sites
.. autofunction:: shephard.proteome.Protein.site
.. autofunction:: shephard.proteome.Protein.add_site
.. autofunction:: shephard.proteome.Protein.add_sites_by_position
.. autofunction:: shephard.proteome.Protein.get_sites_by_range
.. autofunction:: shephard.proteome.Protein.get_sites_by_type
.. autofunction:: shephard.proteome.Protein.get_sites_by_type_and_range



