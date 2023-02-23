Protein
=================


Protein overview
----------------------

Protein objects are the major unit of data storage in SHEPHARD, and each `Proteomes <https://shephard.readthedocs.io/en/latest/proteome.html>`_ is made up of zero or more Protein object. Each Protein object is associated with one amino acid sequence, and a collection of various possible annotations outlined below.

Proteins contain four possible types of metadata:

1. **Attributes**: Arbitrary key-value pairings, where the key must be a string and the value can be any Python object (``str``, ``int``, ``dict``, figure handle, simulation object, lambda function, or any other variable). While attributes with complex datatypes cannot be easily saved within SHEPHARD, simple datatypes (i.e. those that can be cast to strings) can be written and read using the `si_protein_attributes <https://shephard.readthedocs.io/en/latest/interfaces.html#si-protein-attributes>`_ module. 
2. **Tracks**: Vectors that are the same length of the protein, are either numerical or symbolic, and contain residue-specific metadata projected over the whole sequence
3. **Domains**: Continuous sub-regions within the protein
4. **Sites**: Individual site positions within the protein

Each of these can be accessed using associated functions.

Protein indexing
----------------------

In the field of biology, the indexing system used to describe regions and sites **(1)** starts at 1 and **(2)** is inclusive. 

For example, if we had a sequence of `MAEPQRDG` and wanted the region defined by 2-4 this would reflect ``AEP``. In contrast, the Python programming language **(1)** indexes from 0 and **(2)** is exclusive for ranges, so (using slice notation, for example) ``MAEPQRDG[2:4]`` would yield ``EP``. 

Within SHEPHARD, all user-facing tools operate using biology-style indexing. This means you can read in data directly from native databases without worrying about conversion, and use the same number scheme always. Because of this, to sub-select regions of the protein sequence we STRONGLY recommend using the functions :code:`.get_sequence_region(start,end)` or :code:`.get_sequence_context(pos)`, rather than directly slicing the :code:`.sequence` property. 

As an example:

.. code-block:: python

    from shephard import Proteome
    P = Proteome()

    # add a protein with the sequence AACCDDEEFF, the name 'my coold protein'
    # and the unique_ID 'test1'
    P.add_protein('AACCDDEEFF', 'my cool protein', 'test1')

    # slice notation into sequence (BAD) 
    print(P.protein('test1').sequence[1:3])
    >>> AC

    # using get_sequence_region (GOOD)
    print(P.protein('test1').get_sequence_region(1,3))
    >>> AAC

    # using get_residue to excise a specific residue; the second
    # residue should be an A
    print(P.protein('test1').residue(2))
    >>> A

    





.. autoclass:: shephard.protein.Protein

Protein properties
---------------------

.. autofunction:: shephard.protein.Protein.name
.. autofunction:: shephard.protein.Protein.proteome
.. autofunction:: shephard.protein.Protein.sequence
.. autofunction:: shephard.protein.Protein.unique_ID


Sequence functions
---------------------

.. autofunction:: shephard.protein.Protein.residue
.. autofunction:: shephard.protein.Protein.get_sequence_region
.. autofunction:: shephard.protein.Protein.get_sequence_context
.. autofunction:: shephard.protein.Protein.check_sequence_is_valid
.. autofunction:: shephard.protein.Protein.convert_to_valid


Attribute functions
---------------------

.. autofunction:: shephard.protein.Protein.attributes
.. autofunction:: shephard.protein.Protein.attribute
.. autofunction:: shephard.protein.Protein.add_attribute
.. autofunction:: shephard.protein.Protein.remove_attribute

Domain functions
---------------------

.. autofunction:: shephard.protein.Protein.domains
.. autofunction:: shephard.protein.Protein.domain_names
.. autofunction:: shephard.protein.Protein.domain
.. autofunction:: shephard.protein.Protein.domain_types
.. autofunction:: shephard.protein.Protein.add_domain
.. autofunction:: shephard.protein.Protein.add_domains
.. autofunction:: shephard.protein.Protein.build_domains
.. autofunction:: shephard.protein.Protein.get_domains_by_type
.. autofunction:: shephard.protein.Protein.get_domains_by_position
.. autofunction:: shephard.protein.Protein.get_domains_by_position_and_type
.. autofunction:: shephard.protein.Protein.get_domains_by_range


Track functions
---------------------

.. autofunction:: shephard.protein.Protein.tracks
.. autofunction:: shephard.protein.Protein.track
.. autofunction:: shephard.protein.Protein.track_names
.. autofunction:: shephard.protein.Protein.get_track_values
.. autofunction:: shephard.protein.Protein.get_track_symbols
.. autofunction:: shephard.protein.Protein.add_track
.. autofunction:: shephard.protein.Protein.remove_track
.. autofunction:: shephard.protein.Protein.build_track
.. autofunction:: shephard.protein.Protein.build_track_values_from_sequence
.. autofunction:: shephard.protein.Protein.build_track_symbols_from_sequence



Site functions
---------------------

.. autofunction:: shephard.protein.Protein.sites
.. autofunction:: shephard.protein.Protein.site
.. autofunction:: shephard.protein.Protein.site_types
.. autofunction:: shephard.protein.Protein.site_positions
.. autofunction:: shephard.protein.Protein.add_site
.. autofunction:: shephard.protein.Protein.remove_site
.. autofunction:: shephard.protein.Protein.get_sites_by_position
.. autofunction:: shephard.protein.Protein.get_sites_by_range
.. autofunction:: shephard.protein.Protein.get_sites_by_type
.. autofunction:: shephard.protein.Protein.get_sites_by_type_and_range



