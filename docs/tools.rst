Tools
=================
SHEPHARD ships with a collection of stateless helper functions ("tools") for analyzing, constructing, and manipulating annotations. These live under the ``shephard.tools`` sub-package, organized by the object they primarily operate on. Tools never silently mutate the objects passed to them unless explicitly documented to do so; most return new data structures (lists, dictionaries, vectors) that you can then add to a Proteome via the relevant ``interfaces`` function.


Attribute tools
-------------------------

Attribute-associated functions are located in the ``shephard.tools.attribute_tools`` module. These tools help with bulk manipulation of the arbitrary key-value attributes attached to SHEPHARD objects (for example, casting attribute values that were read from a text file as strings into numeric types).

.. autofunction:: shephard.tools.attribute_tools.cast_attributes


Domain tools
-------------------------

Domain-associated functions are located in the ``shephard.tools.domain_tools`` module. These tools enable analysis, construction, and manipulation of Domain objects — for example testing whether two domains overlap (and by how much), filling in the "gaps" between annotated domains, or discretizing a continuous values Track into discrete domains.

.. autofunction:: shephard.tools.domain_tools.domain_overlap
.. autofunction:: shephard.tools.domain_tools.domain_overlap_fraction
.. autofunction:: shephard.tools.domain_tools.domain_overlap_by_position
.. autofunction:: shephard.tools.domain_tools.build_missing_domains
.. autofunction:: shephard.tools.domain_tools.build_domains_from_track_values


Sequence tools
--------------------------

Sequence-associated functions are located in the ``shephard.tools.sequence_tools`` module. These tools enable manipulation or search of sequence information, such as building a single concatenated "mega-string" of many sequences for compositional statistics, or locating (regex-enabled) motif occurrences with biology-style indexing.

.. autofunction:: shephard.tools.sequence_tools.build_mega_string
.. autofunction:: shephard.tools.sequence_tools.find_string_positions


Site tools
--------------------------

Site-associated functions are located in the ``shephard.tools.site_tools`` module. These tools enable manipulation or search of site information, such as computing a sliding-window site-density vector that can be added back as a Track.

.. autofunction:: shephard.tools.site_tools.build_site_density_vector


Track tools
--------------------------

Track-associated functions are located in the ``shephard.tools.track_tools`` module. These tools enable construction and manipulation of Track data — for example binarizing a continuous vector against a threshold, or projecting the Domain architecture of every protein in a Proteome onto a per-residue symbol track.

.. autofunction:: shephard.tools.track_tools.binerize
.. autofunction:: shephard.tools.track_tools.build_track_from_domains


Experimental tools
--------------------------

The ``shephard.tools.experimental`` module contains newer, less-stable analysis helpers. The API of these functions may change in future releases, but they are documented here for completeness. They were used, for example, to quantify the enrichment of post-translational modification sites within disordered regions relative to a within-protein null model.

.. autofunction:: shephard.tools.experimental.get_site_density_in_domain_normalized_by_protein
.. autofunction:: shephard.tools.experimental.get_site_density_percentile_normalized_by_protein
