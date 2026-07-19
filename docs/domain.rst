Domain
=================

Domains represent annotations associated with contiguous subregions along the sequence. Domains are added to proteins using the ``Protein.add_domain()`` function, or using functions in the ``shephard.interfaces.si_domains`` module.

Domains must have a ``domain_type`` as well as a start and end position (using biology-style, 1-based, inclusive indexing). Each domain has a ``domain_name`` that is auto-generated as ``<domain_type>_<start>_<end>`` and is unique within a protein; the ``autoname`` option allows perfectly-overlapping domains of the same type to coexist by appending an incrementer. Many domains can (and typically will) share a common ``domain_type``, which is what you select on. Domains also know the position in the sequence they come from, the underlying residue sequence, and can extract Site and Track information associated with the domain.

Domains for a given protein can be requested using the ``protein.domain(domain_name)`` function. However, it is generally more useful to either request all domains using ``protein.domains`` (which returns a list of all domains in the protein, sorted by start position) or to request specific domains by position, range, type, or some combination of these. Explicit functions for these requests are included in the ``Protein`` object. Finally, all domains (or all domains of a specific type) can be requested from an entire proteome using ``Proteome`` object functions.

Domains can be removed from proteins using the ``Protein.remove_domain()`` function.

.. autoclass:: shephard.domain.Domain


Domain Properties
...................
.. autofunction:: shephard.domain.Domain.start
.. autofunction:: shephard.domain.Domain.end
.. autofunction:: shephard.domain.Domain.protein
.. autofunction:: shephard.domain.Domain.sequence
.. autofunction:: shephard.domain.Domain.domain_type
.. autofunction:: shephard.domain.Domain.domain_name


Domain Functions
.................
.. autofunction:: shephard.domain.Domain.inside_domain
.. autofunction:: shephard.domain.Domain.domain_overlap
.. autofunction:: shephard.domain.Domain.update_domain_name


Domain Attribute Functions
...........................
.. autofunction:: shephard.domain.Domain.attributes
.. autofunction:: shephard.domain.Domain.attribute
.. autofunction:: shephard.domain.Domain.add_attribute
.. autofunction:: shephard.domain.Domain.remove_attribute


Domain Site Functions
...........................
.. autofunction:: shephard.domain.Domain.sites
.. autofunction:: shephard.domain.Domain.site
.. autofunction:: shephard.domain.Domain.site_positions
.. autofunction:: shephard.domain.Domain.get_sites_by_type


Domain Track Functions
...........................
.. autofunction:: shephard.domain.Domain.get_track_values
.. autofunction:: shephard.domain.Domain.get_track_symbols
