Domain
=================

Domains represent annotations associated with subregions along the sequence. Domains are added to proteins using the ``Protein.add_domain()`` function, or using functions in the ``shephard.interfaces.si_domains`` module.

Domains must have a ``domain type``, and a ``domain type`` as well as a start and end position. Each domain generally should have a unique domain name, although this is not strictly enforced. In contrast, many domains can and will be of common domain types. Domains also know the position in the sequence they come from, the underlying residue, and can extract Site and Track information associated with the domain.

Domains for a given protein can be called using the ``protein.domain(domain_name)`` function. However, in general it's more useful to either request all domains using ``protein.domains`` (which returns a list of all domains in protein) or to request specific domains based on their position, location, type, or some combinatin of the two. Explicit functions for these types of requests are included in the ``Protein`` object. Finally, all domains (or all domains of a specific type) can be requested from an entire proteome using ``Proteome`` object functions.

Domains can be removed from proteins using ``Protein.remove_domain()`` function.

.. autoclass:: shephard.domain.Domain


Domain Properties
...................
.. autofunction:: shephard.domain.Domain.start
.. autofunction:: shephard.domain.Domain.end
.. autofunction:: shephard.domain.Domain.protein
.. autofunction:: shephard.domain.Domain.sequence
.. autofunction:: shephard.domain.Domain.domain_type

Domain Functions
.................
.. autofunction:: shephard.domain.Domain.inside_domain


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
.. autofunction:: shephard.domain.Domain.get_sites_by_type


Domain Track Functions
...........................
.. autofunction:: shephard.domain.Domain.get_track_values
.. autofunction:: shephard.domain.Domain.get_track_symbols
