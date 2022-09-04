Site
=================
Sites represent annotations associated with single positions along the sequence. Sites are added to proteins using the ``Protein.add_site()`` function, or using functions in the ``shephard.interfaces.si_sites`` module.

Sites must have a ``site type``, and can be associated with a ``symbol`` or a ``value``. Sites also know the position in the sequence they come from, the underlying residue, and can extract track information associated with the site.

Sites for a given protein can be called using the ``protein.site(<site_position>)`` function. However, in general it's more useful to either request all sites using ``protein.sites`` (which returns a list of all sites in protein) or to request specific sites based on their position, location, type, or some combinatin of the two. Explicit functions for these types of requests are included in the ``Protein`` object. Finally, all sites (or all sites of a specific type) can be requested from an entire proteome using ``Proteome`` object functions.

Sites can be removed from proteins using ``Protein.remove_site()`` function.


.. autoclass:: shephard.site.Site


Site properties
................
.. autofunction:: shephard.site.Site.residue
.. autofunction:: shephard.site.Site.position
.. autofunction:: shephard.site.Site.protein
.. autofunction:: shephard.site.Site.site_type
.. autofunction:: shephard.site.Site.symbol
.. autofunction:: shephard.site.Site.value


Site sequence functions
............................
.. autofunction:: shephard.site.Site.get_local_sequence_context


Site Attribute Functions
...........................
.. autofunction:: shephard.site.Site.attributes
.. autofunction:: shephard.site.Site.attribute
.. autofunction:: shephard.site.Site.add_attribute
.. autofunction:: shephard.site.Site.remove_attribute


Site Track Functions
.......................
.. autofunction:: shephard.site.Site.get_track_values
.. autofunction:: shephard.site.Site.get_track_symbols

