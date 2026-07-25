Site
=================
Sites represent annotations associated with single positions along the sequence. Sites are added to proteins using the ``Protein.add_site()`` function, or using functions in the ``shephard.interfaces.si_sites`` module.

Sites must have a ``site type``, and can be associated with a ``symbol`` (string) and/or a ``value`` (numerical, cast to ``float``). Multiple Sites can coexist at the same residue position, which makes Sites well-suited to representing things like post-translational modifications, mutations, or binding sites. Sites also know the position in the sequence they come from and the underlying residue, and can extract Domain and Track information that overlaps the site.

Sites for a given protein can be requested using the ``protein.site(<site_position>)`` function. However, it is generally more useful to either request all sites using ``protein.sites`` (which returns a list of all sites in the protein) or to request specific sites based on their position, range, type, or some combination of these. Explicit functions for these requests are included in the ``Protein`` object. Finally, all sites (or all sites of a specific type) can be requested from an entire proteome using ``Proteome`` object functions.

A Site's ``value`` and ``symbol`` can be updated in place using ``update_site_value()`` and ``update_site_symbol()`` (each accepts ``None`` to clear the field). Sites can be removed from proteins using the ``Protein.remove_site()`` function.


.. autoclass:: shephard.site.Site


Site properties
................
.. autoproperty:: shephard.site.Site.residue
.. autoproperty:: shephard.site.Site.position
.. autoproperty:: shephard.site.Site.protein
.. autoproperty:: shephard.site.Site.site_type
.. autoproperty:: shephard.site.Site.symbol
.. autoproperty:: shephard.site.Site.value


Site update functions
............................
.. autofunction:: shephard.site.Site.update_site_value
.. autofunction:: shephard.site.Site.update_site_symbol


Site sequence functions
............................
.. autofunction:: shephard.site.Site.get_local_sequence_context


Site Attribute Functions
...........................
.. autoproperty:: shephard.site.Site.attributes
.. autofunction:: shephard.site.Site.attribute
.. autofunction:: shephard.site.Site.add_attribute
.. autofunction:: shephard.site.Site.remove_attribute


Site Domain Functions
.......................
.. autofunction:: shephard.site.Site.get_domains


Site Track Functions
.......................
.. autofunction:: shephard.site.Site.get_track_values
.. autofunction:: shephard.site.Site.get_track_value
.. autofunction:: shephard.site.Site.get_track_symbols
.. autofunction:: shephard.site.Site.get_track_symbol
