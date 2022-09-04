Track
=================

Tracks represent annotations associated with the full length of the protein (i.e. per-residue sequence annotations). Tracks are added to proteins using the ``Protein.add_track()`` function, or using functions in the ``shephard.interfaces.si_tracks`` module.

Tracks must have a ``track name``, as well as either a per-residue symbol annotation or per-residue value annotation, meaning Tracks can be either symbol tracks of value tracks. 

Tracks for a given protein can be called using the ``protein.track(<track name>)`` function.

Tracks can be removed from proteins using ``Protein.remove_track()`` function.

.. autoclass:: shephard.track.Track


Track Properties
................
.. autofunction:: shephard.track.Track.name
.. autofunction:: shephard.track.Track.values
.. autofunction:: shephard.track.Track.symbols
.. autofunction:: shephard.track.Track.protein
.. autofunction:: shephard.track.Track.track_type

Track Functions
...............
.. autofunction:: shephard.track.Track.symbol
.. autofunction:: shephard.track.Track.value
.. autofunction:: shephard.track.Track.symbols_region
.. autofunction:: shephard.track.Track.values_region


Track Attribute Functions
...........................
.. autofunction:: shephard.track.Track.attributes
.. autofunction:: shephard.track.Track.attribute
.. autofunction:: shephard.track.Track.add_attribute
.. autofunction:: shephard.track.Track.remove_attribute


