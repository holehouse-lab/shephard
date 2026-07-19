Track
=================

Tracks represent annotations associated with the full length of the protein (i.e. per-residue sequence annotations). Tracks are added to proteins using the ``Protein.add_track()`` function, built on-the-fly from the sequence with ``Protein.build_track_values_from_sequence()`` / ``Protein.build_track_symbols_from_sequence()``, or read in using functions in the ``shephard.interfaces.si_tracks`` module.

A Track must have a ``track name``, and is **either** a *values* track (a per-residue vector of numbers, cast to ``float``) **or** a *symbols* track (a per-residue list/string of symbols) — never both. A Track has a strict one-to-one mapping with the protein sequence: SHEPHARD validates at construction time that the number of values/symbols exactly matches the protein length, and that a given track name is consistently a values track or a symbols track across the whole Proteome.

Tracks are appropriate when an analysis produces a *continuous* output that can be projected along the sequence (e.g. a disorder score, a hydropathy window, an AlphaFold pLDDT). If, instead, your analysis reveals discrete regions or individual positions, ``Domain`` or ``Site`` annotations are more appropriate, respectively.

Tracks for a given protein can be requested using the ``protein.track(<track name>)`` function (returns the ``Track`` object), or the underlying data accessed directly with ``protein.get_track_values()`` / ``protein.get_track_symbols()``. Sub-regions can be extracted with ``values_region()`` / ``symbols_region()``, and individual positions with ``value()`` / ``symbol()``. All tracks of a given name across an entire proteome can be requested with ``Proteome.get_tracks_by_name()``.

Tracks can be removed from proteins using the ``Protein.remove_track()`` function.

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
