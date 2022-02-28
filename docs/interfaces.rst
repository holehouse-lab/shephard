interfaces
=================

The interfaces package is a sub-package within SHEPHARD that deals with reading in **SHEPHARD** formatted input files. These are files with a specific format that were developed for SHEPHARD. The interfaces mean that, as long as you can write your data in a format that complies with a track, site, domain, or protein_attribute file, you can be sure it will be correctly read into SHEPHARD and is then accessible within the larger framework. 

In all cases, each line in an input file corresponds to a single bit of information that maps to a specific protein, but multiple lines can map to the same protein. 


si_sites
------------
Functions associated with the ``si_sites`` module enable the reading and writing of SHEPHARD Sites files.

.. autofunction:: shephard.interfaces.si_sites.add_sites_from_file
.. autofunction:: shephard.interfaces.si_sites.add_sites_from_dictionary
.. autofunction:: shephard.interfaces.si_sites.write_sites


si_domains
------------
Functions associated with the ``si_domains`` module enable the reading and writing of SHEPHARD Domains files.


.. autofunction:: shephard.interfaces.si_domains.add_domains_from_file
.. autofunction:: shephard.interfaces.si_domains.add_domains_from_dictionary
.. autofunction:: shephard.interfaces.si_domains.add_domain_attributes_from_file
.. autofunction:: shephard.interfaces.si_domains.add_domain_attributes_from_dictionary
.. autofunction:: shephard.interfaces.si_domains.write_domains
.. autofunction:: shephard.interfaces.si_domains.write_domains_from_list


si_tracks
------------
Functions associated with the ``si_tracks`` module enable the reading and writing of SHEPHARD Tracks files.

.. autofunction:: shephard.interfaces.si_tracks.add_tracks_from_dictionary
.. autofunction:: shephard.interfaces.si_tracks.add_tracks_from_file
.. autofunction:: shephard.interfaces.si_tracks.write_all_tracks_separate_files
.. autofunction:: shephard.interfaces.si_tracks.write_all_values_tracks_single_file
.. autofunction:: shephard.interfaces.si_tracks.write_all_symbols_tracks_single_file
.. autofunction:: shephard.interfaces.si_tracks.write_track


si_protein_attributes
----------------------
Functions associated with the ``si_protein_attributes`` module enable the reading and writing of SHEPHARD Protein attribute files.

.. autofunction:: shephard.interfaces.si_protein_attributes.add_protein_attributes_from_dictionary
.. autofunction:: shephard.interfaces.si_protein_attributes.add_protein_attributes_from_file
.. autofunction:: shephard.interfaces.si_protein_attributes.write_protein_attributes


si_proteins
----------------------
Functions associated with the ``si_protein_attributes`` module enable the reading and writing of SHEPHARD Protein files.

.. autofunction:: shephard.interfaces.si_proteins.add_proteins_from_dictionary
.. autofunction:: shephard.interfaces.si_proteins.add_proteins_from_file
.. autofunction:: shephard.interfaces.si_proteins.write_proteins


