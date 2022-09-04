Interfaces
=================

Overview
------------

Interfaces define functions that enable the reading or writing of data into or out of SHEPHARD proteomes. In particular, these functions operate by taking a **Proteome** object and then either annotating the proteins in that Proteome object with Tracks, Domains, Sites, or Attributes, OR writing Tracks, Domains, Sites or Attributes to file.

Example
------------
By way of a simple example, here we use the Domains interface package.

.. code-block:: python

	from shephard.apis import fasta
	from shephard import interfaces
	
	# create a new Proteome from a FASTA file
	small_proteom = fasta.fasta_to_proteome('sequences.fasta')
	
	# use the interfaces package to annotate the Proteome object with the domains 
	# in the file `DNA_binding_domains.tsv`
	interfaces.si_domains.add_domains_from_file(small_proteom, 'DNA_binding_domains.tsv')
	
By using interfaces, we ensure that, as long as you can write your data in a format that complies with a Track, Site, Domain, or Protein_attribute file, you can be sure it will be correctly read into SHEPHARD and is then accessible within the larger framework. 


si_sites
------------
Functions associated with the ``si_sites`` module enable the reading and writing of SHEPHARD Sites files.

.. autofunction:: shephard.interfaces.si_sites.add_sites_from_file
.. autofunction:: shephard.interfaces.si_sites.add_sites_from_dictionary
.. autofunction:: shephard.interfaces.si_sites.write_sites
.. autofunction:: shephard.interfaces.si_sites.write_sites_from_list


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
.. autofunction:: shephard.interfaces.si_tracks.write_tracks
.. autofunction:: shephard.interfaces.si_tracks.write_tracks_from_list


si_protein_attributes
----------------------
Functions associated with the ``si_protein_attributes`` module enable the reading and writing of SHEPHARD Protein attribute files.

.. autofunction:: shephard.interfaces.si_protein_attributes.add_protein_attributes_from_dictionary
.. autofunction:: shephard.interfaces.si_protein_attributes.add_protein_attributes_from_file
.. autofunction:: shephard.interfaces.si_protein_attributes.write_protein_attributes


si_proteins
----------------------
Functions associated with the ``si_protein_attributes`` module enable the reading and writing of SHEPHARD Protein files. While we include this for completeness, our general recommendation is to use FASTA files for protein information, and then write protein attributes out as separate protein attributes files. The reason for this is that this ensures easy readability of both protein sequence information and protein annotation information.

.. autofunction:: shephard.interfaces.si_proteins.add_proteins_from_dictionary
.. autofunction:: shephard.interfaces.si_proteins.add_proteins_from_file
.. autofunction:: shephard.interfaces.si_proteins.write_proteins


