SHEPHARD files
=================

SHEPHARD defines a set of well-defined files for reading in or writing out data into SHEPHARD objects.

The interfaces package is a sub-package within SHEPHARD that deals with reading in **SHEPHARD** formatted input files. These are files with a specific format that were developed for SHEPHARD. The interfaces mean that, as long as you can write your data in a format that complies with a track, site, domain, or protein_attribute file, you can be sure it will be correctly read into SHEPHARD and is then accessible within the larger framework. 

Each line in an input file corresponds to a single piece of information that maps to a specific protein, but multiple lines can map to the same protein. 

The interfaces modules (:code:`si_sites`, :code:`si_domains`, :code:`si_tracks`, :code:`si_protein_attributes` ) contain two user-facing functions; one that lets the user read in data from a file, and another than lets a user read in data from an input dictionary. This means that sites, domains, tracks and attributes can be loaded from disk or in real-time as other pipelines are run.


Sites files
------------

Sites files are tab-separated files where each line represents a different site. Sites represent position-specific information, and are useful for things like mutations, PTMs, binding sites etc.

Each line has five or more columns and uses the following format:

.. code-block:: 
		   
   Unique_ID    position    site type    symbol    value    [   key1:value1    key2:value2     ...keyn:valuen   ]

Note each entry is separated by a tab character.

Here:

* :code:`Unique_ID` is a unique identifier for a protein that will be used to cross reference this site to a specific protein in the Proteome. We recommend uniprot IDs for naturally occurring proteins.
* :code:`position` is the position in the protein this site is found (noting that the first position in a protein is position 1)
* :code:`site type` is a free-form string that describes the type of the site. This can be anything, and is used for site selection in SHEPHARD
* :code:`symbol` is some kind of symbolic information that is ascribed to the site. Because this is a required field, if no symbolic information makes sense simply include a dash here.
* :code:`value` is some kind of numerical information that is ascribed to the site. Because this is a required field, if no numerical information makes sense simply set this to 1
* :code:`key:value` the remainder of the file is made up of key:value pairs which allow arbitrary metadata to be associated with each site. The colon ':' character splits two strings, where the first is a key and second a value. These become accessible via the site's :code:`attributes` properties.


Domain files
------------

Domain files are tab-separated files where each line represents a different domain. Domains represent contiguous regions along the sequence and are useful for defining folded regions, disordered regions, or any other kind of annotation where a discrete region can be assigned.

Each line has four or more columns and uses the following format:

.. code-block:: 
		   
   Unique_ID     start    end    domain_type    [   key1:value1    key2:value2     ...keyn:valuen  ]

Here:

* :code:`Unique_ID` is a unique identifier for a protein that will be used to cross reference this domain to a specific protein in the Proteome. We recommend uniprot IDs for naturally occurring proteins.
* :code:`start` is the start position in the protein where the domain is found (noting that the first position in a protein is position 1)
* :code:`end` is the end position in the protein where the domain is found (noting that the first position in a protein is position 1)
* :code:`domain_type` is a free-form string that describes the type of the domain. This can be anything, and is used for domain selection in SHEPHARD
* :code:`key:value` the remainder of the file is made up of key:value pairs which allow arbitrary metadata to be associated with each domain. The colon ':' character splits two strings, where the first is a key and the second a value. These become accessible via the domain's :code:`attributes` properties. key:value attributes are OPTIONAL


Track files
------------

Track files are tab-separated files where each line represents a distinct track. Tracks represent vectorial information that maps a specific number or symbol to a protein sequence on a per-residue basis. Each track has a variable number of columns, that corresponds to 2 + the number of residues in the protein. Tracks are useful where some type of analysis reveals a continuous output that can be projected along the sequence. If, on the other hand, your analysis reveals discrete regions or positions, Domains or Sites might be more appropriate, respectively.

The format of a track file is as follows: 

.. code-block:: 
		   
   Unique_ID     track_name    v1     v2     v2    ...   vn

Here:

* :code:`Unique_ID` is a unique identifier for a protein that will be used to cross reference this track to a specific protein in the Proteome. We recommend uniprot IDs for naturally occurring proteins. 
* :code:`track_name` is a free-form string that describes the type of the track. This can be anything, and is used for track selection in SHEPHARD


The remainder of the file (:code:`v1`, :code:`v2`, :code:`vn`) are values or symbols that should map to each residue, such that each residue has a value or symbol in the track file. When reading in a track file, an error will occur if the number of symbols/values does not match the protein length.

Tracks must be EITHER all numerical or all symbolic, and when read in the user has to specify which they are. This reflects the fact that Track objects are either symbolic- or value-based. 



Protein attribute files
-------------------------

Domain attribute files are tab-separated files where each line represents a collection of attributes that can be assigned to a given protein. Protein attributes are useful for metadata that has no specific sequence-positional relevance, such as sub-cellular localization, GO annotation, concentration etc.

 The format of a protein attribute file is as follows: 

.. code-block:: 
		   
   Unique_ID    key1:value1   [ key2:value2     ...keyn:valuen   ]

Here:

* :code:`Unique_ID` is a unique identifier for a protein that will be used to cross reference this set of attributes to a specific protein in the Proteome. We recommend uniprot IDs for naturally occurring proteins. 
* :code:`key:value` the remainder of the file is made up of key:value pairs which allow arbitrary metadata to be associated with each protein. The colon ':' character splits two strings, where the first is a key and second a value. These become accessible via the Protein's :code:`attributes` properties. 


