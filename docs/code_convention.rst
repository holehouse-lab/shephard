Code convention
=================
SHEPHARD is designed with the goal of open source contribution and additional extensions. In particular, we imagine the core to be (relatively) stable, but the APIs module offers a means for the development of novel tools to enable analysis to be performed on SHEPHARD Proteomes using third-party tools and libraries.

With this in mind, this file contains a set of recommendations for developing, adding, or contributing to SHEPHARD. This page has several subsections which offer guidance on how to write functions, the philosophical goals of the code layout, and finally some specific guidance on writing new APIs.


Format for function signatures
--------------------------------

First things first - if you want to define a new function, it's important to follow the standard format used for all function docstrings.

The general format for function signatures in SHEPHARD show below.

.. code-block:: python

    def function_name():

	"""
        Summary of function goals, including what data (if any) are changed
        by the function. This should be succinct but clear. Importantly, text
        should not exceed 80 charactes so it's easily renderable inline in a 
        Jupyter notebook. The exception to this is if formatting will break 
        the Sphinx restructured text parsing. For example::

            A code block can be shown inline using two :: characters

        or bullet pointed info can be shared using

        * This will render as a bullet point, but if we add a line break it would disrupt the parsing so bullet-pointed info should be kept on a single line
        
        Finally, it's often useful to discuss pitfals or possible errors
        that can be raised. 

        Parameters
        -----------
    
        <parameter_name> : <type>  (default = <default value>)
            Description of the parameter

        Returns
        -----------
        <return type>
            Description of the return type and what it means
	"""
    

Sometimes a function will be able to take multiple different input values. In this case, we use the following syntax

.. code-block:: python

    """
    Parameters
    -----------
    
    <parameter_name> : <type1>, <type2> (default = <default vaue>)
    """


Code philosophy
---------------------
SHEPHARD is not a sequence/protein analysis package, but instead a package that makes it easy to load, annotate, access, and export sequence analysis. It is divided into several distinct modules, which each serve a specific set of purposes, as outlined below. While we welcome additional contributions to the code, please consider the goals of the different modules before adding and making a pull request. Additionally, we reserve the right to re-organzie submitted code requests if we feel that a feature would be better suited elsewhere.

In general, new code added should consist of stand-alone functions. 

Main data types
.................
The main data types that SHEPHARD supports (Proteomes, Proteins, Sites, Domains, and Tracks) are deliberately kept as bare bones data-structures, as opposed to classes that incorporate specific functionality. Class-associated functionality is limited to ease of access to or organizing data, as opposed to any kind of sequence analysis. More complex functionality for converting data types or performing complex analysis should be coded in the ``shephard.tools`` module (described below), but this too should avoid any kind of actual analysis code and should be limited to (more complex) functions for manipulating data.


Interfaces (``shephard.interfaces``)
.......................................

Modules in the interface module (shephard.interfaces) are strictly controlled and define the I/O functionality for reading/writing SHEPHARD compliant files. Code changes offered here should be bug fixes, but should NOT be used to support additional data types (see ``shephard.apis``). 


APIs (``shephard.apis``)
.......................................

Modules in the APIs (Application Programming Interfaces) module offer the oppertunity to generate code that works with SHEPHARD in a standalone way to annotate SHEPHARD objects. Initially, APIs are limited to specific tools or data-types that are well-defined, but, as we expand SHEPHARD additional APIs can and will be added here to provide additional 'built in' functionality. That said, to limited SHEPHARD's required footprint, dependencies in the apis will not be considered core dependencies and will not be added to the set of required packages for SHEPHARD installation. With this in mind, we will include an import check and warning to trigger download of additional packages or tools as needed. This decision avoids a scenario where SHEPHARD's installation becomes tethered to a large number of distinct and possibly incompatible packages (aka dependency hell).

APIs are likely the main place where new functionality could be contributed to. As such, we have a specific section at the bottom of this document offing a brief guide on how to create an API that meets the expectations for SHEPHARD. Importantly, any new code added must include corresponding tests in the ``shephard.tests`` module, as described below in the **Writing Tests** section.


Tools (``shephard.tools``)
...............................
In general, functions defined in tools modules (shephard.tools) should be stateless and non-mutating. What this means is they should:

1. Take input data only
2. Not change input data passed directly, but instead return a type that can be used to update stateful objects (e.g. Proteomes, Proteins etc).


Miscellaneous modules
-----------------------

In addition to the major module classes outlined above, there are several additional modules that provide generic functionality.


general_utilities (``shephard.general_utilities``)
.....................................................

The general utilities module provides stateless data manipulation functions for doing a variety of non-specific work. This includes data type conversion, simple mathematical operations, and sanity checking functions. Any function that (broadly) carries out a generic Python-associated function can be included here. The functions here could in principle be used by other packages as well, and are in no-way meant to be SHEPHARD specific.


sequence_utilities (``shephard.sequence_utilities``)
.......................................................

The sequence utilities module is, analagous to the general utilities module, a place for a set of stateless functions that perform sequence manipulation. These functions are meant to be limited to SHEPHARD, although in principle like those found in general utilities could be useful outside of SHEPHARD. However, they mostly are included to solve SHEPHARD-specific generic sequence-associated problems. For a broader set of sequence manupulation tools, see the ``shephard.tools.sequence_tools`` module.


sequence_tools (``shephard.tools.sequence_tools``)
......................................................

The sequence tools module contains a set of sequence manipulation functions (where sequences here are just strings) that may be of general use, both inside and outside SHEPHARD. Just as the ``shephard.domain_tools`` is designed to work in a domain-focusse way, the sequence tools module is meant to work with sequence (`str`) focussed way. 


exceptions (``shephard.exceptions``)
........................................

The SHEPHARD exceptions class allows customizable exceptions to be defined. In general we are trying to keep these exceptions somewhat limited in number, but they enable more specific error handling in complex pipelines.


Contributing to SHEPHARD
--------------------------



How to write a new API
..........................
TO DO


How to write tests
..........................

For any new code added to SHEPHARD, a collection of tests to ensure appropriate behaviour when both valid and invalid data are passed is essential. Tests in SHEPHARD are dealt with using PyTest, and test data (to be read in) can be stored in the ``shephard/data/test_data`` directory, which SHEPHARD provides easy access to via an internal function. As an example, below we walk through how to create a new test for some new hypothetical function defined in ``tools.domain_tools`` that takes in a SHEPHARD domain and returns the domain length.


Firstly, we'll define our new function to make this complete transparent

.. code-block:: python

    def calculate_domain_length(domain):
        """
	Function that returns a domain's length
	
	Parameters
	-------------
	domain : Domain
	    Domain object in question

	Returns
	------------
	int
	    Returns the length of the domain
	"""

	return len(domain)


Having defined this function, we move into the directory ``shephard/data/test_data`. We could EITHER create new files for our tests, or take avantage of some of the existing test data. To make things simple and reproducible we'll use two files that already exist:

* ``testset_1.fasta`` - which contains a set of FASTA files in UniProt format, and 

* ``TS1_domains_idr.tsv`` - where TS1 = testset_1, and is a SHEPHARD-compliant Domains file with a set of IDRs

Having established the data we're going to use for our test, we next need to create a test file. To do this we move into ``shephard/tests/`` and we create a new file called ``test_calculate_domain_lengths.py``. A couple of points - in general we recommend creating new files for any new features added, in part because this makes debugging easier. Secondly, test modules MUST start with the name ``test_``.

Once we've created our new file (``test_calculated_domain_lengths.py``) we next are going to create a function that tests our function, and put the following code in that function.

.. code-block:: python

    # file test_calculated_domain_lengths.py

    # must import pytest
    import pytest

    # import additional modules needed for this code
    import shephard
    from shephard.apis import uniprot  
    from shephard import interfaces

    # import the domain_tools module which we're going to test
    from shephard.tools import domain_tools

    def test_first_test():

        ## The first part of this function involves reading in a Proteome object
	#  and annotating with a set of domains
    

        # this line uses shephards internal get_data function to build the full
	# path to the directory where our test data is
	test_data_dir = shephard.get_data('test_data')

	# define filenames - note we use the names we defined above
	fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
	domain_file = '%s/%s' % (test_data_dir, 'TS1_domains_idr.tsv')

	# create a new Proteome object and annotate with the domains
	# in the domain_file
	P = uniprot.uniprot_fasta_to_proteome(fasta_file)
	interfaces.si_domains.add_domains_from_file(P, domain_file)


    # we manually looked in the TS1_domains.idr.tsv file and found a protein
    # where the 1st domain length is KNOWN. The test we're going to write is to
    # to make sure the function can reproduce this 'known' domain length.
    TEST_ID='O00401'

    # get the 1st domain in the protein defined by the TEST_ID ID - note we're indexing
    # into the 0th element in this protein to get the 1st domain.
    domain = P.protein(TEST_ID).domains[0]

    # this is where the test actually happens: we want to use the assert statement 
    # to check that if we run our new function it returns the expected value (20, in 
    # this case).
    assert domain_tools.calculate_domain_length(domain) == 20

The code above takes advantage of the fact that we KNOW the first domain in the protein defined by the ID **O00401** should be 20 residues long. Tests are based on ensuring a function gives an expected outcome, and, as with any expected outcome, you should make sure you know the outcome ahead of time to write a good test! 

As with the file names, the test functions must start with ``test_`` - all other functions will be ignored but can be useful as helper functions.

Finally, to run this test, we can execute ``pytest`` from the command-line. Running:

.. code-block:: bash

    pytest

Or to get information on each individual test run:

.. code-block:: bash

    pytest -v


by itself in the test directory will execute ALL the tests (including your new test). To run the tests only for a specific file, you can use the syntax

.. code-block:: bash

    pytest  test_calculate_domain_lengths.py

or 

.. code-block:: bash

    pytest -v test_calculate_domain_lengths.py


This is the basics for getting tests up and running, but, of course, right now this function is quite vulnerable to invalid input. For example, passing the word 'CAT' would return 3, but 'CAT' is not a valid domain. As such, it behoves us to write user-facing functions that do some sanity checking, and also ensure that sanity check works. Specific examples of this go beyond the scope of this basic tutorial, but pytest provides a syntax for checking a function raises a specific type of exception:


.. code-block:: python

    from shephard.exceptions import ProteinException

    with pytest.raises(ProteinException):
        # some function you EXCEPT to raise a ProteinException

  
In the psuedocode above, we import the ProteinException class (an exception defined by SHEPHARD) and then using the ``with`` syntax run a statement which says, *"if the function raises ProteinException then this passes, else this is considered a failure"*. This type of exception-checking can be EXTREMELY useful for ensuring your code is robust to unexpected input.


