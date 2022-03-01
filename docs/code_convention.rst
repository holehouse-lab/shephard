Code convention
=================
SHEPHARD is designed with the goal of open source contribution and additional extensions. In particular, we imagine the core to be (relatively) stable, but the APIs module offers a means for the development of novel tools to enable analysis to be performed on SHEPHARD Proteomes using third-party tools and libraries.

With this in mind, this file contains a set of recommendations for developing, adding, or contributing to SHEPHARD. This page has several subsections which offer guidance on how to write functions, the philosophical goals of the code layout, and finally some specific guidance on writing new APIs.


Format for function signatures
--------------------------------

First things first - if you want to define a new function, it's important to follow the standard format used for all function docstrings.

The general format for function signatures in SHEPHARD show below.

.. code-block:: bash

    def function_name():

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
    

Sometimes a function will be able to take multiple different input values. In this case, we use the following syntax

.. code-block:: bash

    Parameters
    -----------
    
    <parameter_name> : <type1>, <type2> (default = <default vaue>)


Code philosophy
---------------------
SHEPHARD is not a sequence/protein analysis package, but instead a package that makes it easy to load, annotate, access, and export sequence analysis. It is divided into several distinct modules, which each serve a specific set of purposes, as outlined below. While we welcome additional contributions to the code, please consider the goals of the different modules before adding and making a pull request. Additionally, we reserve the right to re-organzie submitted code requests if we feel that a feature would be better suited elsewhere.

In general, new code added should consist of stand-alone functions. 

Main data types
.................
The main data types that SHEPHARD supports (Proteomes, Proteins, Sites, Domains, and Tracks) are deliberately kept as bare bones data-structures, as opposed to classes that incorporate specific functionality. Class-associated functionality is limited to ease of access to or organizing data, as opposed to any kind of sequence analysis. More complex functionality for converting data types or performing complex analysis should be coded in the ``shephard.tools`` module (described below), but this too should avoid any kind of actual analysis code and should be limited to (more complex) functions for manipulating data.


Interfaces (``shephard.interfaces``)
...................................
Modules in the interface module (shephard.interfaces) are strictly controlled and define the I/O functionality for reading/writing SHEPHARD compliant files. Code changes offered here should be bug fixes, but should NOT be used to support additional data types (see ``shephard.apis``). 


APIs (``shephard.apis``)
.......................
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
.................................
The general utilities module provides stateless data manipulation functions for doing a variety of non-specific work. This includes data type conversion, simple mathematical operations, and sanity checking functions. Any function that (broadly) carries out a generic Python-associated function can be included here. The functions here could in principle be used by other packages as well, and are in no-way meant to be SHEPHARD specific.


sequence_utilities (``shephard.sequence_utilities``)
.................................
The sequence utilities module is, analagous to the general utilities module, a place for a set of stateless functions that perform sequence manipulation. These functions are meant to be limited to SHEPHARD, although in principle like those found in general utilities could be useful outside of SHEPHARD. However, they mostly are included to solve SHEPHARD-specific generic sequence-associated problems. For a broader set of sequence manupulation tools, see the ``shephard.tools.sequence_tools`` module.


sequence_tools (``shephard.tools.sequence_tools``)
....................................
The sequence tools module contains a set of sequence manipulation functions (where sequences here are just strings) that may be of general use, both inside and outside SHEPHARD. Just as the ``shephard.domain_tools`` is designed to work in a domain-focusse way, the sequence tools module is meant to work with sequence (`str`) focussed way. 


exceptions (``shephard.exceptions``)
.................................
The SHEPHARD exceptions class allows customizable exceptions to be defined. In general we are trying to keep these exceptions somewhat limited in number, but they enable more specific error handling in complex pipelines.


Contributing to SHEPHARD
--------------------------



How to write a new API
..........................
TO DO


How to write tests
..........................
TO DO
