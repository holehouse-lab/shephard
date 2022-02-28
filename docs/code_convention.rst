Code convention
=================

Format for function signatures
------------------------------

The general format for function signatures in SHEPHARD show below.

.. code-block:: bash

    Parameters
    -----------
    
    <parameter_name> : <type>, default = <default value>
        Description of the parameter

    Returns
    -----------
    <return type>
        Description of the return type and what it means
    

.. code-block:: bash

    Parameters
    -----------
    
    <parameter_name> : { <type1>, <type2> }, default = <default vaue>
        Curly braces can be used to define if a parameter can be multiple 
	possible types. Text that describes functions or parameters should
	not exceed 80 characters in width for nice formatting when called
	from a Jupyter notebook.


Code philosophy
---------------------





In general, functions defined in the tools module  of the main data classes (Proteome, Proteins, Tracks, Domains, Sites) should be stateless and non-mutating. What this means is they should:

1. Take input data only
2. Not change input data passed directly, but instead return a type that can be used to update stateful objects (e.g. Proteomes, Proteins etc).

The exception to this 
