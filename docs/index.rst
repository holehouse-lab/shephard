.. shephard documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SHEPHARD
=========================================================

SHEPHARD is a Python-based software framework for reading, annotating, and analyzing large protein datasets. It's objective is to make it easy for people to perform reproducible, error-free analyses of large protein datasets with arbitrarily complex sequence annotations.

The general use-case for SHEPHARD would include scenarios such as:

1. Wanting to ask large-scale statistical questions about big protein datasets. 
2. Asking how different types of sequence or protein annotations relate to one another.
3. Easily associating experimentally-generated data with extant sequence information. 


Installation
--------------

SHEPHARD is distributed via the Python packaging index (PyPI). As such, the current public release candidate can be installed using::

	pip install shephard
	
Alternatively, you can install the current bleeding-edge version from GitHub using

.. code-block:: Bash

	pip install shephard@git+git://github.com/holehouse-lab/shephard.git


This should install without issue, and once installed, SHEPHARD is available for import in any Python code you write when executed from within that :code:`conda` environment.


To test if the installation has worked, open up the Python interpreter or a Jupyter notebook and run

.. code-block:: Python

	import shephard
		
	
If this works, you should be good to go!	

A note on Python environments
.................................

To install, we strongly recommend having a :code:`conda` environment set up (the scope of :code:`conda` setup is beyond that of this documentation). :code:`conda` environments let you define a specific Python version and set of local packages that help isolate software tools away from your main system's Python version. For more information there are many `examples of conda introduction tutorials online, like this one here <https://astrobiomike.github.io/unix/conda-intro>`_. 

Demos and examples
--------------------
Once installed, SHEPHARD makes it very easy to work with large protein datasets. As an example, we have a collection of basic notebooks showing functionality located at on the `SHEPHARD GitHub supporting data page <https://github.com/holehouse-lab/shephard-colab>`_. 


Colab notebook
--------------------
If you don't want to install SHEPHARD on your local computer, we also provide several Google colab notebooks with SHEPHARD pre-installed. Check out our `colab notebook repository on GitHub <https://github.com/holehouse-lab/shephard-colab/>`_.


About
---------
SHEPHARD was developed and written by Garrett Ginell and Alex Holehouse in the Holehouse lab. For issues, bugs, or feature requests `please raise an issue on GitHub <https://github.com/holehouse-lab/shephard/issues/>`_.


Documentation index
----------------------


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   overview
   getting_started
   proteome
   protein
   domain
   site
   track
   shephard_file_types
   interfaces
   tools
   apis
   code_convention


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
