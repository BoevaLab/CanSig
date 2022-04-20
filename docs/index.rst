.. CanSig documentation master file, created by
   sphinx-quickstart on Wed Apr 20 13:39:48 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CanSig: discovering cancer signatures
=====================================


Human tumors are highly heterogeneous in their cell composition; specifically, they exhibit heterogeneity in transcriptional states of malignant cells, as has been recently discovered through single-cell RNA sequencing (scRNA-seq). Distinct states of malignant cells have been linked to variability in tumorigenic properties and resistance to anti-cancer treatment. Despite the fact that scRNA-seq data contain all necessary information to uncover shared transcriptional states of malignant cells in tumors, jointly analyzing cells from multiple cancer patients comes with its set of challenges including batch correction and accounting for patient-specific genetic background driving differences between gene expression vectors. We propose CanSig, an easy-to-use approach designed to discover known and de novo shared signatures in cancer single cells. CanSig preprocesses, integrates and analyzes scRNA-seq data to provide new signatures of shared transcriptional states and links these states to known pathways. 


.. note::

   A preprint describing the pipeline and case studies is `now available <https://doi.org/10.1101/2022.04.14.488324>`_.


Getting started
---------------

To install the package,

.. todo::

   When the package is on PyPI, add installation instructions.

.. todo::

   Add the instructions how to analyse a small (e.g., simulated) data set by calling one command.


Tutorials
---------

.. todo::

   Summary of tutorials. In particular, we should show how to use the preprocessing functions.


Contributing
------------

Install the module in editable mode and with the dependencies needed to run the test suite:

.. code-block:: bash

   $ pip install -e ".[test]"

At this point you should be able to run the unit test suite:

.. code-block:: bash

   $ pytest


Developer tools
^^^^^^^^^^^^^^^

As a developer, you will however need more tools. You can install them by running

.. code-block:: bash

   $ pip install -r requirements-dev.txt


.. todo::

   Add the description of code quality tools and pre-commit.
   

Working with documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

We use `Sphinx <https://www.sphinx-doc.org>`_  and `ReStructuredText <https://docutils.sourceforge.io/rst.html>`_. For a reference, see e.g., `this cheatsheet <https://sphinx-tutorial.readthedocs.io/cheatsheet/>`_.

Whenever a Pull Request is merged into the ``main`` branch, a new version of the documentation is automatically built and deployed.
Although docstrings (in Python code) are assembled automatically, high-level documentation (including tutorials), needs to be written manually.
We store the source code of the documentation in the ``docs`` directory. To build it locally, use

.. code-block:: bash

   $ cd docs
   $ sphinx-apidoc -o source/ ../cansig
   $ make html

In the ``_build`` directory, you should see ``index.html`` file. Open it with a webbrowser of your choice.


Indices and tables
------------------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   source/cansig
   source/modules

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

