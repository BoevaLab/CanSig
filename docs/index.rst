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

For the contribution guide and instructions for new developers, see :ref:`contribution-guide`.


Indices and tables
------------------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   contributing
   source/modules
   

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

