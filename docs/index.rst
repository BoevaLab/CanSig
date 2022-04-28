CanSig: discovering cancer signatures
=====================================

.. warning::
   This package is currently *not* ready to be used. We hope it will be available to public use in June.

Human tumors are highly heterogeneous in their cell composition; specifically, they exhibit heterogeneity in transcriptional states of malignant cells, as has been recently discovered through single-cell RNA sequencing (scRNA-seq). Distinct states of malignant cells have been linked to variability in tumorigenic properties and resistance to anti-cancer treatment. Despite the fact that scRNA-seq data contain all necessary information to uncover shared transcriptional states of malignant cells in tumors, jointly analyzing cells from multiple cancer patients comes with its set of challenges including batch correction and accounting for patient-specific genetic background driving differences between gene expression vectors. We propose CanSig, an easy-to-use approach designed to discover known and de novo shared signatures in cancer single cells. CanSig preprocesses, integrates and analyzes scRNA-seq data to provide new signatures of shared transcriptional states and links these states to known pathways. 


.. note::
   A preprint describing the pipeline and case studies is `now available <https://doi.org/10.1101/2022.04.14.488324>`_.


Getting started
---------------

This tutorial assumes you have access to a Unix terminal (e.g., Linux, BSD, MacOS). If you use Windows, install `Windows Subsystem for Linux <https://docs.microsoft.com/en-us/windows/wsl/install>`_ first.

Installation
^^^^^^^^^^^^

Install the package:

.. code-block:: bash

   $ pip install cansig

Example dataset
^^^^^^^^^^^^^^^

We will analyse an artificial single-cell dataset, stored using the HDF5 format. We will download it using an auxiliary script:

.. code-block:: bash

   $ python -m cansig.run.download HDF5

You should see a new directory ``data`` containing the dataset ``example.hdf5``.

.. todo::
   Add the description of this dataset and UMAP visualisations.

Running the analysis
^^^^^^^^^^^^^^^^^^^^

We will run the analysis pipeline with a single command:

.. code-block:: bash

   $ python -m cansig.run.pipeline data/example.hdf5 --dimensions 4 6 --clusters 2 3 5 --model-runs 2 --cluster-runs 2 --output tutorial-output

This command launched the training of several models, used for batch correction and dimension reduction.
We reduce the dimension either to 4 or 6 (``--dimensions 4 6``).
For each of these dimensionalities, we run two independent runs (``--model-runs 2``), what results in 4 models in total.

Then, for each of the four models, we apply clustering into 2, 3, or 5 communities (``--clusters 2 3 5``). Again, we repeat this step twice (``--cluster-runs 2``).

The results are available in the generated directory ``tutorial-output``.

Let's analyse its structure:

* ``summary/``
    * ``nes_plot.pdf``
    * ``pathways.csv``
* ``runs/``
    * ``latent/``
    * ``cluster/``
    * ``gsea/``


.. note::
   Currently, we use `scVI <https://github.com/scverse/scvi-tools>`_ for batch correction and dimension reduction and `Leiden clustering <https://doi.org/10.1038/s41598-019-41695-z>`_.
   For more advanced tutorial, instructing how to use different models, consult :ref:`pipeline-advanced`.

Interpreting the results
^^^^^^^^^^^^^^^^^^^^^^^^

.. todo::
   Basic discussion how to interpret the results.

.. note::
   For more advanced analysis, see the :ref:`interpretation` tutorial.

Tutorials
---------

* To learn more about the pipeline (parallelization, using custom models), see :ref:`pipeline-advanced`.
* For a tutorial how to interpret the results from the biological perspective, see :ref:`interpretation`.
* To learn about the preprocessing module (used to prepare the raw data into the HDF5 format), see the :ref:`preprocessing`.

.. todo::
   Summary of tutorials. Tutorial on interpretation.


Contributing
------------

For the contribution guide and instructions for new developers, see :ref:`contribution-guide`.


Indices and tables
------------------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   preprocessing
   pipeline-advanced
   interpretation
   contributing
   source/modules

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

