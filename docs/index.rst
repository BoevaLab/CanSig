CanSig: discovering cancer signatures
=====================================

.. warning::
   This package is currently *not* ready to be used. We hope it will be available to public use in June.

Human tumors are highly heterogeneous in their cell composition; specifically, they exhibit heterogeneity in transcriptional states of malignant cells, as has been recently discovered through single-cell RNA sequencing (scRNA-seq). Distinct states of malignant cells have been linked to variability in tumorigenic properties and resistance to anti-cancer treatment. Despite the fact that scRNA-seq data contain all necessary information to uncover shared transcriptional states of malignant cells in tumors, jointly analyzing cells from multiple cancer patients comes with its set of challenges including batch correction and accounting for patient-specific genetic background driving differences between gene expression vectors. We propose CanSig, an easy-to-use approach designed to discover known and de novo shared signatures in cancer single cells. CanSig preprocesses, integrates and analyzes scRNA-seq data to provide new signatures of shared transcriptional states and links these states to known pathways. 

.. todo::
   Add a picture of the pipeline here with the figure caption as in the preprint.

.. note::
   A preprint describing the pipeline and case studies is `now available <https://doi.org/10.1101/2022.04.14.488324>`_.


Getting started
---------------

This tutorial assumes you have access to a Unix terminal (e.g., Linux, BSD, MacOS). If you use Windows, install `Windows Subsystem for Linux <https://docs.microsoft.com/en-us/windows/wsl/install>`_ first.

When the Unix terminal is ready, you will need to `install Python <https://wiki.python.org/moin/BeginnersGuide>`_. One popular solution is `Anaconda <https://www.anaconda.com/>`_.

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

You should see a new directory ``data`` with a subdirectory ``pipeline-tutorial/`` containing two files:

* ``data.hdf5``: an artificial single-cell dataset (discussed below)
* ``pathways.gmt``: artificial gene sets, used for Gene Set Enrichment Analysis

.. todo::
   Add the description of this dataset and UMAP visualisations.

Running the analysis
^^^^^^^^^^^^^^^^^^^^

We will run the analysis pipeline with a single command:

.. code-block:: bash

   $ python -m cansig.run.pipeline data/pipeline-tutorial/data.hdf5 --gene-sets data/pipeline-tutorial/pathways.gmt \
                                   --dimensions 4 6 --model-runs 2 \
                                   --clusters 2 3 5  --cluster-runs 2 \
                                   --output tutorial-output 

This command launched the training of several models, used for batch correction and dimension reduction.
We reduce the dimension either to 4 or 6 (``--dimensions 4 6``).
For each of these dimensionalities, we run two independent runs – with different random seeds – (``--model-runs 2``), what results in 4 models in total.

Then, for each of the four models, we apply clustering into 2, 3, or 5 communities (``--clusters 2 3 5``).
Again, we repeat this step twice (``--cluster-runs 2``) using different random seeds.

Note that we specified a custom GMT file with gene sets (as the dataset we are analyzing is artificial). The default database for cancer data is TODO, although many others can be used (e.g., TODO, TODO, and TODO).

.. todo::
   Fill in the TODO items above.

The results will be available in the generated directory ``tutorial-output/`` in a few minutes.

Let's analyse its structure:

* ``summary/``: directory with a high-level summary of all the runs
    * ``nes_heatmap.pdf``: heatmap summarizing Normalized Enrichment Score of the pathways which seem to be differentially expressed
    * ``pathways.csv``: a table with the enrichment score, p-values, and q-values of all pathways
    * ``training.pdf``: plot summarizing the training curves (used to assess the convergence of dimension reduction/batch correction methods)
    * ``training_status.json``: summary whether all the dimension reduction/batch correction methods have converged
* ``runs/``: directory with all generated data, from single runs
    * ``latent/``: each model has a separate directory representing the latent representations inferred by the batch correction/dimension reduction method
    * ``cluster/``: cluster assignments (each of the latent representations file can clustered using different algorithms)
    * ``gsea/``: Gene Set Enrichment Analysis scores for each of the clustering

.. todo::
   This is a sketch. We need to adjust it accordingly, when the implementation is ready.


While we covered the most basic usage of the pipeline, more information can be obtained by running

.. code-block:: bash

   $ python -m cansig.run.pipeline --help

or by consulting the :ref:`pipeline-advanced`.


Interpreting the results
^^^^^^^^^^^^^^^^^^^^^^^^

.. todo::
   Basic discussion how to interpret the results.


For more advanced analysis, including drawing biological insights, see :ref:`interpretation`.


Tutorials
---------

* To learn more about the pipeline (parallelization, using custom models), see :ref:`pipeline-advanced`.
* For a tutorial how to interpret the results from the biological perspective, see :ref:`interpretation`.
* To learn about the preprocessing module (used to prepare the raw data into the HDF5 format), see the :ref:`preprocessing`.


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

