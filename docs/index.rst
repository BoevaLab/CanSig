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

.. note:: CanSig requires Python 3.8 (or a most recent version).

Installation
^^^^^^^^^^^^

Install the package:

.. code-block:: bash

   $ pip install cansig

Example dataset
^^^^^^^^^^^^^^^

To demonstrate the pipeline usage, we will analyse a small esophageal dataset, which has been already preprocessed (for more information about the dataset and the data curation, see :ref:`preprocessing`).

We will download it using an auxiliary script:

.. code-block:: bash

   $ python -m cansig.run.download PIPELINE

You should see a new directory ``data`` with a file ``pipeline-tutorial.h5ad``.

.. todo::
   Add the description of this dataset and UMAP visualisations.

Running the analysis
^^^^^^^^^^^^^^^^^^^^

.. note:: If you are familiar with Python, you may prefer running the analysis in a custom manner. See the :ref:`pipeline-advanced`.

We will run the analysis pipeline with a single command:

.. code-block:: bash

   $ python -m cansig.run.pipeline data/pipeline-tutorial/data.hdf5 
                                   --batch batch
                                   --gene-sets MSigDB_Hallmark_2020 \
                                   --dimensions 4 6 --model-runs 1 \
                                   --clusters 2 3 5  --cluster-runs 1 \
                                   --output tutorial-output 

This command launched the training of several models, used for data integration (batch correction and dimension reduction).
We marked that the batch column in the dataset is named ``batch`` and we will reduce the dimension either to 4 or 6 (``--dimensions 4 6``).
For each of these dimensionalities, we do a single run (``--model-runs 1``). Increasing this number will train more models, differing in random seeds.


.. note::

   In case of additional continuous or discrete variables, which may confound the expression (e.g., the sequencing technology used), one could specify it via the ``--continuous-covariates`` and ``--discrete-covariates`` flags.
   
   For more information, run 

   ``$ python -m cansig.run.pipeline --help``


Then, for each of the four models, we apply clustering into 2, 3, or 5 communities (``--clusters 2 3 5``).
Again, one can repeat this steps, using different random seeds, by increasing ``--cluster-runs``.

We specified the database to be used (``--gene-sets MSigDB_Hallmark_2020``).
This is one of many databases available online. Alternatively, one can provide a path to the GMT file.

.. todo::
   Fill in the TODO items above.

The results will be available in the generated directory ``tutorial-output/`` in a few minutes to hours depending on the size of the dataset.

Let's analyse its structure:

* ``integration/``: each model has a separate directory representing the latent representations inferred by the batch correction/dimension reduction method
   * ``{rundir}/integration-settings.json``: summary of the parameters used for the integration method
   * ``{rundir}/latent-representations.csv``: coordinates in the latent space inferred by the integration method
* ``postprocessing/``: each clustering configuration on each model has a separate directory with the postprocessing information saved
   * ``{rundir}/cluster-labels.csv``: cluster assignments for this run 
   * ``{rundir}/cluster-settings.json``: summary of the parameters used for the clustering step
   * ``{rundir}/gsea-dataframe.csv``: results of GSEA run in this setting (full table, including non-significant pathways and negative enrichment score pathways)
   * ``{rundir}/gsea-settings.json``: summary of the parameters used for the GSEA run
   * ``{rundir}/integration-settings.json``: summary of the parameters used for the integration method
   * ``{rundir}/latent_space_dimred.png``: by default, scatter plot of the PCA on the latent space colored according to batch and clustering (see advanced tutorial for more details)
* ``heatmap.pdf``: heatmap summarizing Normalized Enrichment Score of the pathways which seem to be differentially expressed over all runs indicated


While we covered the most basic usage of the pipeline, more information can be obtained by running

.. code-block:: bash

   $ python -m cansig.run.pipeline --help

or by consulting the :ref:`pipeline-advanced`. 
Options include running your own integration model or controlling the plotting utilities.


Interpreting the results
^^^^^^^^^^^^^^^^^^^^^^^^

.. todo::
   Basic discussion how to interpret the results.


For more advanced analysis, including drawing biological insights, see :ref:`interpretation`.


Tutorials
---------

* To learn more about the pipeline (parallelization, using custom models, plotting), see :ref:`pipeline-advanced`.
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

