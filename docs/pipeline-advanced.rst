.. _pipeline-advanced:

Advanced pipeline usage
=======================

Using custom models
-------------------

If you would like to try another batch correction/dimension reduction method, you can apply it to the data and run postprocessing manually.
For every model you consider, create a directory:

``my-models/``: directory with the results. In the pipeline case its called ``latent/``:

* ``model1-name/``: model name, it can be arbitrary
    * ``params.json``: model parameters, will be used to create a summary
    * ``latent_representations.csv``: for each cell name (index column), the coordinates of the latent codes
* ``model2-name/``: another directory, structured in the same manner
* ...

To help creating such directories, we created a template for the script wrapping your model at TODO

.. todo::
   Put a template for the wrapping script at GitHub (in a new ``templates/`` directory).


When the directory with different latent codes is ready, run:

.. code-block:: bash

   $ python -m cansig.run.postprocessing my-models \
                                   --expression-data original-data.hdf5
                                   --gene-sets data/pipeline-tutorial/pathways.gmt \
                                   --clusters 2 3 5  --cluster-runs 2 \
                                   --output output-dir

.. note::
   You need to specify the original HDF5 file with expression data (using the argument ``--expression-data``) to run the Gene Set Enrichment Analysis.

This will create a directory ``output-dir/`` with the same structure as the original pipeline.
If you do not wish the ``my-models/`` directory to be copied over into ``output-dir/runs/latent``, use the flag ``--no-copy``.

Using the plotting utilities
-------------------

By default, running the pipeline will result in a scatter plot being created for every configuration (i.e. for every pair of latent representation/clustering).
These scatter plots will be colored by default according to (a) the cluster assignments and (b) the batch. 

You have the option to plot using either PCA or UMAP, and you can add plotting according to pre-scored signatures (i.e. present in the adata.obs object).

The plotting utilities can be controlled through three arguments in the command line :

* ``--dim-reduction``: the dimensionality reduction method used for plotting. The options are "pca" or "umap", and defaults to PCA. 
* ``--sigcols``: optionally, the name of the signatures to plot in addition to batch and cluster assignment. 
* ``--disable-plots``: if you do not want plots to be generated, use this flag. This will skip plotting altogether.
.. note::
   The signature column names inputted with ``--sigcols`` must be in the adata.obs.columns.

Example usage to plot UMAP with additional scored signatures named ``signature1`` and ``signature2``.

.. code-block:: bash

   $ python -m cansig.run.pipeline data/pipeline-tutorial/data.hdf5 
                                   --batch batch
                                   --gene-sets MSigDB_Hallmark_2020 \
                                   --dimensions 4 6 --model-runs 1 \
                                   --clusters 2 3 5  --cluster-runs 1 \
                                   --output tutorial-output \
                                   --dim-reduction umap \
                                   --sigcols signature1 signature2

Example usage to disable plotting utilities

.. code-block:: bash

   $ python -m cansig.run.pipeline data/pipeline-tutorial/data.hdf5 
                                   --batch batch
                                   --gene-sets MSigDB_Hallmark_2020 \
                                   --dimensions 4 6 --model-runs 1 \
                                   --clusters 2 3 5  --cluster-runs 1 \
                                   --output tutorial-output \
                                   --disable-plots