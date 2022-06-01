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
----------------------------

By default, running the pipeline will result in a scatter plot being created for every configuration (i.e. for every pair of latent representation/clustering).
These scatter plots will be colored by default according to (a) the cluster assignments and (b) the batch. 

You have the option to plot using either PCA, UMAP, or both (which will result in a UMAP with a PCA inset as seen in the publication). 
You can add plotting according to pre-scored signatures (i.e. present in the adata.obs object).

The plotting utilities can be controlled through three arguments in the command line :

* ``--dim-reduction``: the dimensionality reduction method used for plotting. The options are "pca", "umap", or "both", and defaults to PCA. 
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

Example usage to plot the UMAP with PCA inset with additional scored signatures named ``signature1`` and ``signature2``.

.. code-block:: bash

   $ python -m cansig.run.pipeline data/pipeline-tutorial/data.hdf5 
                                   --batch batch
                                   --gene-sets MSigDB_Hallmark_2020 \
                                   --dimensions 4 6 --model-runs 1 \
                                   --clusters 2 3 5  --cluster-runs 1 \
                                   --output tutorial-output \
                                   --dim-reduction both \
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

Saving signatures and scoring cells
----------------------------

CanSig finds signatures in cancer data: these signatures are defined as the most overexpressed genes in a biologically meaningful part of the subspace (approximated by a cluster).
By default, running the pipeline will result the results of the differential gene expression analysis to be saved for each cluster (for more information on interpretation, see :ref:`interpretation`), as wells as the score for each cell for each of the signatures (defined as the 200 most positively differentially expressed genes in the cluster) and the correlation between the signatures.

The results associated with signatures and scoring can be controlled through the three following command line arguments:
* ``--disable-signatures``: if you do not want any of the results linked to signatures to be saved (ie differential gene expression analysis for each cluster, scores of each signature for each cell and correlation between signatures), then input this flag.
* ``--ngenessig``: number of genes that define a signature associated with a cluster. The cells will be scored using the n top positively differentially expressed genes in the cluster, with n being the minimal value between the total number of differentially expressed genes and ngenessig inputted. Defaults to 200.
* ``--corrmethod``: the method used to compute the correlation between the signatures; can be spearman or pearson. Defaults to pearson.

Example usage to compute scores on cells using the 100 top most differentially expressed genes as a signature and using spearman to compute the correlation between signatures.

.. code-block:: bash

   $ python -m cansig.run.pipeline data/pipeline-tutorial/data.hdf5 
                                   --batch batch
                                   --gene-sets MSigDB_Hallmark_2020 \
                                   --dimensions 4 6 --model-runs 1 \
                                   --clusters 2 3 5  --cluster-runs 1 \
                                   --output tutorial-output \
                                   --ngenessig 100 \
                                   --corrmethod spearman 

Example usage to disable saving any results linked to signatures 

.. code-block:: bash

   $ python -m cansig.run.pipeline data/pipeline-tutorial/data.hdf5 
                                   --batch batch
                                   --gene-sets MSigDB_Hallmark_2020 \
                                   --dimensions 4 6 --model-runs 1 \
                                   --clusters 2 3 5  --cluster-runs 1 \
                                   --output tutorial-output \
                                   --disable-signatures                                 

Running differential CNV analysis
---------------------------------

You have the option to perform differential CNV analysis. With original CNV calls, the will output differential CNV regions between each cluster and the rest, and information about the percentage of gains/losses in the cluster and in the rest.
This module is deactivated by default. There are two main ways to run this analysis: the first assumes that you are using a data object that has been obtained using our preprocessing module (see :ref:`preprocessing`), the second can be run if provided with a external discretized CNV calling, even if the data object has not been obtained through our preprocessing module.

.. note::
   the data provided for the tutorial has been processed using our preprocessing module, and can be thus used for differential CNV analysis assuming so.

Differential CNV analysis for data preprocessed with our module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This calling assumes the data was preprocessing using our preprocessing module (see :ref:`preprocessing`).
This means that the data object you provide will contain the following:

   - "X_cnv" in data.obsm: the CNV called using our preprocessing module
   - "chromosome" in data.var.columns: the chromosome to which the gene belongs
   - "cnv_called" in data.var.columns: if this gene was used for the infercnv call (see
         `cansig._preprocessing` for more details on the CNV calling procedure)
   - "start" in data.var.columns: the start position of the gene on the chromosome
   - "cnv" in data.uns: a summary of the infercnv run
   - "chr_pos" in data.uns["cnv"]: a dictionary containing the mapping between the chromosome and
         the index of the regions in the cnv array

The analysis can be controlled through three arguments in the command line:

* ``--diffcnv``: this flag needs to be added for the differential CNV analysis to be performed. If not indicated, the differential CNV analysis is skipped.
* ``--diffcnv-method``: (optional) the method used to perform the differential CNV analysis. Can be Mann-Whitney U (mwu, default) or a t-test (ttest).
* ``--diffcnv-correction``: if you want to obtain False Discovery Rate (FDR) corrected results, add this flag. It is recommended to use these results rather than uncorrected p-values, as these can result in numerous false discoveries when blindly testing for differential expression (for more information, read https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1716-1)

Example usage to compute the differential CNV analysis with default values (Mann Whitney U test, no FDR correction)

.. code-block:: bash

   $ python -m cansig.run.pipeline data/pipeline-tutorial/data.hdf5 
                                   --batch batch
                                   --gene-sets MSigDB_Hallmark_2020 \
                                   --dimensions 4 6 --model-runs 1 \
                                   --clusters 2 3 5  --cluster-runs 1 \
                                   --output tutorial-output \
                                   --diffcnv

This will result in the following file being added to the ``tutorial-output/`` directory, in addition to all described in the :ref:`index`
* ``postprocessing/``:
   * ``{rundir}/differential-cnvs.csv``: file containing the columns for each cluster cl
      - {cl}\_pvalues: contains the p values of the test cl vs rest
      - {cl}\_perc\_{gains/losses}: contains the percentage of cells in the cluster showing a
            gain/loss at this region
      - {cl}\_rest\_{gains/losses}: contains the percentage of cells in all but the cluster showing a
            gain/loss at this region

Example usage to compute the differential CNV analysis with a t-test and with FDR corrected values (ie q-values)

.. code-block:: bash

   $ python -m cansig.run.pipeline data/pipeline-tutorial/data.hdf5 
                                   --batch batch
                                   --gene-sets MSigDB_Hallmark_2020 \
                                   --dimensions 4 6 --model-runs 1 \
                                   --clusters 2 3 5  --cluster-runs 1 \
                                   --output tutorial-output \
                                   --diffcnv \
                                   --diffcnv-method ttest \
                                   --diffcnv-correction

This will result in the same file as in the previous example with the addition of the columns
      - "{cl}\_qvalues": contains the q values of the test cl vs rest

.. note::
   If trying to run this function as such on a data object that has not been processed with our preprocessing module, this will result in an ValueError

Differential CNV analysis for data not processed with our module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This calling assumes the data was not processed using our module. In this case, you must provide a path to a .csv file that contains pre-called CNV. 
This array must have the following structure:
- first column should contain the cell IDs. The cell IDs must correspond to the cell IDs in the data object provided.
- first row should contain the region IDs. This can correspond to any region you wish - if you have your own mapping, this could also be simply integers corresponding to specific regions.
- values in the cells must be (positive or negative) integers. We thus assume your data has been discretized - running on a CNV array with non integer values will result in spurious results.

.. note::
   We in the tutorial data, we provide the file ``cnv_array.csv`` as an example valid CNV array.

The analysis can be controlled through four arguments in the command line:

* ``--diffcnv``: this flag needs to be added for the differential CNV analysis to be performed. If not indicated, the differential CNV analysis is skipped.
* ``--diffcnv-method``: (optional) the method used to perform the differential CNV analysis. Can be Mann-Whitney U (mwu, default) or a t-test (ttest).
* ``--diffcnv-correction``: if you want to obtain False Discovery Rate (FDR) corrected results, add this flag. It is recommended to use these results rather than uncorrected p-values, as these can result in numerous false discoveries when blindly testing for differential expression (for more information, read https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1716-1)
* ``--cnvarray``: the path to the CNV array as previously described

.. note::
   Forgetting to add the ``--cnvarray`` flag will result in the differential CNV analysis being run on the data object provided, and thus will likely throw an error if this data has not been obtained using our preprocessing module.

Example usage to compute the differential CNV analysis with default values (Mann Whitney U test, no FDR correction)

.. code-block:: bash

   $ python -m cansig.run.pipeline data/pipeline-tutorial/data.hdf5 
                                   --batch batch
                                   --gene-sets MSigDB_Hallmark_2020 \
                                   --dimensions 4 6 --model-runs 1 \
                                   --clusters 2 3 5  --cluster-runs 1 \
                                   --output tutorial-output \
                                   --diffcnv \
                                   --cnvarray data/pipeline-tutorial/cnv_array.csv

This will result in the following file being added to the ``tutorial-output/`` directory, in addition to all described in the :ref:`index`
* ``postprocessing/``:
   * ``{rundir}/differential-cnvs.csv``: file containing the columns for each cluster cl
      - {cl}\_pvalues: contains the p values of the test cl vs rest
      - {cl}\_perc\_{gains/losses}: contains the percentage of cells in the cluster showing a
            gain/loss at this region
      - {cl}\_rest\_{gains/losses}: contains the percentage of cells in all but the cluster showing a
            gain/loss at this region

Example usage to compute the differential CNV analysis with a t-test and with FDR corrected values (ie q-values)

.. code-block:: bash

   $ python -m cansig.run.pipeline data/pipeline-tutorial/data.hdf5 
                                   --batch batch
                                   --gene-sets MSigDB_Hallmark_2020 \
                                   --dimensions 4 6 --model-runs 1 \
                                   --clusters 2 3 5  --cluster-runs 1 \
                                   --output tutorial-output \
                                   --diffcnv \
                                   --diffcnv-method ttest \
                                   --diffcnv-correction \
                                   --cnvarray data/pipeline-tutorial/cnv_array.csv

This will result in the same file as in the previous example with the addition of the columns
      - "{cl}\_qvalues": contains the q values of the test cl vs rest