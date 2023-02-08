.. _pipeline-advanced:

Advanced pipeline usage
=======================

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

Running differential CNV analysis
---------------------------------

You have the option to perform differential CNV analysis. With original CNV calls, this will output differential CNV regions between each cluster and the rest, and information about the percentage of gains/losses in the cluster and in the rest, and the number of patients showing a gain/loss in the region.
This module is deactivated by default. There are two main ways to run this analysis: the first assumes that you are using a data object that has been obtained using our preprocessing module (see :ref:`preprocessing`), the second can be run if provided with a external discretized CNV calling, even if the data object has not been obtained through our preprocessing module.

.. note::
   
   The data provided for the tutorial has been processed using our preprocessing module, and can be thus used for differential CNV analysis assuming so.

Differential CNV analysis for data preprocessed with our module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This calling assumes the data was preprocessing using our preprocessing module (see :ref:`preprocessing`).
This means that the data object you provide will contain the following:

   - "X_cnv" in data.obsm: the CNV called using our preprocessing module
   - "chromosome" in data.var.columns: the chromosome to which the gene belongs
   - "cnv_called" in data.var.columns: if this gene was used for the infercnv call (see ``cansig._preprocessing`` for more details on the CNV calling procedure)
   - "start" in data.var.columns: the start position of the gene on the chromosome
   - "cnv" in data.uns: a summary of the infercnv run
   - "chr_pos" in data.uns["cnv"]: a dictionary containing the mapping between the chromosome and the index of the regions in the cnv array

The analysis can be controlled through three arguments in the command line:

* ``--diffcnv``: this flag needs to be added for the differential CNV analysis to be performed. If not indicated, the differential CNV analysis is skipped.
* ``--subclonalcnv``: (optional) when added, performs the differential CNV analysis using CNV smoothed by subclone rather than on a cell level. This means the CNV of a cell will be that of the subclone it belongs to. This type of call is less noisy but might hide smaller CNV regions or smaller subclone populations that might not have been found with infercnv.
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

This will result in the following file being added to the ``tutorial-output/`` directory, in addition to all the files/directories described on the homepage.

* ``postprocessing/``:
   * ``{rundir}/differential-cnvs.csv``: file containing the columns for each cluster cl
      - {cl}\_pvalues: contains the p values of the test cl vs rest
      - {cl}\_perc\_{gains/losses}: contains the percentage of cells in the cluster showing a gain/loss at this region
      - {cl}\_rest\_{gains/losses}: contains the percentage of cells in all but the cluster showing a gain/loss at this region
      - {cl}\_patients\_{gain/loss}: contains the number of patients that show a gain/loss in this region in this cluster. Specifically, we count a patient as showing a gain/loss in the region if at least one cell in the cluster belongs to this patient and shows a gain/loss.

.. note::
   We use the batch ID as a proxy for the patient in the computation of the number of patients showing a gain/loss. If there are several patients in one batch or several batches per patient, this will count the number of batches showing a gain/loss, not the number of patients.

Example usage to compute the differential CNV analysis with a t-test, smoothing on a subclonal level, and with FDR corrected values (ie q-values)

.. code-block:: bash

   $ python -m cansig.run.pipeline data/pipeline-tutorial/data.hdf5 
                                   --batch batch
                                   --gene-sets MSigDB_Hallmark_2020 \
                                   --dimensions 4 6 --model-runs 1 \
                                   --clusters 2 3 5  --cluster-runs 1 \
                                   --output tutorial-output \
                                   --diffcnv \
                                   --subclonalcnv \
                                   --diffcnv-method ttest \
                                   --diffcnv-correction

This will result in the same file as in the previous example with the addition of the columns
      - "{cl}\_qvalues": contains the q values of the test cl vs rest

.. note::
   If trying to run this function as such on a data object that has not been processed with our preprocessing module, this will result in an ValueError

.. 
  .. Differential CNV analysis for data not processed with our module
  .. ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  This calling assumes the data was not processed using our module. In this case, you must provide a path to a .csv file that contains pre-called CNV. 
  This array must have the following structure:
  - first column should contain the cell IDs. The cell IDs must correspond to the cell IDs in the data object provided.
  - first row should contain the region IDs. This can correspond to any region you wish - if you have your own mapping, this could also be simply integers corresponding to specific regions.
  - values in the cells must be (positive or negative) integers. We thus assume your data has been discretized - running on a CNV array with non integer values will result in spurious results.

..
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

  This will result in the following file being added to the ``tutorial-output/`` directory, in addition to all the files/directories described on the homepage.

  * ``postprocessing/``:
   * ``{rundir}/differential-cnvs.csv``: file containing the columns for each cluster cl
      - {cl}\_pvalues: contains the p values of the test cl vs rest
      - {cl}\_perc\_{gains/losses}: contains the percentage of cells in the cluster showing a gain/loss at this region
      - {cl}\_rest\_{gains/losses}: contains the percentage of cells in all but the cluster showing a gain/loss at this region

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

To help creating such directories, we created auxiliary functions:

.. code-block:: python

   import cansig.filesystem as fs
   
   output_dir = fs.IntegrationDir(output, create=True)
   fs.save_settings(settings=config, path=output_dir.integration_settings)
   fs.save_latent_representations(representations=representations, path=output_dir.latent_representations)


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

Understanding all the pipeline command line options
---------------------------------------------------

The core of the CanSig tool is the command:

.. code-block::

   $ python -m cansig.run.pipeline

which runs the entire analysis pipeline.

This command comes with numerous flags to enable you to control the inputs/outputs of CanSig.
In this section we describe all available flags.

.. note::
   You can get also get information on these flags by running 
   
   .. code-block::

      $ python -m cansig.run.pipeline --help


* ``--batch``: the name of the column in which the batch information is stored. 
  This will typically be the name of the column where the sample ID is stored, as generally each sample is processed separately.
* ``--continuous-covariates``: continuous covariates for which one wishes to condition on in the integration model. This could be the cell cycle score or the percentage of mitochondrial counts, for example.
* ``--discrete-covariates``: similar for ``--continuos-covariates``, but should be used for columns with discrete/categorical variables. This could be the subclonal structure for example.
* ``--gene-sets``: gene set to use for GSEA. The input should be a string (valid for Enrichr) or a ``.gmt`` file. More information on these sets can found on the MSigDB website. 
* ``--n-top-genes N_TOP_GENES``: number of the most highly variable genes to use.
* ``--model-runs``: number of random seeds used for initialization of each integration model. If you are running a integration model with 4 latent dimensions and 3 model runs, this will result in 3 different latent representations for the same number of dimensions.
* ``--cluster-runs``: number of random seeds used for initizalization of each postprocessing run. If you are running postprocessing with 6 clusters and 2 cluster runs, this will result in 2 different clustering partitions for the same number of clusters.
* ``--max-epochs``: maximum number of epochs the integration model will be trained for. 
* ``--dimensions``: list of number of latent dimensions used for integration 
* ``--clusters``: list of number of clusters used for postprocessing
* ``--output``: name of the folder in which the output will be stored (the folder will be created).
* ``--dim-reduction``: the name of the dimensionality reduction method used to plot the latent space (can be PCA, UMAP or both).
* ``--sigcols``: the name of the columns in the ``.obs`` according to which to color the plots 
* ``--disable-plots``: if set, no plots will be created to visualize the latent space (the run will be quicker when this option is on)
* ``--dgex-method``: method used to perform the differential gene expression analysis.
* ``--ngenessig``: number of genes to use to define a signature to score 
* ``--corrmethod``: correlation method used to correlate de novo found signatures
* ``--disable-signatures``: if set, no information linked to de novo signatures found will be saved (the run will require less memory when this option is on). 
* ``--diffcnv``: if set, the differential CNV analysis will be run. For more information, see the part about running the differential CNV analysis on this page.
* ``--diffcnv-method``: the method used to perform the differential CNV.
* ``--subclonalcnv``: if set, the differential CNV analysis will be run using the subclonal inferred CNV representation for each cell, rather than the per cell CNV call.
* ``--diffcnv-correction``: if set, the False Discovery Rate corrected q-value will be computed for the differential CNV analysis.
* ``--cnvarray``: if running the differential CNV analysis on an external array (for those who did not preprocess their data using our preprocessing module), the path to the CNV array used for differential CNV. Using this flag will automatically *disable* running the differential CNV on the AnnData object.
* ``--save-intermediate``: whether the intermediate results should be saved. By default, the results are saved. Turning this off is discouraged, unless the system memory is very limited.
* ``--linkage``: the linkage used for the agglomerative clustering for metasignatures.
* ``--sim-method``: similarity metric to be used to compare signatures.
* ``--threshold``: the threshold above which a metasignature is considered too correlated with another.
* ``--pat-specific-threshold``: the threshold above which a metasignature is considered patient-specific.
* ``--model``, ``--n-latent-batch-effect`` and ``--n-latent-cnv``: settings for the new experimental integration method. Currently we suggest to use the default integration method in the CanSig pipeline.

