.. _formatting:

Input data format
=================

We highly recommend you use our preprocessing module (see :ref:`preprocessing`) to preprocess your data. 
This will ensure the data you input in CanSig is in the correct format and will enable smooth use of functionalities like differential CNV analysis.

However, we support running the pipeline without using our module. 
CanSig takes an .h5ad object as input. This tutorial will walk you through obtaining this format using more a more widespread data format (.csv).

You must have: 

* expression.csv, matrix containing the raw counts (eg UMI), with cell IDs as first column and gene IDs as first row.
* observation.csv, matrix containing information about the cells. Only mandatory information is the batch ID, but this matrix can contain known signatures pre-scored, cell type annotation, number of counts in the cell, etc.

.. note::
    If you have a .txt file, this function should also work, if your file is organized in the same manner as a .csv file, meaning that there is a consistent separator (be it tab or comma) that separates your data into columns.

If you are using CanSig as a model, you will need data for malignant and for non-malignant cells. You should thus have two expression matrices (eg expression-malignant.csv and expression-non-malignant.csv) and two associated observation matrices. You should then run the following procedure both for the malignant and non-malignant matrices.
If you are using scVI as a model, you will only need data for malignant cells. 

.. note::
    You can use any type of gene annotation for your genes as long as your are consistent throughout the pipeline. In particular, GSEA will only work with either Official Gene ID or Entrez. If you have ENSEMBL annotations, we recommend you map them to Official Gene ID before. 

Here is an example of what the beginning of the expression matrix should look like

.. code-block:: python

    	        EEA1  CCDC107  ZNF442  THOC2    GPX3   ASPHD2   GINS1
    cell_1      0.0     1.0	0.0	0.0	3.0	0.0	1.0
    cell_2	0.0	0.0	0.0	0.0	0.0	0.0	0.0
    cell_3	0.0	4.0	0.0	0.0	2.0	0.0	0.0
    cell_4	0.0	0.0	0.0	3.0	0.0	0.0	0.0
    cell_5	0.0	2.0	0.0	0.0	0.0	0.0	0.0
    cell_6	1.0	1.0	0.0	0.0	1.0	0.0	0.0
    cell_7	1.0	1.0	0.0	2.0	0.0	0.0	0.0
    cell_8	2.0	1.0	0.0	0.0	0.0	0.0	0.0

Here is an example of what the beginning of the observation matrix should look like

.. code-block:: python

            batch	n_counts   n_genes  pct_counts_mt
    cell_1  P28T        39525.0    4551     0.2504744
    cell_2  P28T	18224.0    3826	    0.548727
    cell_3  P30T	10079.0    2685     0.00992162

.. note::
    The only required column in the observation matrix of malignant cells is batch. However, if you are also using a non-malignant dataset, there should also be a celltype column.

If your data is formatted as previously explained, then you can simply go into the directory that contains your data and run 

.. code-block:: 

    $ python -m cansig.run.format --expression expression.csv \
                                            --observation observation.csv \
                                            --filename full_data

This will save a file called "full_data.h5ad" in the current directory. 

.. note::
    Make sure you either call the function from the directory where the data is stored or provide the full path to the data in the above call. 
    In which case, the call would look something like this

    .. code-block:: 

        $ python -m cansig.run.format --expression path/to/expression.csv \
                                                --observation path/to/observation.csv \
                                                --filename full_data
.. note::
    If you are calling this function twice, don't forget to change the ``--filename`` argument!

 
