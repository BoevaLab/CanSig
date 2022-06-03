.. _formatting:

Formatting your data to run CanSig
==================================

We highly recommend you use our preprocessing module (see :ref:`preprocessing`) to preprocess your data. 
This will ensure the data you input in CanSig is in the correct format and will enable smooth use of utilities like differential CNV analysis.

However, we support running the pipeline without using our module. 
CanSig takes an .h5ad object as input. This tutorial will walk you through obtaining this format using more a more widespread data format (.csv).

You must have: 

* a .csv matrix containing the raw counts (eg UMI), with cell IDs as first column and gene IDs as first row.
* a .csv matrix containing information about the cells. Only mandatory information is the batch ID, but this matrix can contain known signatures pre-scored, cell type annotation, number of counts in the cell, etc.

.. note::
    You can use any type of gene annotation for your genes as long as your are consistent throughout the pipeline. This means 