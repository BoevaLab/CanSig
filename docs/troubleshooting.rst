.. _troubleshooting:

Troubleshooting
===============

GSEA does not run 
^^^^^^^^^^^^^^^^^

* Is your gene annotation consistent? 
* If your genes are annotated using Entrez in your original expression object, then you should run GSEA using an adapted .gmt file. You can download files for official gene ID and for Entrez on the `MSigDB website <https://www.gsea-msigdb.org/gsea/downloads.jsp#msigdb>`_.
* If your genes are annotated using another system, you should first translate into Official Gene ID or Entrez.

Differential CNV analysis does not work
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Did you run our preprocessing module? If not you need to run the differential CNV analysis in a slightly different way, see :ref:`pipeline-advanced`.
* If you get a message error saying the mapping shapes do not coincide, it is possible the window size you used for preprocessing is too big. When the window size is bigger than the number of regions for a specific chromosome, InferCNVPy returns an object the size of the window, thus losing mappability to specific chromosomal regions (see `Issue 37 <https://github.com/icbi-lab/infercnvpy/issues/37>`_ open in InferCNVPy repository). Try decreasing the window size.
