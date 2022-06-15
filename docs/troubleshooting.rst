.. _troubleshooting:

Troubleshooting / Checklist
===========================
  
On this page you will find a checklist to go through before running CanSig to ensure it will run as smoothly as possible, as well as a troubleshooting section for most common encountered problems. 
If you are confident you tick all the boxes for the checklist and cannot find your issue mentioned in the troubleshooting, please open an issue so we can adress it!

Checklist
^^^^^^^^^

General checklist 

.. todo:: 

  Add in the general checklist 

Preprocessing checklist 

.. todo:: 

  Add in the preprocessing checklist 

Differential CNV checklist 

.. todo:: 

  Add in the differential CNV checklist

Troubleshooting
^^^^^^^^^^^^^^^

1. GSEA will not run 

    * Is your gene annotation consistent? 
      If your genes are annotated using Entrez in your original expression object, then you should run GSEA using an adapted .gmt file. You can download files for official gene ID and for Entrez on the MSigDB website (https://www.gsea-msigdb.org/gsea/downloads.jsp#msigdb).
      If your genes are annotated using another system, you should first translate into Official Gene ID or Entrez.


2. The differential CNV analysis doesn't work 

    * Did you run our preprocessing module? If not you need to run the differential CNV analysis in a slightly different way, see :ref:`pipeline-advanced`.
    * If you get a message error saying the mapping shapes do not coincide, it is possible the window size you used for preprocessing is too big. When the window size is bigger than the number of regions for a specific chromosome, infercnvpy returns an object the size of the window, thus losing mappability to specific chromosomal regions (see Issue open on their github). Try decreasing the window size.
