.. _troubleshooting:

Troubleshooting / Checklist
===========================
  
On this page you will find a checklist to go through before running CanSig to ensure it will run as smoothly as possible, as well as a troubleshooting section for most common encountered problems. 
If you are confident you tick all the boxes for the checklist and cannot find your issue mentioned in the troubleshooting, please open an issue so we can adress it!

Checklist
^^^^^^^^^

General checklist 

Preprocessing checklist 

Differential CNV checklist 

Troubleshooting
^^^^^^^^^^^^^^^

GSEA will not run 
- Is your gene annotation consistent? 
    If your genes are annotated using Entrez in your original expression object, then you should run GSEA using an adapted .gmt file. You can download files for official gene ID and for Entrez on the MSigDB website (https://www.gsea-msigdb.org/gsea/downloads.jsp#msigdb).
    If your genes are annotated using another system, you should first translate into Official Gene ID or Entrez. 
    

.. todo::

   The tutorial to use CanSig from a python script/jupyter notebook.
   It will contain a link to a Colab tutorial