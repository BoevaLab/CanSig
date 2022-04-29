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

