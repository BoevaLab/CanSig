.. CanSig documentation master file, created by
   sphinx-quickstart on Wed Apr 20 13:39:48 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CanSig!
==================



.. toctree::
   :maxdepth: 2
   :caption: Contents:


Contributing
============

Install the module in editable mode and with the dependencies needed to run the test suite:

.. code-block:: bash

   $ pip install -e ".[test]"

At this point you should be able to run the unit test suite:

.. code-block:: bash

   $ pytest


Developer tools
---------------

As a developer, you will however need more tools. You can install them by running

.. code-block:: bash

   $ pip install -r requirements-dev.txt


.. todo::

   Add the description of code quality tools and pre-commit.
   

Building the documentation
--------------------------

We use `Sphinx <https://www.sphinx-doc.org>`_  and `ReStructuredText <https://docutils.sourceforge.io/rst.html>`_. For a reference, see e.g., `this cheatsheet <https://sphinx-tutorial.readthedocs.io/cheatsheet/>`_.

Whenever a Pull Request is merged into the ``main`` branch, a new version of the documentation is automatically built and deployed.
Although docstrings (in Python code) are assembled automatically, high-level documentation (including tutorials), needs to be written manually.
We store the source code of the documentation in the ``docs`` directory. To build it locally, use

.. code-block:: bash

   $ cd docs
   $ make html

In the ``_build`` directory, you should see ``index.html`` file. Open it with a webbrowser of your choice.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

