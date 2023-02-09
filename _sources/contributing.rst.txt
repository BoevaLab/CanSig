.. _contribution-guide:

Contributing
============


Submitting an issue
-------------------

If you find a bug, please submit `a new issue <https://github.com/BoevaLab/CanSig/issues>`_.

To be able to reproduce a bug, we will usually need the following information:
  - Versions of Python packages used (in particular version of this library).
  - A minimal code snippet allowing us to reproduce the bug.
  - What is the desired behaviour in the reported case?
  - What is the actual behaviour?


Submitting a Pull Request
-------------------------
Do:
  - Do consider submitting a draft pull request with a description of proposed changes, before starting actual work.
  - Do follow `PEP8 <https://peps.python.org/pep-0008/>`_ and `Google Style Guide <https://google.github.io/styleguide/pyguide.html>`_.
  - Do write unit tests (``tests`` directory). We use `pytest <https://docs.pytest.org>`_.
  - Do write `Google-style docstrings <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html>`_.
  - Do write high-level documentation (e.g., examples and tutorials, illustrating introduced features) in ``docs``. We use `Sphinx <https://www.sphinx-doc.org>`_  and `ReStructuredText <https://docutils.sourceforge.io/rst.html>`_. For a reference, see e.g., `this cheatsheet <https://sphinx-tutorial.readthedocs.io/cheatsheet/>`_.
  - Do check the Development section.

Don't:
  - Don't include license information. This project is BSD-3 licensed and by submitting your pull request you implicitly and irrevocably agree to use this.
  - Don't implement too many ideas in a single pull request. Multiple features should be implemented in separate pull requests.



For new developers
------------------

Getting started with the code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Clone the repository:

.. code-block:: bash

   $ git clone git@github.com:BoevaLab/CanSig.git
   $ cd CanSig

Install the module in editable mode with the dependencies needed to run the test suite:

.. code-block:: bash

   $ pip install -e ".[test]"

At this point you should be able to run the unit test suite:

.. code-block:: bash

   $ pytest


Tooling
^^^^^^^

As a developer, you will use some tools increasing code quality as well as packages used to build the documentation. You can install them by running

.. code-block:: bash

   $ pip install -r requirements-dev.txt
   $ pre-commit install

We use `black <https://github.com/psf/black>`_ and `flake8 <https://flake8.pycqa.org/en/latest/>`_ to lint the code, `pytype <https://github.com/google/pytype>`_ to check whether the types agree, and `pytest <https://docs.pytest.org>`_ to unit test the code.
These code quality checks are applied to every Pull Request via a workflow in the ``.github/workflows`` directory.


Working with documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

We write `Google-style docstrings <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html>`_ and high-level documentation (e.g., examples and tutorials, illustrating introduced features) in ``docs``.

Although docstrings (in Python code) are assembled automatically, high-level documentation (including tutorials), needs to be written manually. For this, we use `Sphinx <https://www.sphinx-doc.org>`_  and `ReStructuredText <https://docutils.sourceforge.io/rst.html>`_. For a reference, see e.g., `this cheatsheet <https://sphinx-tutorial.readthedocs.io/cheatsheet/>`_.

Whenever a Pull Request is submitted (but not merged yet), we do a test build.
Whenever a Pull Request is merged into the ``main`` branch, a new version of the documentation is automatically built and deployed. Both workflows are in ``.github/workflows``.
(Note that if you want to change the Sphinx theme you will need to modify both workflows as well as the requirements files).

We store the source code of the documentation in the ``docs`` directory. To build it locally, use

.. code-block:: bash

   $ cd docs
   $ sphinx-apidoc -o source/ ../src/cansig
   $ make html

In the ``_build`` directory, you should see ``index.html`` file. Open it with a web-browser of your choice.


Requirements
^^^^^^^^^^^^

Python packages needed to run the CanSig code should be specified inside ``setup.cfg`` (it also includes optional testing dependencies).
Developer tools (described above) are specified in ``requirements-dev.txt``.
We also keep the file ``requirements.txt`` which is supposed to be the union of all possible dependencies.
Our CI pipeline uses this file for caching the dependencies (so that it runs faster).
We check whether the ``requirements.txt`` is the union of the other two files using a simple script ``.ci/check-requirements-cache.py``.

.. note:: Currently this script *does not* parse lines containing both the module and the comment. It is very simple.


Releasing a new version
-----------------------
In this section we briefly discuss how to upload new version to PyPI. It is based on the `official instructions <https://packaging.python.org/en/latest/tutorials/packaging-projects/>`_.

#. Activate the development environment.
#. Install ``twine`` (if you do not have it already) with ``pip install --upgrade twine``.
#. Bump version in the ``setup.cfg`` file. We use `semantic versioning <https://py-pkgs.org/07-releasing-versioning.html>`_.
#. Make sure that the ``dist/`` directory is empty.
#. Run ``python -m build`` and check if the ``dist/`` directory contains only the new version.
#. Upload the package to TestPyPI with ``twine upload --repository testpypi dist/*``.
#. Go to a different directory (which does not contain ``CanSig``) and using a *new environment*, install the package from TestPyPI: ``pip install -i https://test.pypi.org/simple/ cansig``. It is often convenient to reinstall the version from the previous test with ``pip uninstall cansig``.
#. Test whether the package installed properly. For example, run ``python -m cansig.run.pipeline --help`` and check if it does not raise exceptions.
#. If the package seems to work, go back to the development environment (with ``twine`` installed) and upload it to the PyPI via ``twine upload dist/*``.
#. Commit the new ``setup.cfg`` file and submit a PR which bumps the version.

