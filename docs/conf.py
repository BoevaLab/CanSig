# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
#
# We also used https://samnicholls.net/2016/06/15/how-to-sphinx-readthedocs/
#
#
# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath("../src"))


# -- Project information -----------------------------------------------------

project = "CanSig"
copyright = "2022 CanSig team"
author = "CanSig contributors"

# The full version, including alpha/beta/rc tags
release = ""
latex_elements = {"releasename": ""}


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.napoleon",
    "sphinx.ext.githubpages",
    "sphinx.ext.todo",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Extensions --------------------------------------------------------------

# sphinx.ext.todo
# (Support for Todo comments).
todo_include_todos = True

# sphinx.ext.napoleon
# (Support for Google-style docstrings).
napoleon_google_docstring = True
napoleon_use_param = False
napoleon_use_ivar = True

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
# See https://sphinx-themes.org/ for other themes.
# Note that if you change the theme, you may need to modify the GitHub
# Action responsible for deploying the webpage.
html_theme = "pydata_sphinx_theme"


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# Options for the theme, as specified there:
# https://pydata-sphinx-theme.readthedocs.io/en/stable/

html_theme_options = {
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/BoevaLab/CanSig",
            "icon": "fab fa-github-square",
            "type": "fontawesome",
        }
    ]
}
