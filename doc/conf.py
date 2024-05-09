# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

#sys.path.extend([os.path.abspath('../python'),os.path.abspath('../lib'),os.path.abspath('../core/pydis/python')])



# -- Project information -----------------------------------------------------

project = 'OpenDiS'
copyright = '2024, MicroNano'
author = 'MicroNano'

from datetime import datetime
year = datetime.now().year
if year > 2024:
    copyright = f'2024-{year}, {author}'
else:
    copyright = f'2024, {author}'

# The full version, including alpha/beta/rc tags
release = 'x.x.x'


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['myst_parser']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']
