# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'xray_stack_models'
copyright = '2023, Johan Comparat'
author = 'Johan Comparat'
release = '1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.coverage',
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
]
templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

import sys
import os
#import shlex

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

# root package
sys.path.insert(0, os.path.abspath(os.path.join(os.environ['GIT_STMOD'], 'src')) )
sys.path.insert(0, os.path.abspath(os.path.join(os.environ['GIT_STMOD'], 'src', 'models')) )
sys.path.insert(0, os.path.abspath(os.path.join(os.environ['GIT_STMOD'], 'src', 'Mpipelines')) )

# galaxy package
#sys.path.insert(0, os.path.abspath(os.environ['GIT_STMOD']))
#sys.path.insert(0, os.path.abspath(os.path.join(os.environ['GIT_STMOD'],'python')))

# stellar population model
#sys.path.insert(0, os.path.abspath(os.environ['GIT_STMOD']))
#sys.path.insert(0, os.path.abspath(os.path.join(os.environ['GIT_STMOD'],'python')))
