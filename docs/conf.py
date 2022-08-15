# Configuration file for the Sphinx documentation builder.


autodoc_mock_imports = ["pandas", "numpy", "scipy"]


# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------

project = 'beams2D'
copyright = '2022, Aqeel Ahmed'
author = 'Aqeel Ahmed'

# The full version, including alpha/beta/rc tags
release = '0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc", 
    "sphinx.ext.coverage", 
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary", 
    # 'autoapi.extension',
    "sphinx.ext.doctest",
    "sphinx.ext.todo",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.imgmath",
    "nbsphinx",
    "nbsphinx_link"
]


# for readthedoc set to ghostscript
# tikz_proc_suite = 'GhostScript' 
# latex_elements = {
#      'packages': r'\usepackage{tikz}'
# }


nbsphinx_allow_errors = False
nbsphinx_execute = 'never'

source_suffix = [".rst", ".md", ".ipynb"]

imgmath_image_format = 'svg'
#imgmath_image_format = 'png'
# pngmath_latex="/Library/TeX/Distributions/Programs/texbin/latex"
# pngmath_dvipng="/Library/TeX/Distributions/Programs/texbin/"

autodoc_default_options = {
    'ignore-module-all': True
}

add_module_names = False
autoapi_type = 'python'

# autoapi_dirs = ['../beams2D']
autosummary_generate = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']


# -------- Options for HTML output --------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']
