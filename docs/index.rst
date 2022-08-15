.. beams2D documentation master file

Welcome to beams2D's documentation!
===================================

The package is a framework to simulate the dynamics of structurs
using beam elements (e.g. Euler-Bernoulli Beam) in 2 dimensions.

.. autosummary::
   :toctree: generated/
   :template: custom_module.rst

   beams2D


Modules Reference
-----------------

A sub-package (module) :mod:`src` gathers main sources for geometry, mesh and
assembly of the system classes

.. autosummary::
   :toctree: generated/
   :recursive:

   beams2D.src

The src contains the main modules/classes

.. autosummary::
   :toctree: generated/
   :recursive:

   beams2D.src.geom
   beams2D.src.mesh
   beams2D.src.elements
   beams2D.src.assembleSystem
   beams2D.src.analysis


Module :mod:`inputs` is used to read in the 'json' input file.

.. autosummary::
   :recursive:
   :toctree: generated/

   beams2D.inputs

Module :mod:`outputs` gathers output/print commands.

.. autosummary::
   :recursive:
   :toctree: generated/

   beams2D.outputs

Module :mod:`genHelpers` gathers useful functions

.. autosummary::
   :recursive:
   :toctree: generated/

   beams2D.genHelpers

Examples
========

See tutorials directory for detailed inputs/outputs

.. toctree::
    :maxdepth: 2
    :caption: Notebooks

    examples/threeBeams.ipynb

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
