"""Sub-package combining source files relating to geometry, mesh and elements

.. autosummary::
   :toctree:
   
   geom
   mesh
   elements
   assembleSystem
   analysis
"""

from .geom.geom import Geom
from .mesh.mesh import Mesh
from .elements.elementFactory import ElementFactory
from .assembleSystem import Assemble
from .analysis import Modal


__all__ = ["Geom", "Mesh", "ElementFactory", "Assemble", "Modal"]