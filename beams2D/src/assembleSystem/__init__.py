"""
Assemble Module.

Assemble global Degree of Freedoms from the mesh data, build M, C, and K 
using specfic elements

"""

from .assemble import Assemble

__all__ = ["Assemble"]
