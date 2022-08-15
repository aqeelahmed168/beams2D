import beams2D as b2d
import numpy as np
from beams2D import PrettyPrint as pp

inputsB2D = b2d.Inputs(infile="inputs.json")

geom = b2d.Geom(inputsB2D.geom)
geom.process()
# geom.writeVTK()

mesh = b2d.Mesh(inputsB2D.mesh)
mesh.process(geom, printMeshInfo=0)

#%%
# assemble the system using a given element type
sysB2D = b2d.Assemble(mesh)
# build the global stiffness matrix 
sysB2D.buildGlobalMK(inputsB2D.strucParams, elemType="EBB")
