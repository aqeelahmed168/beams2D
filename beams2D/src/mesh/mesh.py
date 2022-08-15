"""Mesh class
==============

.. autoclass:: Mesh
   :members:
   :private-members:

   .. automethod:: __init__
"""

import numpy as np
import pandas as pd
import os

# for type hint (instant of class Geom)
from ..geom.geom import Geom
from ...inputs.inputs import Inputs

class Mesh:
    """Class to build mesh from geom and inputs"""
    def __init__(self, msdf: Inputs):
        """Init method

        Args:
            msdf (Inputs): attribute mesh of the object Inputs
        Returns:
            None
        """
        self.msdf = msdf
    
    def process(self, geom1: Geom, printDebugMessage=0, printMeshInfo=1):
        """Process mesh based on the geom inputs

        Args:
            geom1 (Geom): object of type Geom
            printDebugMessage: flag to show additional outputs

        Returns:
            None
        """
        msdf = self.msdf

        self.nStrucs = geom1.nStrucs

        mInputs = len(msdf.index)
        # check the number of nStrucs is equal to geom
        # set initially to 0, nan is for floats only
        strucElems = np.full(geom1.nStrucs, 0, dtype=np.integer)

        # put all print info in the meshInfo
        meshInfo = "Mesh info:\n"

        # loop over the list of the elements to set the elements
        for i in range(mInputs):
            iStruc = msdf.structure[i]
            strucElems[iStruc] = msdf.nNodes[i]-1

        # set all the remaning 0 elems tstrucs to 1
        strucElems[np.nonzero(strucElems==0)] = 1

        self.strucElems = strucElems

        meshInfo += f"strucElems (number of elements/structure)\n"
        meshInfo += f"{strucElems}\n"
        # print(meshInfo)

        # map of elements and the structure
        maxElmsStruc = max(strucElems)
        nStrucs = geom1.nStrucs
        # 1st row is number of elems and 2nd to last is elem idx if not -9999
        strucElemConn = np.full((nStrucs,maxElmsStruc+1),-9999,dtype=np.integer)
        # assign elems to the strucs
        idxElem = -1
        for i in range(nStrucs):
            nElems = strucElems[i]
            strucElemConn[i,0] = nElems

            for j in range(nElems):
                idxElem = idxElem+1
                strucElemConn[i,j+1] = idxElem

        meshInfo += f"strucElemConn (structure to element connectivity)\n"
        meshInfo += f"{strucElemConn}\n"

        self.strucElemConn = strucElemConn
        self.buildElemNodes(geom1, printDebugMessage=printDebugMessage)

        meshInfo += f"elemNodes (mesh element nodes)\n"
        meshInfo += f"{self.elemNodes}\n"

        meshInfo += f"nodeCoods (mesh node coordinates)\n"
        meshInfo += f"{self.nodeCoords}\n"

        if printMeshInfo:
            print(meshInfo)


    def buildElemNodes(self,geom1, printDebugMessage=0):
        #%% build element nodes
        nStrucs = self.nStrucs
        strucElemConn = self.strucElemConn
        msdf = self.msdf

        strucElems = self.strucElems
        nElems = np.sum(strucElems)
        # the index of nodes for each element
        elemNodes = np.full((nElems,2),-9999,dtype=int)
        maxNodes = nElems*2
        # the node coords for the mesh
        nodeCoords = np.full((maxNodes,2),np.nan)
        # loop over structures
        # then loop over elements (use start/end of struc to set elem nodes)
        # increment global node idx (associate with each elem/struc)
        # if 1 elem -> then just the end start nodes of geom

        debugMessage = 'Building mesh\n!'

        idxNode = -1
        for i in range(nStrucs):
            # get the elements for the each struct
            nElems = strucElemConn[i,0]
            if nElems > 1:
                # operate on 1st element and last element 1st
                iElem = strucElemConn[i,1]
                # print("iStruc %d iElem %d" %(i,iElem))
                debugMessage += f'iStruc {i} iElem {iElem}'
                # set the 1st node of the 1st element
                if elemNodes[iElem,0] == -9999:
                    # the increase node index
                    idxNode = idxNode+1
                    # set the new node to the element
                    elemNodes[iElem,0] = idxNode
                    geomNode = geom1.strucNodes[i,0]
                    # print("\tgeomNode", geomNode)
                    debugMessage += f'\tgeomNode {geomNode}'
                    nodeCoords[idxNode,:] = geom1.nodeCoords[geomNode]
                    # check the connected strcutures to this node
                    idxConnStruc = np.argwhere(geom1.strucNodes[:,0]==geomNode)
                    # print("\tidxConnStruc startNode", idxConnStruc)
                    debugMessage += f'\tidxConnStruc startNode {idxConnStruc}'
                    self.processNodeConnElems(elemNodes,\
                        idxConnStruc,i,iElem,iNode=0,nodeOther="firstNode")
                    # find the strucs with same geom node as end of another elem
                    idxConnStruc = np.argwhere(geom1.strucNodes[:,1]==geomNode)
                    # print("\tidxConnStruc endNode", idxConnStruc)
                    debugMessage += f'\tidxConnStruc endNode {idxConnStruc}'

                    # use iNode of 0 if the starting node is 1st node of the struc
                    # node lastNode means last node of another structure
                    self.processNodeConnElems(elemNodes,\
                        idxConnStruc,i,iElem,iNode=0,nodeOther="lastNode")

                # operate on inner nodes
                totStrucNodes = nElems+1
                # check if node coords for a given struc are in the input.json
                # it should be for nElems > 1
                # find the index of the structure in the input
                coordsStrucIdx = msdf.structure.index[msdf['structure']==i].to_numpy()
                if not len(coordsStrucIdx):
                    raise ValueError('check mesh inputs')
                coordsStrucIdx = coordsStrucIdx.item() # get the scalar out

                # check iStruc has some inputs for the node coords
                if len(msdf.nodeCoords[coordsStrucIdx]) > 0:
                    # print("Using input node coordinates for structure %d" %i)
                    # print(msdf.nodeCoords[coordsStrucIdx])
                    debugMessage += f'Using input node coordinates for structure {i}'
                    debugMessage += f'{msdf.nodeCoords[coordsStrucIdx]}'

                    strucNodeCoords = np.array(msdf.nodeCoords[coordsStrucIdx])
                    # check 1st and last node match with the structure
                    xlocInNodes = strucNodeCoords[:,0]
                    ylocInNodes = strucNodeCoords[:,1]
                else:
                    geomFirstNode = geom1.strucNodes[i,0]
                    geomLastNode = geom1.strucNodes[i,1]
                    geomFNodeCoords = geom1.nodeCoords[geomFirstNode]
                    geomLNodeCoords = geom1.nodeCoords[geomLastNode]
                    # if not given the node postion, linspace the interior nodes
                    xlocInNodes = np.linspace(geomFNodeCoords[0],\
                        geomLNodeCoords[0],num=totStrucNodes)
                    ylocInNodes = np.linspace(geomFNodeCoords[1],\
                        geomLNodeCoords[1],num=totStrucNodes)

                # skip 1st node of 1st elem and the last node of last elem
                for iNode in range(1,totStrucNodes-1):
                    idxNode = idxNode+1
                    # assign the nodes to the elems
                    # elem1 2nd node, elem2 1st node
                    elemNodes[iElem+iNode-1,1] = idxNode
                    elemNodes[iElem+iNode,0] = idxNode
                    # set the node coords
                    nodeCoords[idxNode] = [xlocInNodes[iNode],ylocInNodes[iNode]]

                # operate on last node of the element and last element 1st
                lastElemIdx = strucElemConn[i,0] # 0 is also number of elems
                iElem = strucElemConn[i,lastElemIdx]
                # print("iStruc %d iElem %d" %(i,iElem))
                debugMessage += f'iStruc {i} iElem {iElem}'
                # set the last node of the last element
                if elemNodes[iElem,1] == -9999:
                    # the increase node index
                    idxNode = idxNode+1
                    # set the new node to the element
                    elemNodes[iElem,1] = idxNode
                    geomNode = geom1.strucNodes[i,1]
                    # print("\tgeomNode", geomNode)
                    debugMessage += f'\tgeomNode {geomNode}'
                    nodeCoords[idxNode,:] = geom1.nodeCoords[geomNode]
                    # check the connected strcutures to this node
                    idxConnStruc = np.argwhere(geom1.strucNodes[:,0]==geomNode)
                    # print("\tidxConnStruc startNode", idxConnStruc)
                    debugMessage += f'\tidxConnStruc startNode {idxConnStruc}'
                    self.processNodeConnElems(elemNodes,\
                        idxConnStruc,i,iElem,iNode=1,nodeOther="firstNode")
                    # find the strucs with same geom node as end of another elem
                    idxConnStruc = np.argwhere(geom1.strucNodes[:,1]==geomNode)
                    # print("\tidxConnStruc endNode", idxConnStruc)
                    debugMessage += f'\tidxConnStruc endNode {idxConnStruc}'
                    # use iNode of 0 if the starting node is 1st node of the struc
                    # node lastNode means last node of another structure
                    self.processNodeConnElems(elemNodes,\
                        idxConnStruc,i,iElem,iNode=1,nodeOther="lastNode")


            # if only one element (elem mesh share same nodes as the geom)
            if (nElems == 1):
                iElem = strucElemConn[i,1]
                # print("iStruc %d iElem %d" %(i,iElem))
                debugMessage += f'iStruc {i} iElem {iElem}'
                # check if not already processed element nodes
                # check for the 1st node of the element
                if elemNodes[iElem,0] == -9999:
                    # print("Operating on Elem %d and Node %d" %(iElem,0))
                    debugMessage += f'Operating on Elem {iElem} and Node {0}'
                    idxNode = idxNode+1
                    # set the new node to the element
                    elemNodes[iElem,0] = idxNode
                    geomNode = geom1.strucNodes[i,0]
                    # print("\tgeomNode", geomNode)
                    debugMessage += f'\tgeomNode {geomNode}'
                    nodeCoords[idxNode,:] = geom1.nodeCoords[geomNode]
                    # check which strucs are connected to this node
                    # find the strucs with same geom node as start of another elem
                    idxConnStruc = np.argwhere(geom1.strucNodes[:,0]==geomNode)
                    # print("\tidxConnStruc startNode", idxConnStruc)
                    debugMessage += f'\tidxConnStruc startNode {idxConnStruc}'
                    self.processNodeConnElems(elemNodes,\
                        idxConnStruc,i,iElem,iNode=0,nodeOther="firstNode")
                    # process last node connected
                    idxConnStruc = np.argwhere(geom1.strucNodes[:,1]==geomNode)
                    # print("\tidxConnStruc endNode", idxConnStruc)
                    debugMessage += f'\tidxConnStruc endNode {idxConnStruc}'
                    self.processNodeConnElems(elemNodes,\
                        idxConnStruc,i,iElem,iNode=0,nodeOther="lastNode")
                # use the struc end node
                if elemNodes[iElem,1] == -9999:
                    # print("Operating on Elem %d and Node %d" %(iElem,1))
                    debugMessage += f'Operating on Elem {iElem} and Node {1}'
                    idxNode = idxNode+1
                    elemNodes[iElem,1] = idxNode
                    geomNode = geom1.strucNodes[i,1]
                    nodeCoords[idxNode,:] = geom1.nodeCoords[geomNode]
                    # check which strucs are connected to this node
                    # find the strucs with same geom node as start of another elem
                    idxConnStruc = np.argwhere(geom1.strucNodes[:,0]==geomNode)
                    # print("\tidxConnStruc startNode", idxConnStruc)
                    debugMessage += f'\tidxConnStruc startNode {idxConnStruc}'
                    self.processNodeConnElems(elemNodes,\
                        idxConnStruc,i,iElem,iNode=1,nodeOther="firstNode")
                    # find for the end names
                    idxConnStruc = np.argwhere(geom1.strucNodes[:,1]==geomNode)
                    # print("\tidxConnStruc endNode", idxConnStruc)
                    debugMessage += f'\tidxConnStruc endNode {idxConnStruc}'
                    self.processNodeConnElems(elemNodes,\
                        idxConnStruc,i,iElem,iNode=1,nodeOther="lastNode")

        if printDebugMessage:
            print(debugMessage)

        # delete the nan rows of node coords (all rows with nan)
        nodeCoords = nodeCoords[~np.isnan(nodeCoords).all(axis=1)]
        self.elemNodes = elemNodes
        self.nodeCoords = nodeCoords


    def processNodeConnElems(self,elemNodes, \
        idxConnStruc,iStruc,iElem,iNode=1,nodeOther="firstNode"):
        if nodeOther=="firstNode":
            for j in range(len(idxConnStruc)):
            # check for any other structure then the starting point
                if idxConnStruc[j] != iStruc:
                    # assign the node to the elem
                    jStruc = idxConnStruc[j]
                    # print("\tjStruc",jStruc)
                    jElem = self.strucElemConn[jStruc,1]
                    # print("\tjElem", jElem)
                    elemNodes[jElem,0] = elemNodes[iElem,iNode]
                else:
                    # print("\tskipping...")
                    pass
        if nodeOther=="lastNode":
            for j in range(len(idxConnStruc)):
                # check for any other structure then the starting point
                if idxConnStruc[j] != iStruc:
                    # if end node of struc then set last node of last elem
                    jStruc = idxConnStruc[j]
                    # print("\tjStruc",jStruc)
                    lastElemStruc = self.strucElemConn[jStruc,0] # index 0 is num of elems
                    # print("\tlastElemStruc", lastElemStruc)
                    jElem = self.strucElemConn[jStruc,lastElemStruc] # 0 idx is number of elems (so last idx)
                    # print("\tjElem", jElem)
                    elemNodes[jElem,1] = elemNodes[iElem,iNode]
                else:
                    # print("\tskipping...")
                    pass



    def writeVTK(self):
        #%% write the vtk data 
        # os.mkdir("."+os.path.sep+"vtk")
        vtkdir = os.getcwd()+os.path.sep+"vtk"
        if not os.path.isdir(vtkdir):
            os.mkdir(vtkdir)

        vtkf = open(vtkdir+os.path.sep+"mesh.vtk",'w')

        print('Writing mesh.vtk to %s ...\n' %(vtkdir))

        nodeCoords = self.nodeCoords
        elemNodes = self.elemNodes
        strucElemConn = self.strucElemConn
        nStrucs = self.nStrucs

        vtkf.write('# vtk DataFile Version 4.0\n')
        vtkf.write('vtk output\n')
        vtkf.write('ASCII\n')
        vtkf.write('DATASET UNSTRUCTURED_GRID\n')
        vtkf.write('POINTS %d float\n' %(np.shape(nodeCoords)[0]))

        nNodes = np.shape(nodeCoords)[0]

        for i in range(nNodes):
            vtkf.write('%f %f %f \n' %(nodeCoords[i,0],nodeCoords[i,1],0.0))

        nElems = np.shape(elemNodes)[0]
        # no of cells, (no points for each cell + 1)*nCells
        vtkf.write('CELLS %d %d\n' %(nElems, int(3*nElems)))


        for i in range(nElems):
            # points for each cell (start and end node)
            vtkf.write('%d %d %d\n' %(2, elemNodes[i,0], elemNodes[i,1]))

        # cell type, number 3 for 1d elements
        vtkf.write('CELL_TYPES %d\n' %(nElems))
        for i in range(nElems):
            # points for each cell (start and end node)
            vtkf.write('%d\n' %(3))


        # write the node numbers
        vtkf.write('POINT_DATA %d\n' %(nNodes))
        vtkf.write('SCALARS nodes integer %d\n' %(1))
        vtkf.write('LOOKUP_TABLE default \n')
        for i in range(nNodes):
            vtkf.write('%d\n' %(i))

        # write the cells/elems data
        vtkf.write('CELL_DATA %d\n' %(nElems))
        vtkf.write('SCALARS elements integer %d\n' %(1))
        vtkf.write('LOOKUP_TABLE default \n')
        for i in range(nElems):
            vtkf.write('%d\n' %(i))

        # write the cells of the geom structures
        vtkf.write('SCALARS structures integer %d\n' %(1))
        vtkf.write('LOOKUP_TABLE default \n')
        for i in range(nStrucs):
            for _ in range(strucElemConn[i,0]):
                vtkf.write('%d\n' %(i))

        print('\t... done writing VTK!\n')



