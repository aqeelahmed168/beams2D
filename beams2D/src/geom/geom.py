"""Geom class
=============

.. autoclass:: Geom
   :members:
   :private-members:
   
   .. automethod:: __init__
"""

import numpy as np
import pandas as pd
import os

# for type hint (instant of class Geom)
from ...inputs.inputs import Inputs

class Geom:
    '''Class to build/process geometry'''
    def __init__(self, stdf: Inputs):
        """Init method

        Args:
            stdf (Inputs): attribute geom of the object Inputs
        Returns:
            None
        """
        self.stdf = stdf
    
    def process(self, printDebugMessage=0):
        """Process geometry, build the strucutres/nodes

        Args:
            printDebugMessage: flag to show additional outputs

        Returns:
            None
        """
        stdf = self.stdf
        nStrucs = len(stdf.index)

        strucProcessed = np.zeros(nStrucs,dtype=bool)
        # each structure contains 2 nodes
        strucNodes = np.zeros((nStrucs,2),dtype=int)

        # make list based on max possible nodes (delete later nans)
        nodeCoords = np.full((nStrucs+1,2),np.nan)

        debugMessage = ''

        nodei = int(-1)
        for i in range(nStrucs):
            # print("Operating on structure %d" %(i))
            debugMessage += 'Operating on structure {0}\n'.format(i)
            if not strucProcessed[i]:
                # print("\tStructure %d not processed yet" %(i))
                debugMessage += '\tStructure {0} not processed yet\n'.format(i)
                nodei = nodei+1
                strucNodes[i,0] = nodei
                nodeCoords[nodei,:] = stdf['startNode'][i]
                nodei = nodei+1
                strucNodes[i,1] = nodei
                nodeCoords[nodei,:] = stdf['endNode'][i]
                strucProcessed[i] = True
                # print('\tstructure nodes: ', strucNodes[i,:])
                debugMessage += '\tstructure nodes: {0}\n'.format(strucNodes[i,:])

            for j in range(len(stdf['enConnStrucs'][i])):
                connStruc = stdf['enConnStrucs'][i][j]
                # print("\tstruc:%d, connElem:%d" %(i, connStruc))
                debugMessage += '\tstruc:{0}, connElem:{1}'.format(i, connStruc)
                if strucProcessed[connStruc]:
                    continue
                else:
                    strucProcessed[connStruc] = True
                # check that endNode is indeed the startNode for connected 
                # element (input errors)
                if not np.allclose(stdf['endNode'][i], stdf['startNode'][connStruc]):
                    # print("end node (%f,%f) of structure %d is \
                    #     not the start node (%f,%f) of structure %d\
                    #         " %(stdf['endNode'][i][0], stdf['endNode'][i][1], i,
                    #             stdf['startNode'][connStruc][0], 
                    #             stdf['startNode'][connStruc][1], connStruc))
                    # within bracket line break by \ is not required
                    debugMessage += ("end node ({0},{1}) of structure {2} is \
                        not the start node ({3},{4}) of structure {5}\
                            \n".format(stdf['endNode'][i][0], 
                                stdf['endNode'][i][1], i,
                                stdf['startNode'][connStruc][0], 
                                stdf['startNode'][connStruc][1], connStruc))

                    raise RuntimeError("end node of structure is \
                        not the start node of connected structure")
                # print("\tOperating on connected structure %d" %(connStruc))
                debugMessage += '\tOperating on connected structure {0}\n'.format(connStruc)
                strucNodes[connStruc,0] = strucNodes[i,1]
                # check the end node of connected structure is not the start node 
                # of already processed structure
                # be careful of multiple connected strucutures
                enConnStrucs2ConnStruc = stdf['enConnStrucs'][connStruc]
                enConnElemProcessed = False
                for k in range(len(enConnStrucs2ConnStruc)):
                    enConnStruc2ConnStruc = enConnStrucs2ConnStruc[k]
                    # print("\t\tenConnElem2ConnElem", enConnElem2ConnElem)
                    debugMessage += '\t\enConnStruc2ConnStruc {0}\n'.format(enConnStruc2ConnStruc)
                    # if any of the end node connected element is already defined
                    # break the loop and use that element node
                    if strucProcessed[enConnStruc2ConnStruc]:
                        strucNodes[connStruc,1] = strucNodes[enConnStruc2ConnStruc,0]
                        enConnElemProcessed = True    
                        break
                # if no processed connected element found make a new node
                if not enConnElemProcessed:
                    nodei = nodei+1
                    strucNodes[connStruc,1] = nodei
                    nodeCoords[nodei,:] = stdf['endNode'][connStruc]
                    # nodei = nodei+1
                # print('\t\tconnElem nodes', strucNodes[connStruc,:])
                debugMessage += '\t\tconnElem nodes {0}\n'.format(strucNodes[connStruc,:])
                        

        if printDebugMessage:
            print(debugMessage)

        # delete the nan rows of node coords (all cols are nan)
        nodeCoords = nodeCoords[~np.isnan(nodeCoords).all(axis=1)]

        # put data in the geom object
        self.nStrucs = nStrucs
        self.strucNodes = strucNodes
        self.nodeCoords = nodeCoords

    def writeVTK(self):
        # write the vtk data 
        # os.mkdir("."+os.path.sep+"vtk")
        vtkdir = os.getcwd()+os.path.sep+"vtk"
        if not os.path.isdir(vtkdir):
            os.mkdir(vtkdir)

        vtkf = open(vtkdir+os.path.sep+"geom.vtk",'w')

        print('Wrting geom.vtk to %s ...\n' %(vtkdir))

        nodeCoords = self.nodeCoords
        strucNodes = self.strucNodes

        vtkf.write('# vtk DataFile Version 4.0\n')
        vtkf.write('vtk output\n')
        vtkf.write('ASCII\n')
        vtkf.write('DATASET UNSTRUCTURED_GRID\n')
        vtkf.write('POINTS %d float\n' %(np.shape(nodeCoords)[0]))

        nNodes = np.shape(nodeCoords)[0]

        for i in range(nNodes):
            vtkf.write('%f %f %f \n' %(nodeCoords[i,0],nodeCoords[i,1],0.0))

        nElems = np.shape(strucNodes)[0]
        # no of cells, (no points for each cell + 1)*nCells
        vtkf.write('CELLS %d %d\n' %(nElems, int(3*nElems)))


        for i in range(nElems):
            # points for each cell (start and end node)
            vtkf.write('%d %d %d\n' %(2, strucNodes[i,0], strucNodes[i,1]))

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
        vtkf.write('SCALARS structures integer %d\n' %(1))
        vtkf.write('LOOKUP_TABLE default \n')
        for i in range(nElems):
            vtkf.write('%d\n' %(i))

        print('\t... done writing!\n')



