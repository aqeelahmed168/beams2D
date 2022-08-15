"""Assemble class
=================

Assemble global Degree of Freedoms from the mesh data.
Build M, C, and K using specfic elements

.. autoclass:: Assemble
   :members:

   .. automethod:: __init__
"""

import numpy as np
from scipy.spatial.distance import cdist
from scipy import linalg
import os

# for type hint (instant of class)
from ..mesh.mesh import Mesh
from ...inputs.inputs import Inputs
# for building M, C, K
from ..elements.elementFactory import ElementFactory
# for getting indices on non-zero rows
from ...genHelpers.genHelpers import GenHelpers

class Assemble():

    def __init__(self, mesh: Mesh):
        """Init method

        Args:
            mesh (Mesh): object of type Mesh
        Returns:
            None
        """
        self.nodeCoords =  mesh.nodeCoords
        self.elemNodes = mesh.elemNodes
        self.strucElemConn = mesh.strucElemConn
        self.nStrucs = mesh.nStrucs
        self.K = np.array([])
        self.M = np.array([])
        # compute the DoF
        self._DoF()

    def _DoF(self):
        """Assemble all the dofs.
        This could depend on element type,
        but only beams are used with 3 Dofs for each node

        Args:
            None
        Returns:
            None
        """
        nNodes = np.shape(self.nodeCoords)[0]
        self.nNodes = nNodes
        # each node has 3 dofs (fixed for now)
        maxDofNode = 3 # ux, uy, az

        # global dofs related to nodes
        dofNodes = np.zeros((nNodes,maxDofNode),dtype=int)
        self.dofNodes = dofNodes
        # each node could have more or less then 3 dofs?
        elemNodes = self.elemNodes

        idof = 0
        # assign dofs for to each node
        for i in range(nNodes):
            dofNodes[i,:] = np.arange(idof,idof+maxDofNode,1,dtype=int)
            idof = idof+maxDofNode

        nTotElems = np.shape(elemNodes)[0]
        self.nElems = nTotElems
        # assign the global dofs to elements
        maxDofElem = 6
        nElems = nTotElems
        dofElems = np.zeros((nElems,maxDofElem),dtype=int)

        for i in range(nElems):
            iElemNodes = elemNodes[i,:]
            # might  need a change for different number of nodes/element
            dofElems[i,:] = dofNodes[iElemNodes,:].flatten()

        self.dofElems = dofElems
        self.nDofs = np.max(dofNodes)+1
        # link the element local dofs with the global dofs ...

    def getStrucParams(self, stdf):
        pass

    def buildGlobalMK(self, spdf: Inputs, elemType="EBB", lumpedMass=False):
        """Build the global M/K matrix based on element type

        Args:
            spdf (Inputs): object of type Geom
            elemTypes (str): Name of element type (EBB/EBF etc)
            lumpedMass (bool): Switch whether to used lumped mass matrix
        Returns:
            None
        """

        # get the structural properties/params
        # get the default values
        paramsDefaultIdx = spdf.structure.index[spdf['structure']==-1]
        if not len(paramsDefaultIdx):
            # check the values for strucutre 0 then
            paramsDefaultIdx = spdf.structure.index[spdf['structure']==0]
            # raise RuntimeWarning('Using default values based on strucutre 0')
            if not len(paramsDefaultIdx):
                errMessage = ('Default values not specified, specify in '
                    'StructureParams by setting values for strucutre with '
                        'index -1 or 0 (1st structure)')
                raise RuntimeError(errMessage)
        paramsDefaultIdx = paramsDefaultIdx.to_numpy().item()

        A = spdf.A[paramsDefaultIdx]
        E = spdf.E[paramsDefaultIdx]
        I = spdf.I[paramsDefaultIdx]
        rho = spdf.rho[paramsDefaultIdx]

        self.spdf = spdf

        nTotElems = np.shape(self.elemNodes)[0]
        # assing default values to all
        elemsA = A*np.ones((nTotElems,1))
        elemsE = E*np.ones((nTotElems,1))
        elemsI = I*np.ones((nTotElems,1))
        elemsRho = rho*np.ones((nTotElems,1))

        elemsTheta = np.full((nTotElems,1), np.nan)

        # check if node coords for a given struc are in the input.json
        # it should be for nElems > 1
        # find the index of the structure in the input
        strucElemConn = self.strucElemConn

        # build E/I/A for elems
        for iStruc in range(self.nStrucs):
            nElems = strucElemConn[iStruc,0]
            paramsStrucIdx = spdf.structure.index[spdf['structure']==iStruc]
            # check if the structure is in the list of inputs
            if not len(paramsStrucIdx):
                continue
            else:
                paramsStrucIdx = paramsStrucIdx.to_numpy().item()
            # if not-empty then update the values
            for j in range(nElems):
                iElem = strucElemConn[iStruc,j+1]
                elemsA[iElem] = spdf.A[paramsStrucIdx]
                elemsE[iElem] = spdf.E[paramsStrucIdx]
                elemsI[iElem] = spdf.I[paramsStrucIdx]
                elemsRho[iElem] = spdf.rho[paramsStrucIdx]
                

        self.elemsA = elemsA
        self.elemsE = elemsE
        self.elemsI = elemsI
        self.elemsRho = elemsRho

        elementsK = np.zeros((nTotElems,6,6), dtype=float)
        elementsM = np.zeros((nTotElems,6,6), dtype=float)
        #elemsTransMatrix = np.zeros((nTotElems,6,6), dtype=float)
        # elementsLocalK = np.zeros((nTotElems,6,6))
        # for now only one element type is supported in given simulation
        # could try to add different element type based on inputs
        element1 = ElementFactory.getElement(elemType)


        if elemType=='EBB':
            elementsLocalK = np.zeros((nTotElems,4,4), dtype=float)
            elementsLocalM = np.zeros((nTotElems,4,4), dtype=float)
            elemsTransMatrix = np.zeros((nTotElems,4,6), dtype=float)
        else:
            elementsLocalK = np.zeros((nTotElems,6,6), dtype=float)
            elementsLocalM = np.zeros((nTotElems,6,6), dtype=float)
            elemsTransMatrix = np.zeros((nTotElems,6,6), dtype=float)

        # save reference of the element type
        self.elementType = elemType

        dofNodes = self.dofNodes
        nDofs = np.max(dofNodes)+1
        glbK = np.zeros((nDofs,nDofs))
        glbM = np.zeros((nDofs,nDofs))

        elemNodes = self.elemNodes
        nodeCoords = self.nodeCoords

        for iStruc in range(self.nStrucs):
            nElems = strucElemConn[iStruc,0]
            for j in range(nElems):
                # get the coordinates of the nodes
                # compute the element length L
                jElem = strucElemConn[iStruc,j+1]
                jElemNodes = elemNodes[jElem,:]
                jElemDofs = self.dofElems[jElem,:]

                jElemSNCoords = nodeCoords[jElemNodes[0]]
                jElemENCoords = nodeCoords[jElemNodes[1]]

                Lj = np.linalg.norm(jElemENCoords-jElemSNCoords)
                # get the local K
                EI = elemsE[jElem]*elemsI[jElem]
                EA = elemsE[jElem]*elemsA[jElem]
                rhoA = elemsRho[jElem]*elemsA[jElem]

                jElemLocalK = element1.K(EI,Lj,EA=EA)
                if lumpedMass:
                    jElemLocalM = element1.ML(rhoA,Lj)
                else:
                    jElemLocalM = element1.M(rhoA,Lj)

                elementsLocalK[jElem,:,:] = jElemLocalK
                elementsLocalM[jElem,:,:] = jElemLocalM
                # compute the theta since not sure if aligned or not (radians)
                # w.r.t to unit vector in +ive x [1 0 0]
                jUnitVec = self.unitVector(jElemENCoords-jElemSNCoords)
                cosTheta = np.dot([1,0],jUnitVec)
                sinTheta = np.cross([1,0],jUnitVec)
                jThetaZ = np.arctan2(sinTheta, cosTheta)
                elemsTheta[jElem] = jThetaZ
                jElemT = element1.locToGlbTransformation(jThetaZ)
                jElemK = jElemT.T@jElemLocalK@jElemT
                jElemM = jElemT.T@jElemLocalM@jElemT
                elementsK[jElem,:,:] = jElemK
                elementsM[jElem,:,:] = jElemM
                elemsTransMatrix[jElem,:,:] = jElemT
                # assign to the glbK matrix (hard coded for 6x6 K)
                # sum the contributions
                # https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
                glbK[np.ix_(jElemDofs,jElemDofs)] += jElemK
                glbM[np.ix_(jElemDofs,jElemDofs)] += jElemM
        
        self.elemsTheta = elemsTheta
        self.elementsLocalK = elementsLocalK
        self.elemsTransMatrix = elemsTransMatrix
        self.elementsK = elementsK
        self.elementsM = elementsM
        self.K = glbK
        self.M = glbM

    def applyBoundaryConditions(self, bcdf: Inputs, retainNonZeroDofsOnly=True, tol=1e-10):
        """
        Apply boundary condtions. Remove the fixed DOFs.
        Return index of Global DOF to be removed. Compute  effecive DOFs.
        Compute effective K and M deleting fixed DOFs.

        Args:
            bcdf (Inputs): attribute bcs of object Inputs
            retainNonZeroDofsOnly : Remove DOFs with all rows 0 in K and M
            tol : Tolerance below which element of K/M is considered as zero
        Returns:
            None
        """
        # get the strucutures with input boundary conditions
        # get the nodes with given boundary conditions
        # find the global index of the nodes
        # in case of fixed remove all the DOFs
        # in case of pinned retain the rotational
        # loop over strucutures and mark the global fix/pinned nodes
        fixedNodes = np.array([], dtype=int)
        pinnedNodes = np.array([], dtype=int)

        # build list of fixed/pinned nodes
        for iStruc in range(self.nStrucs):
            bcStrucIdx = bcdf.structure.index[bcdf['structure']==iStruc]
            # check if the structure is in the list of inputs
            if not len(bcStrucIdx):
                continue
            else:
                bcStrucIdx = bcStrucIdx.to_numpy().item()
            # check if the start/end node is specified as fixed/pinned/free
            # get the element and node number
            iStartElem = self.strucElemConn[iStruc,1]
            iStartNode = self.elemNodes[iStartElem,0]
            nElems = self.strucElemConn[iStruc,0]
            iEndElem = self.strucElemConn[iStruc,nElems]
            iEndNode = self.elemNodes[iEndElem,1]

            if bcdf.startNode[bcStrucIdx] == 'free':
                pass
            elif bcdf.startNode[bcStrucIdx] == 'fixed':
                fixedNodes = np.append(fixedNodes, iStartNode)
            elif bcdf.startNode[bcStrucIdx] == 'pinned':
                pinnedNodes = np.append(pinnedNodes, iStartNode)
            else:
                raise RuntimeError('Valid boundary conditions are\
                    free, fixed and pinned')

            if bcdf.endNode[bcStrucIdx] == 'free':
                pass
            elif bcdf.endNode[bcStrucIdx] == 'fixed':
                fixedNodes = np.append(fixedNodes, iEndNode)
            elif bcdf.endNode[bcStrucIdx] == 'pinned':
                pinnedNodes= np.append(pinnedNodes, iEndNode)
            else:
                raise RuntimeError('Valid boundary conditions are\
                    free, fixed and pinned')
        
        self.fixedNodes = fixedNodes
        self.pinnedNodes = pinnedNodes

        # get the dofs for fixed/pinned nodes
        self.fixedDofs = self.dofNodes[fixedNodes,:]
        self.fixedDofs = self.fixedDofs.flatten()
        # for pinned only take 1st two DOFs and set those to fixed
        if len(self.pinnedNodes):
            pinnedDofs = self.dofNodes[pinnedNodes,0:2].flatten()
            self.fixedDofs = np.concatenate((self.fixedDofs,pinnedDofs))


        # self.effK = self.K(np._ix())
        dofList = np.arange(self.nDofs)
        effDofs = np.setdiff1d(dofList, self.fixedDofs)
        self.effDofs = effDofs
        effK = self.K[effDofs,:] # get all cols for effDofs rows
        effK = effK[:,effDofs] # get all rows for effDofs cols
        self.effK = effK
        effM = self.M[effDofs,:] # get all cols for effDofs rows
        effM = effM[:,effDofs] # get all rows for effDofs cols
        self.effM = effM

        if retainNonZeroDofsOnly:
            # find the row of K and M with all zeros
            # check all zero rows
            idxAllZeroRows = GenHelpers.idxOfAllZeroRows(effK, tol=tol)
            # this should match the also the ones in M
            # map the index of effDofs to the index of zero rows
            dofsZeroRows = self.effDofs[idxAllZeroRows]
            effDofs = np.setdiff1d(self.effDofs, dofsZeroRows)
            self.effDofs = effDofs
            effK = self.K[effDofs,:] # get all cols for effDofs rows
            effK = effK[:,effDofs] # get all rows for effDofs cols
            self.effK = effK
            effM = self.M[effDofs,:] # get all cols for effDofs rows
            effM = effM[:,effDofs] # get all rows for effDofs cols
            self.effM = effM

    def applyLoads(self, ldf: Inputs):
        """
        Apply the loads based on the input coordinates of the loads

        Args:
            ldf (Inputs): attribute loads of object Inputs
        Returns:
            None
        """

        # check for the nodes to apply the loads
        # find the index of the mesh node to apply loads
        # put that in the global matrix
        # delete the contrained DOFs (but maintian the global index)

        # find the nodes index based on the coordinate
        # and assign the forces and moments on the global
        staticLoadsVec = np.zeros((self.nDofs,1),dtype=float)
        nStaticLoads = len(ldf)

        for i in range(nStaticLoads):
            # check if the coordinate match any of the mesh node
            iCoords = np.array([ldf.nodeCoords[i]])
            # check if the node exists in the mesh nodes (using scipy)
            c = np.array(cdist(self.nodeCoords, iCoords) == 0)    
            if np.any(c):
                nodeIdx = np.argwhere(c.any(axis=1)).item()
                dofsINode = self.dofNodes[nodeIdx,:]
                # populate the loads vector based on dof
                staticLoadsVec[dofsINode[0]] = ldf.Fx[i]
                staticLoadsVec[dofsINode[1]] = ldf.Fy[i]
                staticLoadsVec[dofsINode[2]] = ldf.Mz[i]
            else:
                errMessage = ('Forcing must be applied on the mesh nodes',
                '\ncheck for the node coordinates of loads')
                raise RuntimeError(errMessage)
        
        self.staticLoads = staticLoadsVec
        self.staticLoadsEff = staticLoadsVec[self.effDofs]

    def computeStaticDeformations(self, useEffDofs=True):
        """
        Compute static loads using static forcing

        Args:
            useEffDofs (bool): switch to use effective dofs
        Returns:
            None
        """

        self.staticDeformations = np.zeros((self.nDofs,1), dtype=float)

        if useEffDofs:
            print('Solving for static deformations')
            x = linalg.solve(self.effK, self.staticLoadsEff)
        else:
            errMessage = ('Static deformation are suppoted only with ',
                '\neffective DOFs, apply by setting BCs')
            raise RuntimeError(errMessage)
        
        self.staticDeformations[self.effDofs] = x


    def writeStaticDeformationsVTK(self, filename='staticDeforms.vtk', scale=1):
        """
        Write VTK file for static deformations with optional scale

        Args:
            filename (str): filename (optional) for output vtk file
            scale: scale for deformation (optional)
        Returns:
            None
        """

        vtkdir = os.getcwd()+os.path.sep+"vtk"
        if not os.path.isdir(vtkdir):
            os.mkdir(vtkdir)

        vtkf = open(vtkdir+os.path.sep+filename,'w')

        print('Writing %s to %s ...\n' %(filename,vtkdir))

        nNodes = self.nNodes
        nodeCoords = self.nodeCoords
        uDef = np.zeros((nNodes,3),dtype=float)
        aDef = np.zeros((nNodes,3),dtype=float)

        # check the dofs with each node 
        # modify the node coordinate for each node based on dofs
        for i in range(nNodes):
            iNodeDofs = self.dofNodes[i,:]
            nodeCoords[i,0] += scale*self.staticDeformations[iNodeDofs[0]]
            nodeCoords[i,1] += scale*self.staticDeformations[iNodeDofs[1]]
            uDef[i,0] = scale*self.staticDeformations[iNodeDofs[0]]
            uDef[i,1] = scale*self.staticDeformations[iNodeDofs[1]]
            aDef[i,2] = scale*self.staticDeformations[iNodeDofs[2]]
            # add impact of rotation around z on coordinates?

        elemNodes = self.elemNodes
        strucElemConn = self.strucElemConn
        nStrucs = self.nStrucs

        vtkf.write('# vtk DataFile Version 4.0\n')
        vtkf.write('vtk output\n')
        vtkf.write('ASCII\n')
        vtkf.write('DATASET UNSTRUCTURED_GRID\n')
        vtkf.write('POINTS %d float\n' %(np.shape(nodeCoords)[0]))

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

        # write the translation deformations
        vtkf.write('VECTORS transDisp(x%d) float \n' %(scale))
        for i in range(nNodes):
            vtkf.write('%f %f %f \n' %(uDef[i,0],uDef[i,1],0.0))

        vtkf.write('VECTORS rotDisp(x%d) float \n' %(scale))
        for i in range(nNodes):
            vtkf.write('%f %f %f \n' %(aDef[i,0],aDef[i,1],aDef[i,2]))

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

    def computeElementReactions(self):
        """Compute for each node of the element internal reaction forces 
        and moments in the element reference frame

        Args:
            None
        Returns:
            None
        """
    
        # check if static deformations are available
        if not len(self.staticDeformations):
            errMsg = ('Compute static deformations before computing reactions')
            raise RuntimeError(errMsg)

        elemReactions = np.zeros((self.nElems,6), dtype=float)
        for i in range(self.nElems):
            elemT = self.elemsTransMatrix[i,:,:]
            dofElem = self.dofElems[i]
            elemDef = self.staticDeformations[dofElem]
            elemK = self.elementsLocalK[i,:,:]
            reacElem = elemK@(elemT@elemDef)
            elemReactions[i,:] = np.ravel(reacElem)

        self.elemReactions = elemReactions

    def glbC():
        pass

    def unitVector(self, v):
        vHat = v/np.linalg.norm(v)
        return vHat

    def rotMatrix(self, theta):
        T = np.array([\
                [np.cos(theta),  np.sin(theta)],\
                [-np.sin(theta), np.cos(theta)],\
                ])
        return T
    

    
    