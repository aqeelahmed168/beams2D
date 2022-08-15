"""Modal class
==============

Use assembled system to perform modal analysis

.. autoclass:: Modal
   :members:

   .. automethod:: __init__
"""

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

# for type hint (instant of class)
from ..assembleSystem import Assemble

class Modal():

    def __init__(self, assembleSys: Assemble):
        """Init method

        Args:
            assembleSys (Assemble): object of type Assemble
        Returns:
            None
        """
        self.nModes =  20 # modes to compute
        self.sysB2D = assembleSys

    def computeModes(self, nModes=20):
        """Compute modes

        Args:
            nModes (optional) number of modes to compute
        Returns:
            None
        """

        # Solve for eigenvalues and eigenvectors
        [omegaSqr, eigVectors] = linalg.eig(self.sysB2D.effK, self.sysB2D.effM)

        omega = np.sqrt(omegaSqr)
        omegaIdx = np.argsort(omega)        
        omega = np.sort(omega)

        self.omega = omega[0:nModes]

        # sort the eigen vectors (1st mode moved to 1st col)
        eigVectors = eigVectors[:, omegaIdx]

        # mode for all DOFs (to have easy index of row of dof)
        # row index is the index of dof
        # col index is mode number
        dataTypeEigVec = eigVectors.dtype
        eigVecFull = np.zeros((self.sysB2D.nDofs, self.sysB2D.nDofs), dtype=dataTypeEigVec)
        eigVecFull[np.ix_(self.sysB2D.effDofs, self.sysB2D.effDofs)] += eigVectors
        # delete the cols for fixed dofs/0 dofs
        eigVecEff = eigVecFull[:, self.sysB2D.effDofs]

        self.eigVecComputed = eigVectors

        self.eigVecFull = eigVecFull
        # retain specific number of modes
        self.eigVecEff = eigVecEff[:,0:nModes]


    def plotMode(self, modeNum=1, scale=1):
        """Plot mode in Jupyter notebook

        Args:
            modeNum (optional) modes number to plot
            scale (optional) scaling of the max deformation
        Returns:
            None
        """

        # plot the mode
        nodeCoords = np.copy(self.sysB2D.nodeCoords)
        nodeDisp = np.zeros_like(self.sysB2D.nodeCoords)

        modeVector = self.eigVecEff[:,modeNum-1]

        # check the dofs with each node 
        # modify the node coordinate for each node based on dofs
        for i in range(self.sysB2D.nNodes):
            iNodeDofs = self.sysB2D.dofNodes[i,:]
            nodeDisp[i,0] = modeVector[iNodeDofs[0]]
            nodeDisp[i,1] = modeVector[iNodeDofs[1]]    

            
        # print(nodeDisp[:,1])
        dispNorm = np.max(np.linalg.norm(nodeDisp, axis=1))
        scale = scale*1/dispNorm

        # check the dofs with each node 
        # modify the node coordinate for each node based on dofs
        for i in range(self.sysB2D.nNodes):
            iNodeDofs = self.sysB2D.dofNodes[i,:]
            nodeCoords[i,0] += scale*modeVector[iNodeDofs[0]]
            nodeCoords[i,1] += scale*modeVector[iNodeDofs[1]]

        plt.plot(self.sysB2D.nodeCoords[:,0], self.sysB2D.nodeCoords[:,1],'.-')
        plt.plot(nodeCoords[:,0], nodeCoords[:,1],'.--')
        plt.show()