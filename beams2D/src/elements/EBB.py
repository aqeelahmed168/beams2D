"""EulerBernoulliBeam class
===========================

.. autoclass:: EulerBernoulliBeam
   :members:
   :private-members:

"""

import numpy as np

from .elementInterface import ElementInterface

class EulerBernoulliBeam(ElementInterface):
    """
    Euler-Bernoulli Beam Element.
    Each node with 2 DOFs, uy, displacement transverse to the axis and, 
    az, rotation around longitudinal axis
    """
    def __init__(self):
        pass

    def K(self, EI, L, EA=0):
        """
        Element stiffness matrix in element reference frame.

        .. math::

            [K] = 
                \\frac{EI}{L^3}
                \\begin{bmatrix} 
                    12  & 6L    & -12  & 6L   \\\ 
                        & 4L^2  & -6L  & 2L^2 \\\ 
                        &       & 12   & -6L  \\\ 
                        &       &      & 4L^2 \\\ 
                \\end{bmatrix}

        Args:
            EI (float): E*I for the element
            L (float): Length of the element
        Returns:
            Element stiffness matrix in element reference frame
        """

        ke = (EI/L**3)*np.array([\
            [12,   6*L,     -12,   6*L   ],\
            [6*L,  4*L**2,  -6*L,  2*L**2],\
            [-12,  -6*L,    12,    -6*L  ],\
            [6*L,  2*L**2,  -6*L,  4*L**2]\
            ], dtype=float)
        return ke
    
    def M(self, rhoA, L):
        """
        Element mass matrix in element reference frame.

        .. math::

            [M] = 
                \\frac{\\rho A L}{420}
                \\begin{bmatrix} 
                    156 & 22L   & 54   & -13L   \\\ 
                        & 4L^2  & 13L  & -3L^2 \\\ 
                        &       & 156  & -22L  \\\ 
                        &       &      & 4L^2 \\\ 
                \\end{bmatrix}

        Args:
            rhoA (float): Mass per unit Length of the element (rho*A)
            L (float): Length of the element
        Returns:
            Element mass matrix in element reference frame
        """
        me = (rhoA*L/420)*np.array([\
            [156,   22*L,    54,    -13*L  ],\
            [22*L,  4*L**2,  13*L,  -3*L**2],\
            [54,    13*L,    156,   -22*L  ],\
            [-13*L, -3*L**2, -22*L, 4*L**2 ]\
            ])
        return me


    def ML(self, rhoA, L):
        """
        Element lumped mass matrix in element reference frame.

        .. math::

            \\text{For beam element with DOF } \\{u_{y_1}, a_{z_1}, u_{y_2}, a_{z_2}\\}
            \\\\
            [M]_{beam} = 
                \\frac{\\rho A L}{2}
                \\begin{bmatrix} 
                    1   & 0    & 0    & 0  \\\ 
                    0   & 0    & 0    & 0  \\\ 
                    0   & 0    & 1    & 0  \\\ 
                    0   & 0    & 0    & 0  \\\ 
                \\end{bmatrix}
                \\\\

        Args:
            rhoA (float): Mass per unit Length of the element (rho*A)
            L (float): Length of the element
        Returns:
            Element mass matrix in element reference frame
        """
        me = (rhoA*L/2)*np.array([\
            [1,  0,   0,   0],\
            [0,  0,   0,   0],\
            [0,  0,   1,   0],\
            [0,  0,   0,   0]\
            ])
        return me


    def locToGlbTransformation(self, theta: float) -> np.ndarray:
        """

        Transformation matrix to move from element reference frame to 
        global reference frame

        * Element DOF = {uy1', az1', uy2', az2'}'
        * Global DOF = {ux1, uy1, az1, ux2, uy2, az2}'

        The transformation matrix

        .. math::

            [T] = 
            \\begin{bmatrix} 
                \\cos(\\theta)  & \\sin(\\theta) & 0 & 0 & 0               & 0  \\\ 
                0               & 0              & 1 & 0 & 0               & 0  \\\ 
                0               & 0              & 0 & -\\sin(\\theta) & \\cos(\\theta) & 0 \\\ 
                0               & 0              & 0 & 0 & 0               & 1  \\\ 
            \\end{bmatrix}

        Args:
            theta: angle in radians between the element local axis and the global axis (x-y)
        Returns:
            np.ndarray: Transformation matrix for the element

        """
        c = np.cos(theta)
        s = np.sin(theta)
        T = np.array([\
                [-s, c, 0,  0, 0, 0],\
                [0,  0, 1,  0, 0, 0],\
                [0,  0, 0, -s, c, 0],\
                [0,  0, 0,  0, 0, 1],\
                ])
        return T 

