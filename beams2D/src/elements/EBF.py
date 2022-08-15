"""EulerBernoulliFrame class
============================

.. autoclass:: EulerBernoulliFrame
   :members:
   :private-members:

"""

import numpy as np

from .elementInterface import ElementInterface

class EulerBernoulliFrame(ElementInterface):
    """
    Element with Euler-Bernoulli Beam including axial DOFs.
    Each node with 3 DOFs, uy, displacement transverse to the axis,
    ux, displacement along the longitudinal axis, and
    az, rotation around longitudinal axis
    """
    def __init__(self):
        pass

    def K(self, EI, L, EA=0.0) -> np.ndarray:
        """
        Element stiffness matrix in element reference frame

        .. math::

            [K] = 
            \\frac{EI}{L^3}
            \\begin{bmatrix} 
                \\frac{EA}{L}\\frac{L^3}{EI}  & 0   & 0    & -\\frac{EA}{L}\\frac{L^3}{EI} & 0   & 0    \\\ 
                0                             & 12  & 6L   & 0                             & -12 & 6L   \\\ 
                0                             & 6L  & 4L^2 & 0                             & -6L & 2L^2 \\\ 
                -\\frac{EA}{L}\\frac{L^3}{EI} & 0   & 0    & \\frac{EA}{L}\\frac{L^3}{EI}  & 0   &  0   \\\ 
                0                             & -12 & -6L  & 0                             & -12 & 6L   \\\ 
                0                             & 6L  & 2L^2 & 0                             & -6L & 4L^2 \\\ 
            \\end{bmatrix}


        Args:
            EI (float): E*I for the element
            L (float): Length of the element
            EA (float): E*A for the element (required non-zero for EBF)
        Returns:
            Element stiffness matrix in element reference frame
        """

        if (EA < 1e-10):
            raise RuntimeError("EA must be non zero (>1e-10) for EBF element")
    
        alpha = EI/L**3

        ke = np.array(
            [\
                [EA/L,  0,         0,            -EA/L, 0,          0           ],\
                [0,     12*alpha,  6*L*alpha,    0,     -12*alpha,  6*L*alpha   ],\
                [0,     6*L*alpha, 4*L**2*alpha, 0,     -6*L*alpha, 2*L**2*alpha],\
                [-EA/L, 0,         0,            EA/L,  0,          0           ],\
                [0,     -12*alpha, -6*L*alpha,   0,     12*alpha,   -6*L*alpha  ],\
                [0,     6*L*alpha, 2*L**2*alpha, 0,     -6*L*alpha, 4*L**2*alpha]\
            ],dtype=float)
        
        return ke


    def M(self, rhoA, L) -> np.ndarray:
        """
        Element mass matrix (consistent) in element reference frame.

        .. math::

            \\text{For bar element with DOF } \\{u_{x_1}, u_{x_2}\\}
            \\\\
            [M]_{bar} = 
                \\frac{\\rho A L}{6}
                \\begin{bmatrix} 
                    2 & 1  \\\ 
                    1 & 2  \\\ 
                \\end{bmatrix}
                \\\\
            \\text{For beam element with DOF } \\{u_{y_1}, a_{z_1}, u_{y_2}, a_{z_2}\\}
            \\\\
            [M]_{beam} = 
                \\frac{\\rho A L}{420}
                \\begin{bmatrix} 
                    156 & 22L   & 54   & -13L  \\\ 
                        & 4L^2  & 13L  & -3L^2 \\\ 
                        &       & 156  & -22L  \\\ 
                        &       &      & 4L^2  \\\ 
                \\end{bmatrix}
                \\\\
            \\text{For frame element with DOF } 
            \\{u_{x_1}, u_{y_1}, a_{z_1}, u_{x_2}, u_{y_2}, a_{z_2}\\}
            \\\\
            [M] = [M]_{bar} + [M]_{beam} = 
                \\frac{\\rho A L}{420}
                \\begin{bmatrix} 
                    140  & 0     & 0     & 70   & 0   & 0     \\\ 
                    0    & 156   & 22L   & 0    & 54   & -13L  \\\  
                    0    & 22L   & 4L^2  & 0    & 13L  & -3L^2 \\\ 
                    70   & 0     & 0     & 140  & 0    & 0     \\\ 
                    0    & 54    & 13L   & 0    & 156  & -22L  \\\ 
                    0    & -13L  & -3L^2 & 0    & -22L & 4L^2  \\\ 
                \\end{bmatrix}

        Args:
            rhoA (float): Mass per unit Length of the element (rho*A)
            L (float): Length of the element
        Returns:
            Element mass matrix in element reference frame
        """

        me = (rhoA*L/420.0)*np.array(
        [\
            [140,   0,         0,            70,    0,     0         ],\
            [0,     156,       22*L,         0,     54,    -13*L     ],\
            [0,     22*L,      4*L**2,       0,     13*L,  -3*L**2   ],\
            [70,    0,         0,            140,   0,     0         ],\
            [0,     54,        13*L,         0,     156,   -22*L     ],\
            [0,     -13*L,     -3*L**2,      0,     -22*L, 4*L**2    ]\
        ],dtype=float)
        
        return me

    def ML(self, rhoA, L) -> np.ndarray:
        """
        Element lumped mass matrix in element reference frame.

        .. math::

            \\text{For bar element with DOF } \\{u_{x_1}, u_{x_2}\\}
            \\\\
            [M]_{bar} = 
                \\frac{\\rho A L}{2}
                \\begin{bmatrix} 
                    1 & 0  \\\ 
                    0 & 1  \\\ 
                \\end{bmatrix}
                \\\\
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
            \\text{For frame element with DOF } 
            \\{u_{x_1}, u_{y_1}, a_{z_1}, u_{x_2}, u_{y_2}, a_{z_2}\\}
            \\\\
            [M] = [M]_{bar} + [M]_{beam} = 
                \\frac{\\rho A L}{2}
                \\begin{bmatrix} 
                    1    & 0    & 0   & 0   & 0    & 0 \\\ 
                    0    & 1    & 0   & 0   & 0    & 0 \\\  
                    0    & 0    & 0   & 0   & 0    & 0 \\\ 
                    0    & 0    & 0   & 1   & 0    & 0 \\\ 
                    0    & 0    & 0   & 0   & 1    & 0 \\\ 
                    0    & 0    & 0   & 0   & 0    & 0 \\\ 
                \\end{bmatrix}

        Args:
            rhoA (float): Mass per unit Length of the element (rho*A)
            L (float): Length of the element
        Returns:
            Element mass matrix in element reference frame
        """

        me = (rhoA*L/2)*np.array(
            [\
                [1,     0,     0,    0,    0,    0 ],\
                [0,     1,     0,    0,    0,    0 ],\
                [0,     0,     0,    0,    0,    0 ],\
                [0,     0,     0,    1,    0,    0 ],\
                [0,     0,     0,    0,    1,    0 ],\
                [0,     0,     0,    0,    0,    0 ],\
            ],dtype=float)
        return me


    def locToGlbTransformation(self, theta) -> np.ndarray:
        """
        Transformation matrix to move from element reference frame to 
        global reference frame

        * Element DOF = {ux1', uy1', az1', ux2', uy2', az2'}'
        * Global DOF = {ux1, uy1, az1, ux2, uy2, az2}'

        The transformation matrix

        .. math::

            [T] = 
            \\begin{bmatrix} 
                \\cos(\\theta)  & \\sin(\\theta) & 0 & 0 & 0               & 0  \\\ 
                -\\sin(\\theta) & \\cos(\\theta) & 0 & 0 & 0               & 0  \\\ 
                0               & 0              & 1 & 0 & 0               & 0  \\\ 
                0               & 0              & 0 & \\cos(\\theta)  & \\sin(\\theta) & 0 \\\ 
                0               & 0              & 0 & -\\sin(\\theta) & \\cos(\\theta) & 0 \\\ 
                0               & 0              & 0 & 0 & 0               & 1  \\\ 
            \\end{bmatrix}

        Args:
            theta: angle in radians between the element local axis (x'-y') 
                and the global axis (x-y)
        Returns:
            np.ndarray: Transformation matrix for the element

        """
        c = np.cos(theta)
        s = np.sin(theta)

        T = np.array([\
                [c,  s, 0,  0, 0, 0],\
                [-s, c, 0,  0, 0, 0],\
                [0,  0, 1,  0, 0, 0],\
                [0,  0, 0,  c, s, 0],\
                [0,  0, 0, -s, c, 0],\
                [0,  0, 0,  0, 0, 1],\
                ],dtype=float)
        return T