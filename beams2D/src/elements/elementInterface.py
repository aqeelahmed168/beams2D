"""Element Interface Class
==========================

Abstract class for elements

.. autoclass:: ElementInterface
   :members:
   :private-members:

"""

from abc import ABCMeta, abstractmethod

class ElementInterface(metaclass=ABCMeta):
    "The Element Interface"
    
    @staticmethod 
    @abstractmethod
    def K(EI: float, L: float, EA: float=0.0):
        """Abstract method K

        Args:
            EI (float): E*I for the element
            L (float): Length of the element
            EA (float): Optional E*A for the element (depending on the element)
        Returns:
            Element stiffness matrix in element reference frame
        """
        pass

    @staticmethod 
    @abstractmethod
    def M(rhoA: float, L: float):
        """Abstract method M

        Args:
            massPerUnitLength (float): Mass per unit length of the element
            L (float): Length of the element
        Returns:
            None
        """
        pass

    @staticmethod 
    @abstractmethod
    def ML(rhoA: float, L: float):
        """Abstract method M (Lumped Mass Matrix)

        Args:
            massPerUnitLength (float): Mass per unit length of the element
            L (float): Length of the element
        Returns:
            None
        """
        pass

    @staticmethod
    @abstractmethod
    def locToGlbTransformation(theta: float):
        """Abstract method to define the transformation from element local
        coordinates to global coordinates

        Args:
            theta (float): angle in radians betweem the element local reference 
                frame and the global reference frame
        Returns:
            None
        """
        pass
