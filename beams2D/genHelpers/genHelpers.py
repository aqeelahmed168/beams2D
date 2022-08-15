"""GenHelpers class
===================

.. autoclass:: GenHelpers
   :members:
   :private-members:
"""

import numpy as np

class GenHelpers():
    """
    Class to handle general operations on the numpy arrays
    """
    def __init__(self):
        """
        Initiate the object of type GenHelpers (not required generally)
        The GenHelpers class contains static functions moslty
        """
        pass

    @staticmethod
    def idxOfAllNonZeroRows(K: np.ndarray, tol=1e-15):
        """
        Get index of rows where all col values are non-zero (entire row of non zeros)

        Args:
            M: Input np.ndarray
            tol: Limit below which the float is considered as 0

        Returns:
            np.ndarray of idx
        """
   
        zeroRowsMask = (np.absolute(K) < 1e-15).all(axis=0)
        nonZeroRowsIdx = np.argwhere(zeroRowsMask==False).flatten()
        # zeroRowsIdx = np.argwhere(zeroRowsMask==True).flatten()

        return nonZeroRowsIdx


    @staticmethod
    def idxOfAllZeroRows(K: np.ndarray, tol=1e-10):
        """
        Get index of rows where all col values are zero (entire row of zeros)

        Args:
            M: Input np.ndarray
            tol: Limit below which the float is considered as 0

        Returns:
            np.ndarray of idx
        """
   
        zeroRowsMask = (np.absolute(K) < tol).all(axis=0)
        # nonZeroRowsIdx = np.argwhere(zeroRowsMask==False).flatten()
        zeroRowsIdx = np.argwhere(zeroRowsMask==True).flatten()


        return zeroRowsIdx

