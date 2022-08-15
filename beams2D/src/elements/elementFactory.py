
"""Element Factory Class
========================

.. autoclass:: ElementFactory
   :members:
   :private-members:

"""

from .EBB import EulerBernoulliBeam
from .EBF import EulerBernoulliFrame


class ElementFactory:
    """
    Class to instantiate the element based on type
    """
    @staticmethod
    def getElement(elementType):
        "Static method to get element type"
        if elementType == 'EBB':
            return EulerBernoulliBeam()
        if elementType == 'EBF':
            return EulerBernoulliFrame()
        return RuntimeError("Please specify valid element type!")

