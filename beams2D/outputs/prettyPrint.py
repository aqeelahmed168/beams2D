"""PrettyPrint class
====================

.. autoclass:: PrettyPrint
   :members:
   :private-members:
"""

import pandas as pd
from ..inputs.inputs import Inputs

class PrettyPrint():
    """
    Class to print matrices/vectors in notebooks
    """
    def __init__(self):
        """
        Initiate the object of type Outputs (not required generally)
        The Outputs class contains static functions moslty
        """
        pass

    @staticmethod
    def printMatrix(M, fontsize=10, text=''):
        """
        Print matrix in jupyter lab (notebooks)

        Args:
            M: Matric/vector to print
            fontsize: (optoinal) font size
            text: (optional) title of the output

        Returns:
            None
        """
        df = pd.DataFrame(M)
        ft = str(fontsize)
        display(df.style.set_table_attributes('style="font-size:'+ ft + 'px"').set_caption(text))

    @staticmethod
    def printInputs(inputs: Inputs, name='', fontsize=12, title=None):
        """
        Print inputs object attributes with name

        Args:
            inputs (Inputs): Inputs object created using inputs.json file
            name: name of the input sub-group ('geom, 'mesh', 'bcs')
            title: (optional) title of the output

        Returns:
            None
        """
        if not name:
            raise RuntimeError("name must be valid attribute")
        df = getattr(inputs,name)
        ft = str(fontsize)
        display(df.style.set_table_attributes('style="font-size:'+ ft + 'px"').set_caption(title))

