"""Inputs class
===============

.. autoclass:: Inputs
   :members:
   :private-members:
   
   .. automethod:: __init__
"""


import pandas as pd
import os
import json5
import json

class Inputs():
    """
    Class to parse/gather/process inputs
    """
    def __init__(self, infile: str=''):
        """
        Initiate the object of type Inputs

        Args:
            infile (str): input .json file
        Returns:
            None
        """
        
        self.mesh = pd.DataFrame()
        self.geom = pd.DataFrame()
        self.strucParams = pd.DataFrame()
        self.staticPointLoads = pd.DataFrame()
        self._load(infile=infile)

    def _load(self, infile: str=''):
        """
        Load the inputs from the input json file

        Args:
            infile (str): input .json file
        Returns:
            None
        """

        if not infile:
            indir = os.getcwd()
            infile = indir + os.path.sep + "inputs.json"

        if os.path.isfile(infile):
            # read using json5 to handle comments
            with open(infile) as f:
                jsonDat = json5.load(f)
            # move back to standar json string
            jsonDatString = json.dumps(jsonDat)
            # load the input data in a dataframe
            indf = pd.read_json(jsonDatString)
            # indf = pd.read_json(infile)
            self.infile = infile
        else:
            raise RuntimeError(f'Input file {infile} is missing!')

        idxStrucs = indf.index[indf.Name=="Structures"]
        stdf = pd.DataFrame(indf.loc[idxStrucs].Info.iloc[0])

        idxMesh = indf.index[indf.Name=="Mesh"]
        msdf = pd.DataFrame(indf.loc[idxMesh].Info.iloc[0])

        self.mesh = msdf
        self.geom = stdf

        idxStrucParams = indf.index[indf.Name=="StructureParams"]
        spdf = pd.DataFrame(indf.loc[idxStrucParams].Info.iloc[0])
        self.strucParams = spdf

        idxBCs = indf.index[indf.Name=="BoundaryConditions"]
        bcdf = pd.DataFrame(indf.loc[idxBCs].Info.iloc[0])
        self.bcs = bcdf

        self.staticPointLoads = \
            self.inputsToSubDataFrame(indf, "StaticPointLoads")

    def reload(self):
        '''
        Reload the inputs, using the same file
        '''
        self.load(infile=self.infile)

    def inputsToSubDataFrame(self, indf, inputName: str):
        """
        Extract a seperate dataframe for required input from the main input
        dataframe (created by parsing the input json file)
        """
        # get the index of the required sub input
        idxInput = indf.index[indf.Name==inputName]
        if not len(idxInput):
            return pd.DataFrame()
        inputSubdf = pd.DataFrame(indf.loc[idxInput].Info.iloc[0])
        return inputSubdf
