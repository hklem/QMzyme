###############################################################################
# Code written by Heidi Klem while at
# Colorado State University as a graduate student
# in the Paton and McCullagh groups and at the
# National Institute of Standards and Technology
# as an NRC Postdoc (Fed).
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
Class in charge of creating the QM enzyme models given a starting structure. 
This module can be used to call CalculateModel to generate QM calculation 
input files.

Required input
...............

    *   Prepared starting structure in PDB format.
"""

from QMzyme import MDAnalysisWrapper
from QMzyme.ModelBuilder import ModelBuilder
from QMzyme.RegionBuilder import RegionBuilder
import os

class GenerateModel:

    def __init__(self, starting_structure, name=None):
        """
        Constructor method for GenerateModel class.

        :param pdb_file: PDB file name.  
        :type pdb_file: str, required
        :param id: Name to associate structure with, defaults to file name. Will be used in naming of any files produced.
        :type id: str
        """
        if name is None:
            name = os.path.basename(starting_structure).split('.')[0]

        self.name = name
        self.structure_file = starting_structure

        model_builder = ModelBuilder.init_model(starting_structure, name)
        self.starting_structure = starting_structure
        self._dict__ = model_builder.get_model().__dict__

        # Raise warning if hydrogens are not present- can cause issues with truncate methods
        Hs_present = False
        for atom in self.starting_structure.atoms:
            if atom.element == 'H':
                Hs_present = True
                break
        if Hs_present == False:
            raise UserWarning("Structure does not contain hydrogens. Please pre-process this structure. We recommend using tLeap. See Documentation.")
    
    # def __repr__(self):
    #     return f"<QMzyme id={self.name}>"

    def set_catalytic_center(
        self,
        selection):        
        '''
        Function to define the center of the QMzyme model. This is
        typically the ligand/substrate. To add onto an already defined catalytic center 
        run set_catalytic_center again with the new selection and set 'overwrite=False'. 

        :param selection: Selection of atoms to be made- based on `MDAnalysis selection command language <https://docs.mdanalysis.org/stable/documentation_pages/selections.html>`_.
        :type selection: str, required`

        :param overwrite: To clear any existing catalytic_center definition. Set to False if you would like to append current catalytic_center. Defaults to True.
        :type overwrite: bool, default: overwrite=True
        '''
        catalytic_center = self.starting_structure.select_atoms(selection)

        self.regions.append()
        self.init_region('catalytic_center', catalytic_center) 
        #self.set_region('catalytic_center', catalytic_center)