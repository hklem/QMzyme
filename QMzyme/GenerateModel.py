###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
Module in charge of generating the QMzymeModel given a starting structure. 
"""

from QMzyme.RegionBuilder import RegionBuilder
from QMzyme.QMzymeModel import QMzymeModel
import os
import QMzyme.MDAnalysisWrapper as MDAwrapper
from QMzyme.utils import translate_selection
from QMzyme.truncation_schemes import truncation_schemes

class GenerateModel(QMzymeModel):

    def __init__(self, *args, name=None, universe=None, **kwargs):
        """
        GenerateModel can be instantiated with an MDAnalysis Universe directly,
        or any combination of parameters that MDAnalysis.core.universe.Universe
        accepts to create a Universe i.e., (example.prmtop, example.dcd, dt=5).
        See https://userguide.mdanalysis.org/stable/universe.html for details.

        :param name: Name of QMzymeModel.
        :type name: str, default=None
        :param universe: MDAnalysis Universe object.
        :type universe: MDAnalysis.core.universe.Universe, default=None
        """
        if universe is None:
            universe = MDAwrapper.init_universe(*args, **kwargs)
        self.universe = universe
        if name is None:
            name = os.path.basename(self.universe.filename).split('.')[0]
        model = QMzymeModel(name, universe)
        self.__dict__.update(model.__dict__)


    def __repr__(self):
        return f"<ModelBuilder: Current QMzymeModel built from {self.starting_structure} contains {self.n_regions} region(s)>"


    def set_catalytic_center(self, selection):
        """
        Method to create a QMzymeRegion called 'catalytic_center'. Accepted input
        includes (i) str that can be interpreted by the MDAnalysis selection 
        command, (ii) an MDAnalysis.core.groups.AtomGroup, (iii) a QMzyme.QMzymeRegion.
        """
        self.set_region('catalytic_center', selection)
        return self.regions[-1]


    def set_region(self, region_name='no_name', selection=None):
        """
        Method to a QMzymeRegion. Accepted input includes (i) str that can be 
        interpreted by the MDAnalysis selection command, (ii) an 
        MDAnalysis.core.groups.AtomGroup, (iii) a QMzyme.QMzymeRegion.
        """
        selection = translate_selection(selection, self.universe)
        region_builder = RegionBuilder(region_name)
        #region = region_builder.init_atom_group(selection).get_region()
        region_builder.init_atom_group(selection)
        region = region_builder.get_region()
        self.add_region(region)
        return self.regions[-1]
    

    def truncate_region(self, region, truncation_scheme='CA_terminal'):
        """
        Method to truncate a QMzymeRegion. This will create a new region, and leave
        the original region unchanged.
        """
        new_region = truncation_schemes[truncation_scheme](region)
        self.add_region(new_region)
        return new_region

    def remove_region(self, region_index):
        """
        Method to remove a region from the QMzymeModel.
        """
        del self.regions[region_index]
