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
from QMzyme.selection_schemes import selection_schemes
from QMzyme.CalculateModel import write_QM, write_QMMM

class GenerateModel(QMzymeModel):

    def __init__(self, *args, name=None, universe=None, frame=0, **kwargs):
        """
        GenerateModel can be instantiated with an MDAnalysis Universe directly,
        or any combination of parameters that MDAnalysis.core.universe.Universe
        accepts to create a Universe i.e., (example.prmtop, example.dcd, dt=5).
        See https://userguide.mdanalysis.org/stable/universe.html for details.

        :param name: Name of QMzymeModel.
        :type name: str, default=None
        :param universe: MDAnalysis Universe object.
        :type universe: MDAnalysis.core.universe.Universe, default=None
        :param frame: If trajectory was provided, specify a frame to base coordinates on
        :type frame: int, default=0
        """
        if universe is None:
            universe = MDAwrapper.init_universe(*args, **kwargs)
        if frame != 0:
            universe.trajectory[frame]
        self.universe = universe
        if name is None:
            name = os.path.basename(self.universe.filename).split('.')[0]
        model = QMzymeModel(name, universe)
        self.__dict__.update(model.__dict__)


    def __repr__(self):
        return f"<QMzymeModel built from {self.universe} contains {self.n_regions} region(s)>"


    def set_catalytic_center(self, selection):
        """
        Method to create a QMzymeRegion called 'catalytic_center'. Accepted input
        includes (i) str that can be interpreted by the MDAnalysis selection 
        command, (ii) an MDAnalysis.core.groups.AtomGroup, (iii) a QMzyme.QMzymeRegion.
        """
        self.set_region(layer=None, selection=selection, region_name='catalytic_center',)
        #return self.regions[-1]


    def set_region(self, layer, selection, region_name="None"):
        """
        Method to a QMzymeRegion. Accepted input includes (i) str that can be 
        interpreted by the MDAnalysis selection command, (ii) an 
        MDAnalysis.core.groups.AtomGroup, or (iii) a QMzyme.QMzymeRegion.
        """
        selection = translate_selection(selection, self.universe)
        region_builder = RegionBuilder(region_name)
        region_builder.init_atom_group(selection)
        region = region_builder.get_region()
        self.add_region(region)
        #return self.regions[-1]
    

    def truncate_region(self, region, truncation_scheme='CA_terminal'):
        """
        Method to truncate a QMzymeRegion. This will create a new region, and leave
        the original region unchanged.
        """
        new_region = truncation_schemes[truncation_scheme](region)
        self.add_region(new_region)
        #return new_region

    def selection_scheme(self, selection_scheme: str='distance_cutoff', **kwargs):
        """
        Method to call a selection scheme to identify what atoms to include in the QMzymeRegion.
        If selection_scheme='distance_cutoff', the keyword argument cutoff must also be supplied.
        """
        new_region = selection_schemes[selection_scheme](self, kwargs['cutoff'])
        self.add_region(new_region)


    def calculate(self, filename):
        write_file = {
            "QMQM": write_QM()
        }
        layers = {}
        methods = []
        for region in self.regions:
            if region.layer == None:
                continue
            if region.method == None:
                raise UserWarning(f"Calculation method not set for {region}.")
            layers[region.layer] = region
            methods.append(region.method[type])
        method = "".join(methods)

    def _check_layers(self):
        regions = [r for r in self.regions if r.layer != None]
        layers = []
        for r in regions:
            if r.layer in layers:
                raise UserWarning(f"Region {r} is already set as the {r.layer} layer.")

