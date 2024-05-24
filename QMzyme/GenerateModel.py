###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
Module in charge of generating the QMzymeModel given a starting structure. 
"""

from QMzyme.QMzymeModel import QMzymeModel
import os
import QMzyme.MDAnalysisWrapper as MDAwrapper
from QMzyme.utils import make_selection
from QMzyme.TruncationSchemes import CA_terminal
from QMzyme.CalculateModel import CalculateModel
from QMzyme.Writers import Writer


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
        CalculateModel._reset()
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
        self.set_region(selection=selection, name='catalytic_center')


    def set_region(self, selection, name=None, **kwargs):
        """
        Method to a QMzymeRegion. Accepted input includes (i) str that can be 
        interpreted by the MDAnalysis selection command, (ii) an 
        MDAnalysis.core.groups.AtomGroup, or (iii) a QMzyme.QMzymeRegion.
        """
        region = make_selection(selection, model=self, name=name, **kwargs)
        self.add_region(region)
    
    
    def truncate_region(self, region, scheme=CA_terminal, name=None):
        """
        Method to truncate a QMzymeRegion. This will create a new region, and leave
        the original region unchanged.
        """
        s = scheme(region=region, name=name)
        region = s.return_region()
        self.add_region(region)


    def write_input(self, filename=None, memory=None, nprocs=None):
        writer = "".join(CalculateModel.calculation)
        Writer(filename, writer, memory, nprocs)

