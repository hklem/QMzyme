from QMzyme.RegionBuilder import RegionBuilder
from QMzyme.QMzymeModel import QMzymeModel
import os
import QMzyme.MDAnalysisWrapper as MDAwrapper
from QMzyme.utils import translate_selection
from QMzyme.truncation_schemes import truncation_schemes

class GenerateModel(QMzymeModel):
    "The Director, building a complex representation."
    def __init__(self, file, name=None):
        if name is None:
            name = os.path.basename(file).split('.')[0]
        universe = MDAwrapper.init_universe(file)
        model = QMzymeModel(name, universe)
        self.__dict__.update(model.__dict__)


    def __repr__(self):
        return f"<ModelBuilder: Current QMzymeModel built from {self.starting_structure} contains {self.n_regions} region(s)>"


    def set_catalytic_center(self, selection):
        self.set_region('catalytic_center', selection)
        return self.regions[-1]


    def set_region(self, region_name, selection):
        selection = translate_selection(selection, self.starting_structure)
        region_builder = RegionBuilder(region_name)
        #region = region_builder.init_atom_group(selection).get_region()
        region_builder.init_atom_group(selection)
        region = region_builder.get_region()
        self.add_region(region)
        return self.regions[-1]
    

    def truncate_region(self, region, truncation_scheme='CA_terminal'):
        new_region = truncation_schemes[truncation_scheme](region)
        return new_region
    

    def remove_region(self, region_index):
        del self.regions[region_index]


