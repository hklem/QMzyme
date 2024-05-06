"""
Concrete builder class to construct QMzymeModel product.
"""

from QMzyme.MDAnalysisWrapper import init_universe
from QMzyme.QMzymeModel import QMzymeModel
from QMzyme.RegionBuilder import RegionBuilder
from QMzyme.QMzymeRegion import QMzymeRegion

class ModelBuilder:
    # def __init__(self, name, universe):
    #     self.model = QMzymeModel(name, universe)

    def __init__(self):
        self.model = None

    def init_model(self, name, universe):
        self.model = QMzymeModel(name, universe)
        return self.model

    def init_region(self, region_name, atom_group):
        if atom_group.__class__ == QMzymeRegion:
            region = atom_group
        else:
            region = RegionBuilder(region_name, atom_group).get_region()
        self.model.add_region(region)
        return self.model.regions[-1]
        #return self

    def get_model(self):
        return self.model

    
