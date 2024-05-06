from QMzyme.ModelBuilder import ModelBuilder
from QMzyme.QMzymeModel import QMzymeModel
import os
import QMzyme.MDAnalysisWrapper as MDAwrapper
from QMzyme.utils import translate_selection
from QMzyme.truncation_schemes import truncation_schemes
import copy

class GenerateModel(QMzymeModel):
    "The Director, building a complex representation."
    def __init__(self, file, name=None):
        if name is None:
            name = os.path.basename(file)
        universe = MDAwrapper.init_universe(file)
        #model = ModelBuilder(model_name, universe).get_model()
        self.__model_builder = ModelBuilder()
        model = self.__model_builder.init_model(name, universe)
        self.__dict__.update(model.__dict__)

        # model = ModelBuilder(model_name, universe).get_model()
        # self.model = model
        # self = ModelBuilder(model_name, universe).get_model()
    
    def set_catalytic_center(self, selection):
        #selection = translate_selection(selection, self)
        self.set_region('catalytic_center', selection)
        return self.regions[-1]
    
    def set_region(self, region_name, selection):
        selection = translate_selection(selection, self)
        self.__model_builder.init_region(region_name, selection)
        #self.regions.append(selection)
        # self.init_region(region_name, selection)
        return self.regions[-1]
    
    def truncate_region(self, region, truncation_scheme='CA_terminal'):
        new_region = truncation_schemes[truncation_scheme](region)
        return new_region
