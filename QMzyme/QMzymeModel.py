###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
Product of the ModelBuilder class.
"""

class QMzymeModel:
    # def __init__(self, name, starting_structure, regions=[]):
    def __init__(self, name, universe):
        self.name = name
        self.starting_structure = universe
        self.filename = universe.filename
        self.regions = []

    def __repr__(self):
        return f"<QMzymeModel {self.name} built from {self.starting_structure} contains {self.n_regions} region(s)>"

    @property
    def n_regions(self):
        return len(self.regions)
    
    def add_region(self, region):
         setattr(self, region.name, region)
         self.regions.append(region)

    def get_region_names(self):
         return [r.name for r in self.regions]
    
    def get_region(self, region_name=None):
        # if region_name in self.get_region_names():
        #     return self.regions[self.get_region_names().index(region_name)]
        try:
            return getattr(self,region_name)
        except:
            raise UserWarning(f"Region with name {region_name} does not exist. "+
                              f"Existing regions are: {self.get_region_names()}")
        
    def has_region(self, region_name):
        # return region_name in self.get_region_names()
        return hasattr(self, region_name)

         
        