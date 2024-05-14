###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
Product of the ModelBuilder class.
"""

import os

class QMzymeModel:
    # def __init__(self, name, starting_structure, regions=[]):
    def __init__(self, name, universe):
        self.name = name
        self.universe = universe
        self.filename = universe.filename
        self.regions = []

    def __repr__(self):
        return f"<QMzymeModel {self.name} built from {self.universe} contains {self.n_regions} region(s)>"

    @property
    def n_regions(self):
        return len(self.regions)
    
    def add_region(self, region):
         if hasattr(self, region.name):
             raise UserWarning(f"Region with name {region.name} already exists in QMzymeModel {self.name}. Please use a different region name or remove the existing region via remove_region({region.name}).")
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
    
    def pymol_visualize(self, region_colors=['cyan', 'orange', 'purple', 'green']):
        """
        Creates a QMzymeModel_visualize.py script that you can load into PyMol.
        """
        lines = ''
        starting_structure = self.name
        self.universe.atoms.write(f"{self.name}.pdb")
        file = os.path.abspath(f'{self.name}.pdb')
        lines += f"cmd.load('{file}', '{self.name}')\n"
        #lines += f"cmd.color('gray70', self.name)\n"
        lines += f"cmd.orient()\n"
        lines += f"cmd.scene('Starting Structure', 'store')\n"
        n_atoms = [region.n_atoms for region in self.regions]
        ordered_regions = [x for _, x in sorted(zip(n_atoms, self.regions))]
        for i,region in enumerate(ordered_regions):
            region.write(f'{region.name}.pdb')
            file = os.path.abspath(f'{region.name}.pdb')
            lines += f"cmd.load('{file}', '{region.name}')\n"
            lines += f"cmd.show_as('sticks', '{region.name}')\n"
            #lines += f"cmd.color('cyan','Catalytic_Center')\n"
            lines += f"cmd.hide('everything', '{self.name}')\n"
            lines += f"cmd.zoom('visible')\n"
            lines += f"cmd.orient('visible')\n"
            lines += f"cmd.scene('{region.name}', 'store')\n"

        with open ('QMzymeModel_visualize.py', 'w+') as f:
            f.write(lines)

    def remove_region(self, region_name):
        """
        Method to remove a region from the QMzymeModel.
        :param region_name: Name of the region to be removed.
        :type region_name: str, required 
        """
        del self.regions[region_name]
        delattr(self, region_name)





            



         
        