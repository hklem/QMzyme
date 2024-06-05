###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
Product of the ModelBuilder class.
"""

import os
from QMzyme.CalculateModel import CalculateModel

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
             raise UserWarning(f"Region with name {region.name} already exists in QMzymeModel {self.name}."+
                               "Please use a different region name or remove the existing region via remove_region({region.name}).")
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
    
    def pymol_visualize(self, filename:str=None):
        """
        Creates a QMzymeModel_visualize.py script that you can load into PyMol.
        """
        lines = ''
        lines += f"cmd.bg_color('white')\n"
        starting_structure = self.name
        self.universe.atoms.write(f"{self.name}_universe.pdb")
        file = os.path.abspath(f'{self.name}_universe.pdb')
        lines += f"cmd.load('{file}', '{self.name}')\n"
        #lines += f"cmd.color('gray70', self.name)\n"
        lines += f"cmd.set('surface_color', 'gray')\n"
        lines += f"cmd.set('transparency', 0.75)\n"
        lines += f"cmd.zoom('visible')\n"
        lines += f"cmd.orient('visible')\n"
        lines += f"cmd.scene('Starting Structure', 'store')\n"
        lines += f"cmd.hide('everything', '{self.name}')\n"

        for region in self.regions:
            region.write(f'{region.name}.pdb')
            file = os.path.abspath(f'{region.name}.pdb')
            lines += f"cmd.load('{file}', '{region.name}')\n"
            lines += f"cmd.hide('cartoon', '{region.name}')\n"
            lines += f"cmd.show_as('sticks', '{region.name}')\n"
            lines += f"cmd.zoom('visible')\n"
            lines += f"cmd.orient('visible')\n"
            lines += f"cmd.scene('{region.name}', 'store')\n"
            lines += f"cmd.hide('everything', '{region.name}')\n"


        if CalculateModel.calc_type != None:
            region = CalculateModel.calculation[CalculateModel.calc_type]
            region.write(f'{region.name}.pdb')
            file = os.path.abspath(f'{region.name}.pdb')
            lines += f"cmd.load('{file}', '{region.name}')\n"
            lines += f"cmd.hide('cartoon', '{region.name}')\n"
            lines += f"cmd.color('gray85', '{region.name} and elem c')\n"
            lines += f"cmd.color('oxygen','{region.name} and elem o')\n"
            lines += f"cmd.color('slate', '{region.name} and elem n')\n"
            lines += f"cmd.color('gray98', '{region.name} and elem h')\n"
            lines += f"cmd.color('sulfur', '{region.name} and elem s')\n"
            lines += f"cmd.show_as('sticks', '{region.name} and segid QM')\n"
            lines += f"cmd.show_as('lines', '{region.name} and (not segid QM)')\n"
            fixed = [str(i+1) for i, atom in enumerate(region.atoms) if atom.is_fixed]
            fixed_sel = f"id {'+'.join(fixed)}"
            if len(fixed) > 0:
                lines += f"cmd.create('fixed_atoms', '{region.name} and {fixed_sel}')\n"
                lines += f"cmd.hide('cartoon', 'fixed_atoms')\n"
                lines += f"cmd.set('sphere_scale', 0.15, 'fixed_atoms')\n"
                lines += f"cmd.set('sphere_color', 'black', 'fixed_atoms')\n"
                #lines += f"cmd.set('sphere_transparency', 0.7, 'fixed_atoms')\n"
                lines += f"cmd.show_as('spheres', 'fixed_atoms')\n"
            #lines += f"cmd.select('residue_labels', '{region.name}')\n"
            lines += f"cmd.create('residue_labels', '{region.name}')\n"
            lines += f"cmd.hide('everything', 'residue_labels')\n"
            lines += f"cmd.set('label_size', 14)\n"
            lines += f"cmd.label('n. ha and residue_labels', 'resn+resi')\n"

            lines += f"cmd.zoom('visible')\n"
            lines += f"cmd.create('model_surface', '{region.name}')\n"
            lines += f"cmd.show_as('surface', 'model_surface')\n"
            lines += f"cmd.orient('visible')\n"
            lines += f"cmd.scene('{region.name}', 'store')\n"

            lines += f"cmd.set('cartoon_transparency', 0.6)\n"
            #lines += f"cmd.show('surface', '{region.name}')\n"
            lines += f"cmd.show('cartoon', '{self.name}')\n"
            lines += f"cmd.zoom('visible')\n"
            lines += f"cmd.orient('visible')\n"

        # if CalculateModel.calculation != {}:
        #     for calc, region in CalculateModel.calculation.items():
        #         region.write(f'{region.name}_{calc}.pdb')
        #         file = os.path.abspath(f'{region.name}_{calc}.pdb')
        #         lines += f"cmd.load('{file}', '{calc}')\n"
        #         if 'QM' in calc and calc != "QM2":
        #             lines += f"cmd.show_as('sticks', '{calc}')\n"
        #         else:
        #             lines += f"cmd.show_as('lines', '{calc}')\n"
        #         #lines += f"cmd.color('cyan','Catalytic_Center')\n"
        #         lines += f"cmd.hide('everything', '{self.name}')\n"
        #         lines += f"cmd.zoom('visible')\n"
        #         lines += f"cmd.orient('visible')\n"
        #         lines += f"cmd.scene('{calc}', 'store')\n"

        # else:
        #     for i,region in enumerate(ordered_regions):
        #         region.write(f'{region.name}.pdb')
        #         file = os.path.abspath(f'{region.name}.pdb')
        #         lines += f"cmd.load('{file}', '{region.name}')\n"
        #         lines += f"cmd.show_as('sticks', '{region.name}')\n"
        #         #lines += f"cmd.color('cyan','Catalytic_Center')\n"
        #         lines += f"cmd.hide('everything', '{self.name}')\n"
        #         lines += f"cmd.zoom('visible')\n"
        #         lines += f"cmd.orient('visible')\n"
        #         lines += f"cmd.scene('{region.name}', 'store')\n"

        if filename == None:
            filename = f'QMzymeModel_{self.name}_visualize.py'
        elif not filename.endswith('.py'):
            filename = filename+'.py'
        with open (filename, 'w+') as f:
            f.write(lines)

    def remove_region(self, region_name):
        """
        Method to remove a region from the QMzymeModel.
        :param region_name: Name of the region to be removed.
        :type region_name: str, required 
        """
        #del self.regions[region_idx]
        delattr(self, region_name)
        n_regions = len(self.regions)
        for i in range(n_regions):
            if self.regions[i].name == region_name:
                del self.regions[i]
                break