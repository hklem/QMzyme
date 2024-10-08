###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

import os
import pickle
from QMzyme.CalculateModel import CalculateModel
import QMzyme.MDAnalysisWrapper as MDAwrapper
from QMzyme.data import protein_residues, residue_charges

class QMzymeModel:
    """
    Base class for :class:`QMzyme.GenerateModel`. Contains methods to create and 
    modify a ``QMzymeModel`` instance. A QMzymeModel an be instantiated with an 
    MDAnalysis Universe directly, or any combination of parameters that 
    MDAnalysis.core.universe.Universe accepts to create a Universe i.e., 
    (example.prmtop, example.dcd, dt=5).
    See https://userguide.mdanalysis.org/stable/universe.html for details.

    :param name: Name to give to the QMzymeModel. This is used for default file naming 
        purposes throughout the QMzyme package. If not provided, it will default to
        the base name of the universe filename attribute. 
    :type name: str, optional
    :param universe: MDAnalysis Universe object.
    :type universe: `MDAnalysis.Universe <https://userguide.mdanalysis.org/stable/universe.html>`_, default=None
    

    :param name: Name of QMzymeModel
    :type name: str, default=None
    :param universe: MDAnalysis Universe object
    :type universe: `MDAnalysis.Universe <https://userguide.mdanalysis.org/stable/universe.html>`_, default=None
    :param frame: If trajectory was provided, specify a frame to base coordinates on
    :type frame: int, default=0
    :param pickle_file: Provide name/path+file of previously pickled QMzymeModel object to inialize
    :type pickle_file: str, default=None
    """
    def __init__(self, *args, name, universe, select_atoms='all', frame=0, pickle_file=None, **kwargs):
        if pickle_file is not None:
            with open(pickle_file, 'rb') as f:
                self.__dict__ = pickle.load(f).__dict__
            return
        elif universe is None:
            universe = MDAwrapper.init_universe(*args, frame=frame, **kwargs)
        self.universe = MDAwrapper.universe_selection(universe, select_atoms)
        self.select_atoms = select_atoms
        if name is None:
            name = os.path.basename(universe.filename).split('.')[0]
        self.name = name
        self.frame = frame
        self.filename = universe.filename
        self.regions = []

        if not hasattr(self.universe.atoms, "charges"):
            print("\nCharge information not present. QMzyme will try to guess "+
                  "region charges based on residue names consistent with AMBER naming "+
                  "conventions (i.e., aspartate: ASP --> Charge: -1, aspartic acid: ASH --> Charge: 0.). "+
                  "See QMzyme.data.residue_charges for the full set.")
            unk_res = []
            for res in self.universe.residues:
                if res.resname not in residue_charges:
                    if unk_res == []:
                        print("\n\tNonconventional Residues Found")
                        print("\t------------------------------")
                    if res.resname not in unk_res:
                        unk_res.append(res.resname)
                        print(f"\t{res.resname} --> Charge: UNK, defaulting to 0")
            if unk_res != []:
                print("\nYou can update charge information for nonconventional residues by running "+
                      "\n\t>>>QMzyme.data.residue_charges.update("+"{"+"'3LETTER_RESNAME':INTEGER_CHARGE}). "+
                      "\nNote your changes will not be stored after you exit your session. "+
                      "It is recommended to only alter the residue_charges dictionary. "+
                      "If you alter the protein_residues dictionary instead that could cause "+
                      "unintended bugs in other modules (TruncationSchemes).\n")

    def __repr__(self):
        return f"<QMzymeModel {self.name} built from {self.universe} contains {self.n_regions} region(s)>"

    @property
    def n_regions(self):
        return len(self.regions)
    
    def add_region(self, region):
        if region.n_atoms == 0:
            raise UserWarning(f"Region contains no atoms and will not be created.")
        if hasattr(self, region.name):
            raise UserWarning(f"Region with name {region.name} already exists in QMzymeModel {self.name}."+
                               "Please use a different region name or remove the existing region via remove_region({region.name}).")
        setattr(self, region.name, region)
        self.regions.append(region)

    def get_region_names(self):
        return [r.name for r in self.regions]
    
    def get_region(self, region_name=None):
        try:
            return getattr(self,region_name)
        except:
            raise UserWarning(f"Region with name {region_name} does not exist. "+
                              f"Existing regions are: {self.get_region_names()}")
        
    def has_region(self, region_name):
        # return region_name in self.get_region_names()
        return hasattr(self, region_name)

    def remove_region(self, region_name):
        """
        Removes a QMzymeRegion from the QMzymeModel.
        
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
    
    def pymol_visualize(self, filename:str=None):
        """
        Creates a QMzymeModel_visualize.py script that you can load into PyMol.

        :param filename: Name of PyMol .py file. If not specified, the name 
            attribute of the QMzymeModel will be used.
        :type filename: str, optional
        
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

        if filename == None:
            filename = f'QMzymeModel_{self.name}_visualize.py'
        elif not filename.endswith('.py'):
            filename = filename+'.py'
        with open (filename, 'w+') as f:
            f.write(lines)

    def store_pickle(self, filename=None):
        """
        The pickle file will by default be named after the QMzymeModel.name 
        attribute, which by default is the base name of the file originally 
        used to initialize the model. You can also specify a filename by 
        passing the argument 'filename'.

        :param filename: Name of the pickle file generated. If no extension 
        is provided the '.pkl' extension will be added to the str
        :type filename: str
        """
        if filename == None:
            filename = self.name+'.pkl'
        if len(filename.split('.')) < 2:
            filename += '.pkl'
        with open(filename, 'wb') as f:
            pickle.dump(self, f)
