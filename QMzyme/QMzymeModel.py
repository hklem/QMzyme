###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@lsu.edu
###############################################################################

import os
import pickle
from QMzyme.CalculateModel import CalculateModel
import QMzyme.MDAnalysisWrapper as MDAwrapper
from QMzyme.data import protein_residues, residue_charges
import MDAnalysis
from QMzyme.RegionBuilder import RegionBuilder

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
        """
        Add a QMzymeRegion to the model.

        Registers the region as an attribute on the model instance and appends it 
        to the tracked regions list.

        :param region: The QMzymeRegion instance to add to the model.
        :type region: :class:`~QMzyme.QMzymeRegion.QMzymeRegion`
        """
        if region.n_atoms == 0:
            raise UserWarning(f"Region contains no atoms and will not be created.")
        if hasattr(self, region.name):
            raise UserWarning(f"Region with name {region.name} already exists in QMzymeModel {self.name}."+
                               "Please use a different region name or remove the existing region via remove_region({region.name}).")
        setattr(self, region.name, region)
        self.regions.append(region)

    def get_region_names(self):
        """
        Get the names of all regions currently registered to the model.

        :returns: A list of region names.
        :rtype: list of str
        """
        return [r.name for r in self.regions]
    
    def get_region(self, region_name=None):
        """
        Retrieve a specific QMzymeRegion by its name.

        :param region_name: The name of the region to fetch.
        :type region_name: str, default=None

        :returns: The requested QMzyme region object.
        :rtype: :class:`~QMzyme.QMzymeRegion.QMzymeRegion`
        """
        try:
            return getattr(self,region_name)
        except:
            raise UserWarning(f"Region with name {region_name} does not exist. "+
                              f"Existing regions are: {self.get_region_names()}")
        
    def has_region(self, region_name):
        """
        Check if a region name is already registered as an attribute on the model.

        :param region_name: The name of the region to look for.
        :type region_name: str

        :returns: True if the model has an attribute matching region_name, False otherwise.
        :rtype: bool
        """
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
    
    def import_region(self, region_file, name):
        """
        Loads a new structure file to the existing universe, shifts its
        atom IDs to prevent overlapping with the current model's universe.

        :param lig_file: Path to the structure file (e.g., PDB).
        :type lig_file: str
        :param name: Name of the returned region.
        :type name: str
        :returns: The built replacement ligand region.
        :rtype: :class:`~QMzyme.QMzymeRegion.QMzymeRegion`
        """
        if name is None:
            raise UserWarning("Please specify name of the new region.")

        # Load the replacement ligand into an MDAnalysis universe
        add_universe = MDAnalysis.Universe(region_file)

        # Shift atom IDs using this region's parent universe max ID to prevent conflict
        if hasattr(self, '_universe') and self._universe is not None:
            max_protein_id = self._universe.atoms.ids.max()
            add_universe.atoms.ids += max_protein_id
        elif hasattr(self, 'universe') and self.universe is not None:
            max_protein_id = self.universe.atoms.ids.max()
            add_universe.atoms.ids += max_protein_id

        # Create the new QMzymeRegion using the RegionBuilder
        region = RegionBuilder(
            name=name, 
            atom_group=add_universe.select_atoms("all")
        ).get_region()

        self.add_region(region)

        return region

    def pymol_visualize(self, filename:str=None, model_surface:bool=True, output_dir=None):
        """
        Creates a QMzymeModel_visualize.py script that you can load into PyMol.

        :param filename: Name of PyMol .py file. If not specified, the name 
            attribute of the QMzymeModel will be used.
        :type filename: str, optional
        
        :model_surface: Turning this into False will reduce the GPU load for PyMOL, sometimes preventing GL error.
        :type model_surface: boolean, optional

        :output_dir: Name of the output directory. If not specified. the python sript and pdb files generated from
            pymol_visualize() will be placed at the base directory.
        :type output_folder: str, optional
        """
        
        original_dir = None

        if output_dir is not None:
            os.makedirs(output_dir, exist_ok=True)
            original_dir = os.getcwd()
            os.chdir(output_dir)
        
        try:
            lines = ''
            lines += f"cmd.bg_color('white')\n"
            starting_structure = self.name
            self.universe.atoms.write(f"{self.name}_universe.pdb")
            file = os.path.abspath(f'{self.name}_universe.pdb')
            lines += f"cmd.load(r'{file}', '{self.name}')\n"
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
                lines += f"cmd.load(r'{file}', '{region.name}')\n"
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
                lines += f"cmd.load(r'{file}', '{region.name}')\n"
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
                if model_surface is True:
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

        finally:
            if original_dir is not None:
                os.chdir(original_dir)

    def print_overview(self):
        """
        Prints a formatted overview of the model and its regions, including 
        atom/residue counts, designated methods, and creation parameters.
        """
        print("-" * 29)
        print(f"Model Overview: {self.name} ")
        print("-" * 29)
        print(f"  - total atoms: {self.universe.atoms.n_atoms}")
        print(f"  - total residues: {self.universe.residues.n_residues}")
        print(f"  - total regions: {len(self.regions)}")
        print("-" * 29)
        print("Region Overview")
        print("-" * 29)

        for region in self.regions:
            print(f"Region Name: {region.name}")
            
            # Use built-in QMzymeRegion attributes
            print(f"  - atoms: {region.n_atoms}")
            print(f"  - residues: {region.n_residues}")
            
            # Safely retrieve the method, defaulting to None if it doesn't exist
            method = getattr(region, 'method', None)
            print(f"  - method: {method}")
            
            # If it does not have any creation_params, state that information is not avaialable!
            try:
                attr = region.creation_attr
                for key, value in attr.items():
                    if key in ['total_atoms', 'total_residues']:
                        continue
                    if key == 'selection_scheme' and 'parent_region' in attr:
                        continue
                    print(f"  - {key}: {value}")
            except AttributeError:
                print("  - selection_scheme: information not available")

            # Derive truncation_scheme from residue-level truncation_params
            residue_schemes = []
            for scheme in region._residue_truncation_attr.values():
                if scheme is not None and scheme not in residue_schemes:
                    residue_schemes.append(scheme)
            if residue_schemes:
                print(f"  - truncation_scheme: {residue_schemes}")

            print("-" * 29)

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
