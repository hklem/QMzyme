###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@lsu.edu
###############################################################################

"""
Module containing functions to define a QMzymeRegion based on some logic/workflow.
"""

import QMzyme.MDAnalysisWrapper as MDAwrapper
import QMzyme
from QMzyme.RegionBuilder import RegionBuilder
from QMzyme.QMzymeRegion import QMzymeRegion
from QMzyme.QMzymeModel import QMzymeModel
from QMzyme.TruncationSchemes import BetaCarbon
from QMzyme.CalculateModel import QM_Method
from typing import Type
import pickle
import pandas as pd

import abc

class SelectionScheme(abc.ABC):
    """
    SelectionScheme is the abstract base class used to prescribe concrete selection scheme
    sub classes. Fork QMzyme to build your own concrete selection scheme class and submit
    a Pull Request (PR) on github to have your scheme added to QMzyme! You will need to create
    comprehensive tests before the PR is accepted. See the 
    `documentation on contributing <https://qmzyme.readthedocs.io/en/latest/Contributing/index.html>`_ 
    for more information.

    Below is a template that you can use to implement your own selection scheme class:

    class InformativeName(SelectionScheme):
        This is an example docstring.

        Include a detailed description of how the scheme works in the class level doc string. Also
        include any parameters the `__init__()` method of your class accepts.

        :Parameters:
            :model: (:class:`~QMzyme.QMzymeModel.QMzymeModel`) QMzymeModel to provide starting structure that selection will be performed on. When
                using the main :class:`~QMzyme.GenerateModel.GenerateModel` class, the QMzyme model is automatically 
                passed as an argument to the selection scheme. It is recommended you use the Universe (universe attribute)
                representing the starting structure to perform the selection on. 
            :param name: (str, required) Name of the region generated.

        The return should always be a QMzyme region.

        :Returns:
            :class:`~QMzyme.QMzymeRegion.QMzymeRegion`
        
        .. note::

            Include any notes you want users to be aware of.

    """
    def __init__(self, model, name):
        """
        Assign any key word arguments as attributes to self. Then in your
        `select_atoms()` method you can pull any necessary args from
        self attributes, instead of relying on passing them.

        Every concrete scheme `__init__()` method should include this line 
        at the very end:

        .. code:: python

            super().__init__(model, name)

        This will automatically run your `select_atoms()` method and return the resulting region.
        """
        self.name = name
        self.model: QMzymeModel = model
        self.region: QMzymeRegion
        self.select_atoms()
        self.reference()
        if self.reference is not None:
            print(f"Use of this selection scheme requires citing the following reference(s): \n \t{self.reference}")
        self.return_region()

    @abc.abstractmethod
    def select_atoms(self):
        """
        Write your code to perform the selection. 

        At the end of your code you should set `self.region = {region}`. 

        The product of your selection scheme needs to be a QMzymeRegion 
        in order for it to work with `GenerateModel().set_region()`.

        This method is automatically called in the ``super().__init__(model, name)``
        line of your `__init__()` method.
        """
        ...
    
    def method_name(self):
        """
        You can add whatever other methods you want in your class, but 
        you should call those methods as necessary in `__init__()` otherwise
        your scheme will be automated in `GenerateModel.set_region()`
        """
        pass

    def return_region(self):
        """
        This method belongs to the base class and is automatically called in 
        the ``super().__init__(model, name)`` line of your `__init__()` method. All
        you have to do is make sure you have created a class attribute called `region`.
        """
        self.region.rename(self.name)
        return self.region

    @abc.abstractmethod
    def reference(self):
        """
        This method needs to be included in your class. All it should do is create an
        attribute called `reference` that provides a citable reference of the scheme, to 
        give credit where credit is due. The reference will be automatically printed when
        the class is instantiated. This is taken care of in the the ``super().__init__(model, name)`` 
        line of your `__init__()` method.
        
        Example:

        .. code:: python

            self.reference = "1. Alegre‐Requena, J. V., Sowndarya S. V., S., Pérez‐Soto, R., Alturaifi, T. M. & Paton, R. S. AQME: Automated quantum mechanical environments for researchers and educators. WIREs Comput Mol Sci 13, e1663 (2023)."

        In some cases, there might not be a direct reference (see DistanceCutoff class), but 
        there might be relevant work a user might be interested in. Please only refer to the 
        work of interest in the class doc string, not in the reference method. 
        
        If there are no references, please only include the line:

        .. code:: python

            self.reference = None

        """
        ...


class DistanceCutoff(SelectionScheme):
    """
    The DistanceCutoff class performs a selection simply based on the distance of 
    atoms from a pre-defined catalytic_center region. Users must first call 
    ``GenerateModel().set_catalytic_center(kwargs)`` in order to then run DistanceCutoff via
    ``GenerateModel().set_region(selection=DistanceCutoff, name={str}, cutoff={int})``. 
    
    This scheme is known to require rather large QM regions to achieve agreement 
    with experiment (Ex., Kulik HJ, Zhang J, Klinman JP, Martínez TJ. How Large 
    Should the QM Region Be in QM/MM Calculations? The Case of Catechol 
    O-Methyltransferase. J Phys Chem B. 2016 Nov 10;120(44):11381-11394. 
    doi: 10.1021/acs.jpcb.6b07814.).


    :param model: QMzymeModel to provide starting structure that selection 
        will be performed on.
    :type model: :class:`~QMzyme.QMzymeModel.QMzymeModel`, required.

    :param name: Name of the region generated.
    :type name: str, required.

    :param cutoff: Numerical value to define cutoff.
    :type cutoff: float, required.

    :param include_whole_residues: Informs code whether or not to only include
        atoms within the cutoff, or include whole residues if they have at least one
        atom within the cutoff. 
    :type include_whole_residues: bool, default=True.


    :returns: :class:`~QMzyme.QMzymeRegion.QMzymeRegion`

    .. note::

        Users are encouraged to evaluate the resulting region. There may be situations where
        a charged residue is within the cutoff distance, however, its charge partner is not. 
        Such situations can drastically alter the chemistry of the model! Maybe someone could
        write up a less generic distance based selection scheme that would take such situations
        into consideration. Or modify the current class to include an argument 
        `include_charge_partners=True`. 
        
    """
    def __init__(self, model, name, cutoff, include_whole_residues=True):
        """
        """
        if name is None:
            name = f'cutoff_{cutoff}'
        self.cutoff = cutoff
        self.include_whole_residues = include_whole_residues
        super().__init__(model, name)

    def select_atoms(self):
        """
        """
        if not self.model.has_region('catalytic_center'):
            raise UserWarning("You must first define a catalytic_center. See method `set_catalytic_center()`.")
        
        neighbors = MDAwrapper.get_neighbors(
        self.model.universe.select_atoms('all'),
        self.model.get_region('catalytic_center')._atom_group, self.cutoff)
    
        neighbors_byres = neighbors.residues.sorted_unique.atoms

        # Convert MDAnalysis AtomGroup to QMzymeRegion
        qmz_neighbors = RegionBuilder(name=f"cutoff_{self.cutoff}", atom_group=neighbors).get_region()
        qmz_neighbors_byres = RegionBuilder(name=f"cutoff_{self.cutoff}", atom_group=neighbors_byres).get_region()

        for atom in qmz_neighbors_byres.atoms:
            if atom.id in qmz_neighbors.ids:
                atom.set_neighbor(True)
            else:
                atom.set_neighbor(False)

        region = qmz_neighbors_byres
        if self.include_whole_residues is False:
            region = qmz_neighbors

        setattr(region, 'catalytic_center', self.model.get_region('catalytic_center'))
        self.region = region

    def reference(self):
        """
        """
        self.reference = None


class ChargeShiftAnalysis(SelectionScheme):
    """
    The ChargeShiftAnalysis Selection Scheme is based on the CSA approach developed by 
    Prof. Heather Kulik (see (1) Kulik, Heather J.; Zhang, Jianyu; Klinman, Judith P.; Martinez, Todd J.(2016) 
    How Large Should the QM Region Be in QM/MM Calculations? The Case of Catechol O-Methyltransferase. Journal of Physical 
    Chemistry B, 120(44). and (2) Karelina, M., & Kulik, H. J. (2017). Systematic quantum 
    mechanical region determination in QM/MM simulation. Journal of chemical theory 
    and computation, 13(2)). In this scheme, two initial large single-point QM calculations
    are performed with atomic population analysis with and without the substrate present (Holo and Apo, respectively). 
    The rank, or importance of including a residue in the QM region is determined by how much the residue-summed partial
    charges shift in the presence of the substrate. Residues above a threshold charge shift value will be included in 
    the returned QMzymeRegion.
    
    This class requires two steps, and therefore two separate calls. The first call will create the Apo and Holo QM calculation 
    input files. The user will then need to run those calculations. After the calculations are complete, the user can reload 
    their saved serialized QMzyme.GenerateModel object (".pkl" format) and then call the ChargeShiftAnalysis class, this time 
    uploading the Apo and Holo QM output files and specifying a charge shift cutoff to return the QMzymeRegion. 
    
    :param model: QMzymeModel to provide starting structure that selection will be performed on.
    :type model: :class:`~QMzyme.QMzymeModel.QMzymeModel`, required.

    :param name: Name of the region generated.
    :type name: str, required.

    :param method: Method used to write QM input file.
    :type method: :class:`~QMzyme.CalculateModel.QM_Method`, required.

    :param min_atoms: Minimum number of atoms for first QM input file.
    :type min_atoms: int, default=900.

    :param max_atoms: Maximum number of atoms for first QM input file.
    :type max_atoms: int, default=1000.

    :param include_whole_residues: Informs code whether or not to only include
        atoms within the cutoff, or include whole residues if they have at least one
        atom within the cutoff. 
    :type include_whole_residues: bool, default=True.
    
    :param alanine_mutation: Selects specific residues to mutate into alanine.
        The input should be the format for MDAnalysis residue selection.
    :type alanine_mutation: str.
    
    :param memory: Memory for QM input file.
    :type memory: str, default=None.

    :param nprocs: Number of processors for QM input file.
    :type nprocs: int, default=None.

    :param CSA_elbow: Chooses elbow method for CSA selection.
    :type CSA_elbow: bool, default=False.

    :param holo_output_files: Location of holo output file.
    :type holo_output_files: str, default=None.

    :param apo_output_files: Location of apo output file.
    :type apo_output_files: str, default=None.
    
    :param pop: Method for population analysis.
    :type pop: str, default="hirshfeld".

    :param charge_threshold: Partial charge threshold for CSA selection.
    :type charge_threshold: float, required.

    :returns: :class:`~QMzyme.QMzymeRegion.QMzymeRegion`
    
    
    :Usage:
        Users must first set a catalytic center and define a QM_Method:
        
        .. code-block:: python
                import QMzyme
                model = QMzyme.GenerateModel("filename.pdb")
                model.set_catalytic_center(kwargs)
                qm_method = QMzyme.CalculateModel.QM_Method(kwargs)
                
        Then they can call the ChargeShiftAnalysis Selection Scheme:
        
        .. code-block:: python
                model.set_region(selection=QMzyme.SelectionSchemes.ChargeShiftAnalysis, method=qm_method, name={str})
            
        This will create the Holo and Apo input files. The user will then run these QM calculations. 
        Once the calculations have completed, the user can start back where they left off by loading in the saved .pkl file 
        containing the QMzyme.GenerateModel object they instantiated in the previous step:
        
        .. code-block:: python
                import QMzyme
                model = QMzyme.GenerateModel("model.pkl") 
                model.set_region(selection=QMzyme.SelectionSchemes.ChargeShiftAnalysis, holo_output_files={str}, apo_output_files={str}, pop={str}, charge_threshold={float})
                


    .. note::

        Users are encouraged to evaluate the resulting region. There may be situations where a
        charged residue is included because it falls within the distance and atom number range,
        but its charge partner is excluded. Such imbalances can drastically alter the chemistry
        of your model! Additionally, please note that the minimum and maximum atom limits set
        within the class apply to the number of atoms before truncation; therefore, the final
        atom count might not fall within that initial range. It is highly recommended to verify
        the structure of your initial input file to ensure it is of high quality.

        In addition, it is crucial to set the method as QM_Method. This setting is strictly required
        for performing the initial population analysis and defining the QMzyme region using charge
        shift analysis. Any other method selection will return an error message.

        For the second part of the ChargeShiftAnalysis class, you must have the .pkl file generated from the
        first run. Without this file, the ChargeShiftAnalysis class cannot perform the charge shift analysis.
        
    """

    def __init__(self, model, name, method: QM_Method, min_atoms=900, max_atoms=1000, include_whole_residues=True, alanine_mutation=[], memory=None, nprocs=None, holo_output_files=None, apo_output_files=None, pop="hirshfeld", charge_threshold=float):

        self.min_atoms = min_atoms
        self.max_atoms = max_atoms
        self.model = model
        self.cutoff = 5
        self.method = method
        self.include_whole_residues = include_whole_residues
        self.alanine_mutation = alanine_mutation
        self.holo_output_file = holo_output_files
        self.apo_output_file  = apo_output_files
        self.memory = memory
        self.nprocs = nprocs
        self.pop = pop
        self.charge_threshold = charge_threshold

        # Store the parameters that you want the regions to "remember"
        self.holo_selection_attr = {
            'selection_scheme': self.__class__,
            'min_atoms': min_atoms,
            'max_atoms': max_atoms
        }

        self.apo_selection_attr = {
            'selection_scheme': self.__class__,
            'min_atoms': min_atoms,
            'max_atoms': max_atoms,
            'alanine_mutation': alanine_mutation
        }

        # Function that selects catalytic center atoms that are not protein or water.
        self.select_cat_residues()
        catalytic_center = self.model.get_region('catalytic_center')

        # Raise user warning if catalytic_center is missing.
        if catalytic_center is None:
            raise UserWarning("You must first define a catalytic_center. See method `set_catalytic_center()`.")
        
        # Raises user warning if catalytic center does not contain non-protein and non-water atoms for apo form.
        if self.cat_center_atoms is None or len(self.cat_center_atoms) == 0:
            raise UserWarning("The catalytic center must contain non-protein and non-water atoms.")
        
        if method is None or not isinstance(method, QMzyme.QM_Method):
            raise UserWarning("You must set method as QM_Method for partial charge calculation.")

        # Checking if there is any QM_output_files provided by the user
        if holo_output_files is not None and apo_output_files is not None:
            if charge_threshold is not float:
                print("Calculating partial charge differences from provided QM output files.")
                name = "CSA_cutoff_region"
                super().__init__(model, name)
                QMzyme.CalculateModel._reset()
                self.method.assign_to_region(self.region)

            elif charge_threshold is None or not isinstance(charge_threshold, (float)):
                raise UserWarning("You must set charge_threshold as float or int for partial charge calculation.")

        else:
            print("Creating Holo and Apo input files for partial charge calculations.")
            self.create_apo_holo_regions()
            self.name = "CSA_initial_selection_conditions"

            self.full_method = self.method.__class__.__new__(self.method.__class__)
            self.full_method.__dict__.update(self.method.__dict__)
            self.full_method.starting_file = f"{self.model.name}"

            self.full_method.assign_to_region(self.region)
            self.reference()
            if self.reference is not None:
                print(f"Use of this selection scheme requires citing the following reference(s): \n \t{self.reference}")
            self.region.reset_creation_attr()

    def select_cat_residues(self):
        """
        Docstring for select_cat_residues
        
        :param self: Description
        """
        catalytic_center = self.model.get_region('catalytic_center')
        u = self.model.universe
        indices = [atom.index for atom in catalytic_center.atoms]
        cat_atoms = u.atoms[indices]

        # Add warning when there is presence of protein atoms
        cat_protein_atoms = cat_atoms.select_atoms("protein")
        if len(cat_protein_atoms) > 0:
            print(f"The catalytic center contains protein amino acid residues {cat_protein_atoms}")

        cat_center_atoms = cat_atoms.select_atoms("not protein and not water")
        self.cat_center_atoms = cat_center_atoms

    def create_apo_holo_regions(self):
        """
        """
        # Counter to avoid infinite loops
        counter = 0
        
        # Examine what the total number of atom is!
        all_atoms = self.model.universe.select_atoms('all')

        # Check if the total number of atoms within the system is less than max atom limit
        if len(all_atoms.atoms) <= self.max_atoms:
            holo_region = RegionBuilder(name="CSA_holo", atom_group=all_atoms).get_region()
            CSA_initial_region = RegionBuilder(name="CSA_initial", atom_group=all_atoms).get_region()
            # This sets a region attribute to a QMzymeRegion object. So a potential subregion...? But it does not change the catalytic center region of the original model.
            setattr(holo_region, 'catalytic_center', self.model.get_region('catalytic_center'))

            # Assign method to region and truncate
            self.method.assign_to_region(holo_region)
            self.model.truncate()

        # If # atoms within the whole enzyme is more than max atoms, need to truncate
        else:

            # Repeat this loop until the correct truncation or region error
            while True:
                # Safety counter to avoid infinite loops
                counter += 1
                if counter > 20:
                    raise UserWarning("Exceeded maximum number of iterations to find suitable cutoff.")

                # Setting neighbors using MDAwrapper
                neighbors = MDAwrapper.get_neighbors(
                self.model.universe.select_atoms('all'),
                self.model.get_region('catalytic_center')._atom_group, self.cutoff)
            
                neighbors_byres = neighbors.residues.sorted_unique.atoms

                # Convert MDAnalysis AtomGroup to QMzymeRegion
                qmz_neighbors = RegionBuilder(name=f"cutoff_{self.cutoff}", atom_group=neighbors).get_region()
                qmz_neighbors_byres = RegionBuilder(name=f"cutoff_{self.cutoff}", atom_group=neighbors_byres).get_region()

                # Mark neighbor atoms
                for atom in qmz_neighbors_byres.atoms:
                    if atom.id in qmz_neighbors.ids:
                        atom.set_neighbor(True)
                    else:
                        atom.set_neighbor(False)

                # Define region based on whole residues or not
                holo_region = RegionBuilder(name="CSA_holo", atom_group=neighbors_byres).get_region()

                # Check if the number of atoms is within the min and max atoms set
                if self.min_atoms <= len(holo_region.atoms) <= self.max_atoms:
                    CSA_initial_region = RegionBuilder(name="CSA_initial", atom_group=neighbors_byres).get_region()
                    # This sets a region attribute to a QMzymeRegion object. So a potential subregion...?
                    # But it does not change the catalytic center region of the original model
                    setattr(holo_region, 'catalytic_center', self.model.get_region('catalytic_center'))

                    # Create a fresh instance of the same class without calling __init__
                    holo_method = self.method.__class__.__new__(self.method.__class__)

                    # Copy the attributes over (this is a shallow clone)
                    holo_method.__dict__.update(self.method.__dict__)

                    # Assign the unique clone to the region
                    holo_method.assign_to_region(holo_region)

                    # Now, save region and creation_param for the pickle file
                    self.model.set_region(holo_region)
                    holo_region.set_creation_attr(self.holo_selection_attr)

                    self.model.truncate()

                    break

                # If not, increase or decrease the cutoff accordingly
                if self.min_atoms > len(holo_region.atoms):
                    self.cutoff += 1

                elif self.max_atoms < len(holo_region.atoms):
                    self.cutoff -= 0.25

        # Write holo input file
        kwargs = {}
        if self.memory is not None:
            kwargs["memory"] = self.memory
        if self.nprocs is not None:
            kwargs["nprocs"] = self.nprocs
        kwargs["reset_calculation"] = True
        self.model.write_input(**kwargs)

        # Now that we have a selected holo region, it is time to truncate it to create a fully treated holo region!
        self.select_cat_residues()
        cat_center_selection = RegionBuilder(name="cat_center", atom_group=self.cat_center_atoms).get_region()
        
        # Subtract catalytic center from holo region to create apo region
        apo_sub_region = holo_region.subtract(cat_center_selection)
        apo_builder = RegionBuilder(name="CSA_apo", universe=self.model.universe)
        
        for atom in apo_sub_region.atoms:
            apo_builder.init_atom(atom)
            
        apo_region = apo_builder.get_region()
        apo_region.universe = self.model.universe

        # Examining if alanine mutation is requested by user
        if self.alanine_mutation is not None and len(self.alanine_mutation) > 0:
            BetaCarbon(apo_region, alanine_mutation=self.alanine_mutation, name="apo_Ala_mutation")

        # Set name of apo region and assign QM menthod
        apo_method = self.method.__class__.__new__(self.method.__class__)
        apo_method.__dict__.update(self.method.__dict__)
        apo_method.assign_to_region(apo_region)

        del self.model.truncated
        apo_region.name = "CSA_apo"
        self.model.set_region(apo_region)
        apo_region.set_creation_attr(self.apo_selection_attr)
        self.model.truncate()

        # Method and truncation for apo region should be the same as holo region, so I think we can just write input
        self.model.write_input(**kwargs)

        del self.model.truncated

        with open("CSA_holo_apo.pkl", "wb") as file:
            pickle.dump(self.model, file)

        self.region = CSA_initial_region


    def select_atoms(self):
        import cclib
        """
        """
        # List to hold selected atoms based on CSA cutoff
        CSA_selected_atoms = []
        CSA_selected_resid = []
        delta_resid_charges = {}
        holo_resid_charges = {}
        apo_resid_charges = {}
        all_atoms = self.model.universe.select_atoms("all")

        # Check if the required CSA regions exist within the regions list
        required_regions = ["CSA_holo_truncated", "CSA_apo_truncated"]

        # self.model.regions is a LIST of objects. We extract the .name attribute from each
        existing_region_names = [region.name for region in self.model.regions]
        missing_regions = [r for r in required_regions if r not in existing_region_names]

        if not missing_regions:
            # To retrieve the actual objects, we find them by looking at the name in the list
            CSA_holo = next(r for r in self.model.regions if r.name == "CSA_holo_truncated")
            CSA_apo = next(r for r in self.model.regions if r.name == "CSA_apo_truncated")
            print(f"Successfully retrieved {', '.join(required_regions)} from the QMzymeModel.")
        else:
            raise UserWarning(
                f"Required regions {missing_regions} not found in the QMzymeModel."
                "Please run the Charge Shift Analysis (CSA) preprocessing to generate these regions.")

        # Read the SP calculation outputs for both holo and apo regions
        holo_output = cclib.io.ccread(self.holo_output_file)
        apo_output = cclib.io.ccread(self.apo_output_file)

        CSA_holo = self.model.CSA_holo_truncated
        CSA_apo = self.model.CSA_apo_truncated

        # Raising error if length or atomic number does not match between .pkl and output file
        for atoms, atomnos in ((CSA_holo.atoms, holo_output.atomnos),(CSA_apo.atoms,  apo_output.atomnos),):
            if len(atoms) != len(atomnos):
                raise UserWarning(".pkl region and output files do not match. Check your files.")

            for a, n in zip(atoms, atomnos):
                if a.atomic_number != n:
                    raise UserWarning(".pkl region and output files do not match. Check your files.")

        # Extract hirshfeld charges for each atom
        holo_atom_charges = holo_output.atomcharges[self.pop]
        apo_atom_charges = apo_output.atomcharges[self.pop]

        for atom, charge in zip(CSA_holo.atoms, holo_atom_charges):
            setattr(atom, "Partial Charge", charge)

        for atom, charge in zip(CSA_apo.atoms, apo_atom_charges):
            setattr(atom, "Partial Charge", charge)

        # Get the residue IDs for each atom in both holo and apo regions for per residue calculation
        holo_atom_resids = [atom.resid for atom in CSA_holo.atoms]
        apo_atom_resids = [atom.resid for atom in CSA_apo.atoms]
        
        # Sum per-residue charges for holo
        for i, resid in enumerate(holo_atom_resids):
            holo_resid_charges[resid] = holo_resid_charges.get(resid, 0.0) + holo_atom_charges[i]

        # Sum per-residue charges for apo
        for i, resid in enumerate(apo_atom_resids):
            apo_resid_charges[resid] = apo_resid_charges.get(resid, 0.0) + apo_atom_charges[i]

        # Find residues that exist in both holo and apo
        common_resids = set(holo_resid_charges.keys()) & set(apo_resid_charges.keys())

        # Loop over residues and calculate delta charges
        for resid in common_resids:
            delta_resid_charges[resid] = abs(holo_resid_charges[resid] - apo_resid_charges[resid])

        # Select atoms based on charge difference cutoff. Cutoff is set during class instantiation
        else:
            for resid, delta_charges in delta_resid_charges.items():
                if delta_charges >= self.charge_threshold:
                    CSA_selected_resid.append(resid)
                    resname = next(
                        (atom.resname for atom in CSA_holo.atoms if atom.resid == resid),
                        "UNK")
            for atom in all_atoms:
                if atom.resid in CSA_selected_resid:
                    CSA_selected_atoms.append(atom)

        # Get the indices of the selected atoms
        indices = [atom.index for atom in CSA_selected_atoms]

        # Build an AtomGroup from the universe
        CSA_atom_group = self.model.universe.atoms[indices]
        CSA_region = RegionBuilder(name="CSA_cutoff_region", atom_group=CSA_atom_group).get_region()
        region = CSA_region.combine(self.model.catalytic_center)

        # Making a csv file of atom, resid, and delta charge data
        holo_atom_data = []
        for i, atom in enumerate(CSA_holo.atoms):
            resid = atom.resid
            holo_atom_data.append({
                "holo_atom_id": i + 1,
                "holo_atom_charge": holo_atom_charges[i]
            })

        df_holo_atom = pd.DataFrame(holo_atom_data)
        df_holo_atom.to_csv(f"holo_charges_{self.pop}.csv", index=False)

        holo_resid_data = []
        for i, atom in enumerate(CSA_holo.atoms):
            resid = atom.resid
            holo_resid_data.append({
                "holo_residue_id": resid,
                "holo_residue_name": atom.resname,
                "holo_residue_charge": holo_resid_charges[resid]
            })

        df_holo_resid = pd.DataFrame(holo_resid_data)
        df_holo_resid.to_csv(f"holo_residue_charges_{self.pop}.csv", index=False)

        apo_atom_data = []
        for i, atom in enumerate(CSA_apo.atoms):
            resid = atom.resid
            apo_atom_data.append({
                "apo_atom_id": i + 1,
                "apo_atom_charge": apo_atom_charges[i]
            })

        df_apo_atom = pd.DataFrame(apo_atom_data)
        df_apo_atom.to_csv(f"apo_charges_{self.pop}.csv", index=False)

        apo_resid_data = []
        for i, atom in enumerate(CSA_apo.atoms):
            resid = atom.resid
            apo_resid_data.append({
                "apo_residue_id": resid,
                "apo_residue_name": atom.resname,
                "apo_residue_charge": apo_resid_charges[resid]
            })
        
        df_apo_resid = pd.DataFrame(apo_resid_data)
        df_apo_resid.to_csv(f"apo_residue_charges_{self.pop}.csv", index=False)

        delta_data = []
        for resid, delta in delta_resid_charges.items():
            resname = next(
                (atom.resname for atom in CSA_holo.atoms if atom.resid == resid),
                next((atom.resname for atom in CSA_apo.atoms if atom.resid == resid), "UNK")
            )
            delta_data.append({
                "delta_charge_id": resid,
                "delta_residue_name": resname,
                "delta_charge_difference": delta
            })

        df_delta = pd.DataFrame(delta_data)
        df_delta.to_csv(f"delta_charges_{self.pop}.csv", index=False)

        # Removing methods from the model
        self.model.CSA_holo_truncated.method = None
        self.model.CSA_holo.method = None
        self.model.CSA_apo_truncated.method = None
        self.model.CSA_apo.method = None
        
        self.region = region
        
    def reference(self):
        """
        """
        self.reference = "(1) Kulik, Heather J.; Zhang, Jianyu; Klinman, Judith P.; Martinez, Todd J.(2016) How Large Should the QM Region Be in QM/MM Calculations? The Case of Catechol O-Methyltransferase. Journal of Physical Chemistry B, 120(44). and (2) Karelina, M., & Kulik, H. J. (2017). Systematic quantum mechanical region determination in QM/MM simulation. Journal of chemical theory and computation, 13(2)."
