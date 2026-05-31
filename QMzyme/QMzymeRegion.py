###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@lsu.edu
###############################################################################

from typing import TYPE_CHECKING, Any, Dict, Generic, List, Optional, TypeVar
from QMzyme.QMzymeAtom import QMzymeAtom
from QMzyme.converters import region_to_atom_group
import warnings
import numpy as np
import copy
from QMzyme import MDAnalysisWrapper as MDAwrapper
from QMzyme.data import residue_charges, backbone_atoms


_QMzymeAtom = TypeVar("_QMzymeAtom", bound="QMzymeAtom")

class QMzymeRegion:
    """
    Product of the RegionBuilder class. A QMzymeRegion is composed of QMzymeAtom objects. 

    .. seealso:: `QMzymeRegion.QMzymeResidue`
    """
    
    def __init__(self, name, atoms: list, atom_group = None, universe = None):
        self.name = name
        self.atoms = atoms
        self.atoms = self.sorted_atoms()
        self._atom_group = atom_group
        if atom_group is not None:
            universe = atom_group.universe
        self._universe = universe
        if not hasattr(self, "method"):
            self.method = None
        self._selection_attr = {}
        self._creation_attr = {}

    def __repr__(self):
        return f"<QMzymeRegion {self.name} contains {self.n_atoms} atom(s) and {self.n_residues} residue(s)>"
    
    def __sub__(self, other, name=None):
        """
        Remove any atoms with IDs matching atom IDs in `other` region.
        """
        diff_array = np.setdiff1d(np.array(self.ids), np.array(other.ids))
        atoms = [self.get_atom(id=id) for id in diff_array]
        if name == None:
            name = f"{self.name}_{other.name}_subtracted"
        return QMzymeRegion(name=name, 
                            atoms=atoms, 
                            universe=self._universe)
    
    def __add__(self, other, name=None):
        """
        Combine two QMzymeRegions. If an atom with the same ID appears in both regions, the attributes
        from the first (self) region will be kept.
        """
        if name == None:
            name = f"{self.name}_{other.name}_combined"
        region = QMzymeRegion(name=name, 
                              atoms=other.atoms, 
                              universe=self._universe)
        for atom in [self.get_atom(id=id) for id in self.ids]:
            region.add_atom(atom, override_same_id=True)
        return region

    def __eq__(self, other):
        if other == None:
            return
        eq_ids = np.array_equal(self.ids, other.ids)
        names = [atom.name for atom in self.atoms]
        names_other = [atom.name for atom in other.atoms]
        eq_names = np.array_equal(names, names_other)

        return (eq_ids, eq_names) == (True, True)

    @property
    def ids(self):
        """
        A list of atom numbers/ids from the original starting structure. An atom id
        of an atom from the starting structure should not change.

        """
        return [atom.id for atom in self.atoms]
    
    @property
    def ix_array(self):
        """
        A list of atom indices starting from 0. If the order of atoms changes the ix 
        assigned to an atom will change. See ids as an alternative.

        .. seealso:: :func:`~QMzyme.QMzymeRegion.QMzymeRegion.ids`
        """
        return [ix for ix in range(self.n_atoms)]
        #return np.array([atom.ix for atom in self.atoms])
    
    @property
    def resids(self):
        """
        :returns: Sorted residue IDs for residues within this region.
        :rtype: list
        """
        return sorted(list(set([atom.resid for atom in self.atoms])))
    
    @property
    def n_atoms(self):
        """
        :returns: Number of atoms within this region.
        :rtype: int
        """
        return len(self.atoms)
    
    @property
    def n_residues(self):
        """
        :returns: Number of residues within this region.
        :rtype: int
        """
        return len(self.residues)
    
    @property 
    def residues(self):
        """
        :returns: Residues in this region.
        :rtype: List[:class:`~QMzyme.QMzymeRegion.QMzymeResidue`]
        """
        residues = []
        for resid in self.resids:
            atoms = [atom for atom in self.atoms if atom.resid == resid]
            resname = atoms[0].resname
            res = QMzymeResidue(resname, resid, atoms, region=self)
            residues.append(res)
        return residues
    
    @property
    def positions(self):
        """
        :returns: Positions of all atoms in region.
        :rtype: NumPy array
        """
        coordinates = np.empty((self.n_atoms, 3), dtype=np.float32)
        for i, atom in enumerate(self.atoms):
            coordinates[i] = atom.position
        return coordinates

    @property
    def atom_group(self):
        """
        Returns the :class:`~MDAnalysis.core.groups.AtomGroup` instance. If the region was not 
        bult from an AtomGroup initially, it will be converted to one. Note, in this case, the universe
        of that AtomGroup will 
        """
        if self._atom_group is None or self._atom_group.n_atoms != self.n_atoms:
            return region_to_atom_group(self)
        return self._atom_group
    
    @property
    def segids(self):
        if hasattr(self.atoms[0], 'segid'):
            return [atom.segid for atom in self.atoms]
        else:
            return
    
    @property
    def creation_attr(self):
        """
        Filters _selection_attr into _creation_attr.
        Ultimately, it creates a property to the region of how the region was created.
        :rtype: List[:class:`~QMzyme.QMzymeRegion.creation_attr`]
        """
        # Clear previous state to avoid doubling up
        self._creation_attr = {}

        # Creates selection_scheme first
        if 'selection_scheme' in self._selection_attr:
            v = self._selection_attr['selection_scheme']
            self._creation_attr['selection_scheme'] = v.__name__ if hasattr(v, '__name__') else str(v)

        # Loop through everything else
        for k, v in self._selection_attr.items():

            # Skip method since it is saved elsewhere
            if k == 'method' or k == 'selection_scheme':
                continue
            
            # Append the list of other values within the attribute of set_region
            else:
                self._creation_attr[k] = v
        
        return self._creation_attr
    
    def reset_creation_attr(self):
        """
        Resets the selection and creation parameters for the region.
        """
        self._selection_attr = {}
        self._creation_attr = {}

    def set_creation_attr(self, attr=None, **kwargs):
        """
        Manually sets or updates the selection parameters for the region. This is
        particularly useful for assigning parameters to internally generated regions
        that need to inherit the creation history of their parent selection scheme.
        """
        
        if not hasattr(self, '_selection_attr'):
            self._selection_attr = {}
            
        if attr is not None:
            self._selection_attr.update(attr)
            
        if kwargs:
            self._selection_attr.update(kwargs)

    def get_atom(self, id):
        """
        Selects the atom with set atom ID.

        :param id: The id of an atom you want to select.
        :type id: float
        """
        for i in self.atoms:
            if i.id == id:
                return i
            
    def has_atom(self, id):
        """
        Examines if specific atom exists within the region.

        :param id: The ID of an atom you want to check for existance in the region.
        :type id: int
        """
        if id in self.ids:
            return True
        return False
    
    def has_residue(self, resid):
        """
        Examines if specific residue exists within the region.

        :param resid: The residue ID you want to check for existance in the region.
        :type resid: int
        """
        if resid in self.resids:
            return True
        return False
    
    def add_atom(self, atom: _QMzymeAtom, override_same_id=False):
        """
        Adds additional QMzymeAtom to the QMzymeRegion.

        :param atom: The atom you want to add to the QMzymeRegion. 
        :type atom: :class:`~QMzyme.QMzymeAtom.QMzymeAtom`, required
        :param overrisde_same_id: An argument to decide if the atoms with same IDs are replaced.
        :type overrisde_same_id: bool, optional

        .. warning:: Ths will modify the QMzymeRegion directly.
        """
        self.atoms.append(atom)
        self.atoms = self.sorted_atoms(override_same_id=override_same_id)

    def remove_atom(self, atom: _QMzymeAtom):
        """
        Removes QMzymeAtom from the QMzymeRegion.

        :param atom: The atom you want to remove from the QMzymeRegion. 
        :type atom: :class:`~QMzyme.QMzymeAtom.QMzymeAtom`, required
        :param overrisde_same_id: An argument to decide if the atoms with same IDs are replaced.
        :type overrisde_same_id: bool, optional

        .. warning:: Ths will modify the QMzymeRegion directly.
        """
        self.atoms.remove(atom)

    def sorted_atoms(self, override_same_id=False):
        """
        Returns a list of atoms sorted by their IDs, with optional duplicate handling.
        This method checks the most recently added atom's ID against prexisting atoms' IDs.
        If a duplicate ID is found and override_same_id=False, warning is raised.
        If override_same_id=True, it will replace older atom and replace it with newer atom
        with same ID.

        :param overrisde_same_id: An argument to decide if the atoms with same IDs are replaced.
        :type overrisde_same_id: bool, optional

        :return: A list of QMzymeAtom objects sorted numerically by their atom IDs.
        :rtype: list
        """
        if self.atoms != [] and self.atoms[-1].id in self.ids[:-1]:
            if override_same_id == False: 
                #self.remove_atom(self.atoms[-1])
                raise UserWarning(f"Atom {self.atoms[-1]} cannot be added to region because another atom with the same id already exists: {self.get_atom(self.atoms[-1].id)}.")
            else:
                atom = self.get_atom(self.atoms[-1].id)
                self.remove_atom(atom)
        atoms = self.atoms
        ids = [atom.id for atom in self.atoms]
        return [x for _, x in sorted(zip(ids, atoms))]
    
    def get_residue(self, resid):
        """
        Retrieve a residue from the QMzymeRegion.

        :param resid: The resid of the QMzymeResidue within QMzymeRegion.
        :type resid: int
        :return: The corresponding QMzymeResidue object if found, otherwise None.
        :rtype: QMzymeResidue
        """
        for res in self.residues:
            if res.resid == resid:
                return res
            
    def rename(self, name):
        """
        Rename the region.

        :param name: The new name for the QMzymeRegion.
        :type name: str
        """
        self.name = name

    def write(self, filename=None, format='pdb'):
        """
        Converts QMzymeRegion to an `MDAnalysis AtomGroup <https://userguide.mdanalysis.org/stable/atomgroup.html>`_,
        then uses the AtomGroup method to write file. 
        
        :param filename: Name of file. If not specified, the region name attribute will be used. If the region name 
            is an empty string the name will become 'noname'. Note, if you specify a file format in filename (i.e., filename='name.xyz'),
            but the value of the 'format' argument does not match, the value assigned to format will be used.
        :type filename: str, optional

        :param format: Format of created file. 
        :type format: str, default='pdb'
        
        .. note:: Any file ending supported by the MDAnalysis AtomGroup write method is
            supported as long as all the necessary information to write that file type is present in the region. 
        """
        from QMzyme.utils import check_filename
        warnings.filterwarnings('ignore')
        # Housekeeping
        if filename is None:
            if self.name == '':
                self.name='noname'
            filename = f"{'_'.join(self.name.split(' '))}.{format}"
        filename = check_filename(filename, format)
        #ag = self.convert_to_AtomGroup()
        #ag.write(filename)
        self.atom_group.write(filename)
        return filename

    # def convert_to_AtomGroup(self):
    #     return MDAwrapper.build_universe_from_QMzymeRegion(self)
    
    def set_fixed_atoms(self, ids: list= None, atoms: list=None):
        """
        Example Usage: 
        To fix all alpha carbons
        >>> ids = get_ids(attribute='type', value='CA') 
        >>> set_fixed_atoms(ids=ids)

        :param ids: Atom ids in QMzymeRegion to fix.
        :type ids: List[int], optional, default=None

        :param atoms: Atoms in QMzymeRegion to fix.
        :type atoms: List[:class:`~QMzyme.QMzymeAtom.QMzymeAtom`], default=None

        .. note:: Must specify either ids or atoms (if both are 
            specified, atoms will be used). During calculation file writing if
            an atom has attribute ``is_fixed=True``, that atom will be frozen.
        """
        if atoms is not None:
            for atom in atoms:
                setattr(atom, "is_fixed", True)
        elif ids is not None:
            for atom in self.atoms:
                if atom.id in ids:
                    #atom.set_fixed(value = True)
                    setattr(atom, "is_fixed", True)

    def get_ids(self, attribute: str, value):
        """
        :Example: 
        To fix all alpha carbons
        >>> ids = get_ids(attribute='type', value='CA') 
        >>> set_fixed_atoms(ids=ids)

        :param attribute: QMzymeAtom object attribute name.
        :type attribute: str, required

        :param value: Value of interest for the corresponding attribute. 
        :type value: Depends on the attribute, required

        :returns: atom ids.
        :rtype: List[int]

        .. seealso:: :func:`~QMzyme.QMzymeRegion.QMzymeRegion.get_atoms`
        """
        ids = []
        for atom in self.atoms:
            if getattr(atom, attribute) == value:
                ids.append(atom.id)
        return ids
    
    def get_atoms(self, attribute: str, value):
        """
        :Example: 
            To see what atoms are in residue with resid 14.
            > atoms = get_atoms(attribute='resid', value=14) 

        :param attribute: QMzymeAtom object attribute name.
        :type attribute: str, required

        :param value: Value of interest for the corresponding attribute. 
        :type value: Depends on the attribute, required

        :returns: atom ids.
        :rtype: List[int]

        List[:class:`~QMzyme.QMzymeAtom.QMzymeAtom`]
         
        .. seealso:: :func:`~QMzyme.QMzymeRegion.QMzymeRegion.get_ids`
        """
        atoms = []
        for atom in self.atoms:
            if getattr(atom, attribute) == value:
                atoms.append(atom)
        return atoms

    def get_indices(self, attribute: str, value):
        """
        Get the indices of QMzymeAtoms that match a specific attribute value.

        :param attribute: The name of the atom attribute to filter by (e.g., 'resname', 'Val').
        :type attribute: str
        :param value: The value the attribute must match.
        :type value: str

        :return: A list of integer indices representing QMzymeAtom IDs.
        :rtype: list[int]
        """
        ids = self.get_ids(attribute, value)
        return self.get_ix_array_from_ids(ids)

    def get_ix_array_from_ids(self, ids):
        """
        :Example: 
        >>> ixs = get_ids(attribute='type', value='CA')
        """
        ix_array = []
        for ix, atom in enumerate(self.atoms):
            if atom.id in ids:
                ix_array.append(ix)
        return ix_array
    
    def check_missing_attr(self, attr):
        """
        Verify that all atoms in the region possess a specific attribute.

        Checks all QMzymeAToms within the QMzymeRegion for the presence of
        specific attribute. If any QMzymeAtoms are missing the attribute, a warning is raised.

        :param attr: The name of the attribute.
        :type attr: str
        :raises UserWarning: If one or more QMzymeAtoms lack the attribute or if the attribute is None.
        """
        missing = []
        for atom in self.atoms:
            if not hasattr(atom, attr) or getattr(atom, attr) == None:
                missing.append(atom)
        if missing != []:
            raise UserWarning(f"The following atoms are missing {attr} information: {missing}")
        
    def set_method(self, method):
        """
        This method designates a specific QM method to the QMzymeRegion. This is used by the
        :module:`~QMzyme.CalculateModel` module when the region is passed to a Calculation Method class.

        :param method: The method that is going to be applied to QMzymeRegion.
        :type method: dict
        """
        if type(method) != dict:
            method = method.__dict__
        self.method = method

    def set_charge(self, charge):
        """
        Sets specific charge to QMzymeRegion
        :param charge: The charge of the QMzymeRegion.
        :type charge: int
        """
        self.charge = charge
        try:
            self.method["charge"] = charge
        except:
            pass

    def guess_charge(self, verbose=True):
        """
        Guesses charge based on the residue_charges information in QMzyme\configuration\__init__.py.
        QMzyme contains charge information of standard AMBER amino acid residues.
        If non-AMBER residues are present, it will raise an error. To update the charge of the
        unknown residue, the user can use QMzyme.data.residue_charges.update({'unknown residue name': int})

        :param verbose: Returns print statements including warning and estimated charge.
        :type verbose: bool
        """
        if hasattr(self.atoms[0], "charge"):
            self.read_charges(verbose)
            return
        txt = ''
        txt += f"\nEstimating total charge for QMzymeRegion {self.name} based on residue naming conventions in QMzyme.data.residue_charges..."
        unk_res = []
        chrg = 0
        for res in self.residues:
            if res.resname not in residue_charges:
                unk_res.append(res)
                txt+=f"\n\t{res} --> Charge: UNK, defaulting to 0"
            else: 
                q = residue_charges[res.resname.upper()]
                chrg += q
                txt+=f"\n\t{res} --> Charge: {q}"
        self.set_charge(chrg)
        if unk_res == []:
            txt+=f"\nQMzymeRegion {self.name} has an estimated charge of {chrg}."
        else:    
            txt+=f"\nWARNING: Charge estimation may be inaccurate due to presence of residue(s) with unknown charge: {unk_res}." 
            txt+=f"QMzymeRegion {self.name} has an estimated total charge of {chrg}."
        if verbose == True:
            print(txt)
        else:
            print(txt.split('\n')[-1])


    def read_charges(self, verbose=True):
        """
        Calculates total charge of QMzymeResidue by using tolopogy attribute 'charge'.

        :param verbose: Returns print statements including warning and estimated charge.
        :type verbose: bool
        """
        txt=''
        txt+=f"\nCalculating total charge for QMzymeRegion {self.name} based on charges read from topology attribute 'charge'..."
        chrg = 0
        #for atom in self.atoms:
        #    chrg += atom.charge
        #chrg = round(chrg)
        self.set_charge(chrg)
        for res in self.residues:
            q = 0
            for atom in res.atoms:
                q+=atom.charge
            txt+=f"\n\t{res} --> Charge: {round(q)}"
            chrg += round(q)
            #print(f"\n\t{res} --> Charge: {q}, rounded {round(q)}")
        chrg = round(chrg)
        self.set_charge(chrg)
        txt+=f"\nQMzymeRegion {self.name} has a total charge of {chrg}."
        if verbose == True:
            print(txt)


    def combine(self, other, name=None):
        """
        Combine region with another region. 

        :param other: Region containing atoms to combine with self.
        :type other: :class:`~QMzyme.QMzymeRegion.QMzymeRegion`

        :param name: Name of returned region. If no name is specified the name will be 
            'self.name_other.name_combined.'.
        :type name: str, default=None
        
        :returns: Combined region.
        :rtype:  :class:`~QMzyme.QMzymeRegion.QMzymeRegion`

        .. note:: Self region is not altered, and returned region will not contain 
            properties or attributes of self region not described from the QMzymeAtom level 
            (i.e., atom_group) except for _universe. Atoms found in both self and other will 
            only be copied from self. This is important if you have different attribute values 
            on an atom that appears in both regions (i.e., is_fixed=True and is_fixed=False).
        """
        combined_atoms = copy.copy(self.atoms)
        for atom in other.atoms:
            # if not atom.is_within(self):
            # if not atom.id in self.ids:
            #     combined_atoms.append(atom)
            # for atom2 in self.atoms:
            #     if atom != atom2:
            #         combined_atoms.append(atom)
            if atom not in self.atoms:
                if atom.id in self.ids:
                    atom = copy.copy(atom)
                    atom.id = np.max([a.id for a in combined_atoms])+1
                combined_atoms.append(atom)
        if name == None:
            name = f"{self.name}_{other.name}_combined"
        combined_region = QMzymeRegion(name=name, atoms=combined_atoms, universe=self._universe)
        setattr(combined_region, "universe", self.atom_group.universe)
        return combined_region
    
    def subtract(self, other, name = ''):
        """
        Creates a new QMzymeRegion that does not contain any atoms found in other. 

        :param other: Region containing atoms to remove from self.
        :type other: :class:`~QMzyme.QMzymeRegion.QMzymeRegion`

        :param name: Name of returned region.
        :type name: str, default=''
        
        :returns: New region with atoms found in other region removed.
        :rtype:  :class:`~QMzyme.QMzymeRegion.QMzymeRegion`

        .. note:: Self region is not altered, and returned region will not contain 
            properties or attributes of self region not described from the QMzymeAtom level 
            (i.e., atom_group) except for _universe.
        """
        atoms = []
        for atom in self.atoms:
            if not atom.is_within(other):
                atoms.append(atom)
        region = QMzymeRegion(name=name, atoms=atoms, universe=self._universe)
        return region
    
    def get_overlapping_atoms(self, other):
        """
        :param other: Other QMzymeRegion to measure overlap with.
        :type other: :class:`~QMzyme.QMzymeRegion.QMzymeRegion`
        
        :returns: Atoms in self that are also found in other.
        :rtype: :class:`~QMzyme.QMzymeAtom.QMzymeAtom`
        """
        atoms = []
        for atom in self.atoms:
            if atom in other.atoms:
                atoms.append(atom)
        return atoms
    
    def set_atom_segid(self, segid):
        """
        Sets segment id attribute to QMzymeAtom objects within QMzymeResidue.

        :param segid: segment id of the QMzymeAtom
        :type segid: int
        """
        for atom in self.atoms:
            atom.segid = segid
    
    # def guess_bonds():
    #     """
    #     Method under development.
    #     """
    #     pass

    def summarize(self, filename=None):
        """
        Creates a summary of QMzymeRegion, which includes resid, resname, charge, removed atoms, and fixed atoms
        of QMzymeResidues within QMzymeRegion. If a file name is specified, the summary will be created as a .txt
        file. Else, it will return a dictionary with summary information.

        :param filename: Name of the output file for the summarized information.
        :type filename: str
        """
        summary = {
            "Resid": [],
            "Resname": [],
            "Charge": [],
            "Removed atoms": [],
            "Fixed atoms": [],
        }
        for res in self.residues:
            summary["Resid"].append(res.resid)
            summary["Resname"].append(res.resname)
            if not hasattr(res, "charge"):
                res.guess_charge(verbose=False)
            summary["Charge"].append(res.charge)
            summary["Removed atoms"].append(res.removed_atoms)
            summary["Fixed atoms"].append([a.name for a in res.get_atoms('is_fixed', True)])
        summary["Segids"] = [res.atoms[0].segid for res in self.residues]
        if filename == None:
            return summary
        if not filename.endswith('.txt'):
            filename += '.txt'
        with open(filename, 'w') as f:
            print(summary, file=f)

    def align_to(self, other, self_selection='all', other_selection='all', update_region=True):
        """
        This method aligns the selection of QMzymeAtom objects within self QMzymeRegion to the
        specified selection of QMzymeAtom objects from a QMzymeRegion. The computed RMSD of pre-
        and post-alignment is printed. It is important to align the smaller subset of QMzymeAtom
        objects (self) to the QMzymeAtom objects of larger QMzymeRegion (other) to align the smaller
        subset to the larger QMzymeRegion. This can be ensured by specifying the smaller model before
        the “align_to” function. This will assign the smaller model to “self”.

        :param other: QMzymeRegion of the other QMzymeAtoms for alignment.
        :type other: QMzymeRegion
        :param self_selection: Selection of QMzymeAtoms based on MDAnalysis selection nomenclatire.
        :type self_selection: str
        :param other_selection: Selection of QMzymeAtoms based on MDAnalysis selection nomenclatire.
        :type other_selection: str
        :parm update_region: Updates the aligned QMzymeAtom objects to the QMzymeRegion.
        :type update_region: bool
        """
        from QMzyme.utils import compute_translation_and_rotation, kabsch_transform, rmsd
        mobile = self.atom_group.select_atoms(self_selection)
        target = other.atom_group.select_atoms(other_selection)
        rmsd_before_alignment = rmsd(mobile.positions, target.positions)
        if len(mobile.atoms) != len(target.atoms):
            raise UserWarning("The same number of atoms must be selected for alignment. Please adjust selections.")
        t, r = compute_translation_and_rotation(mobile.positions, target.positions)
        aligned_positions = kabsch_transform(self.positions, t, r)
        mobile_aligned_positions = kabsch_transform(mobile.positions, t, r)
        rmsd_after_alignment = rmsd(mobile_aligned_positions, target.positions)
        print(f"RMSD before alignment: {rmsd_before_alignment} \u00C5")
        print(f"RMSD after alignment: {rmsd_after_alignment} \u00C5")
        if update_region is True:
            self._atom_group.positions = aligned_positions
            for i, atom in enumerate(self.atoms):
                atom.position = aligned_positions[i]
        else:
            return aligned_positions
        
    def find_nearby_residues(self, other, dist):
        """
        Finds residues in 'other' region that have at least one atom within 'dist' 
        Angstroms of any atom in 'self' region.
        
        :param other: The target region to search within for nearby residues.
        :type other: :class:`~QMzyme.QMzymeRegion.QMzymeRegion`

        :param dist: The distance cutoff in Angstroms.
        :type dist: float

        :returns: List of unique residues from the other region within the cutoff distance.
        :rtype: list of :class:`~QMzyme.QMzymeResidue.QMzymeResidue`
        
        Returns:
        list of unique QMzymeResidue objects from the 'other' region.
        """

        nearby_residues = []

        for self_res in self.residues:
            for other_res in other.residues:
                if self_res.resid == other_res.resid:
                    continue
                
                contact_found = False

                for self_atom in self_res.atoms:
                    for other_atom in other_res.atoms:
                        d = np.linalg.norm(self_atom.position - other_atom.position)

                        if d <= dist:
                            contact_found = True
                    
                    if contact_found:
                        break
                
                if contact_found:
                    print (f"{other_res} is within {dist} \N{ANGSTROM SIGN} of {self_res}")

                    if other_res not in nearby_residues:
                        nearby_residues.append(other_res)

        return nearby_residues

    def store_calculation_results(self, calculation_file):
        """
        Once you run a calculation on a QMzymeRegion you can then store the calculation 
        results. This method uses `cclib <https://cclib.github.io/index.html>`_ to parse the output, 
        so any calculation file that cclib can read can be stored to your QMzymeRegion, and the 
        stored dictionary will contain all the `data that cclib parses <https://cclib.github.io/data.html>`_.
        """
        import cclib
        data = cclib.io.ccread(calculation_file)
        if not hasattr(self, 'calculations'):
            self.calculations = {}
        if data.metadata['methods'] == []:
            data.metadata['methods'] = method
        self.calculations[calculation_file] = data


class QMzymeResidue(QMzymeRegion):
    """
    QMzymeResidue is a subclass of QMzymeRegion representing a molecular unit of the system as defined in the
    starting topology file. 

    Required Parameters
    --------------------
    
    :param resname: Three letter residue name: ex., 'VAL'
    :type resname: str

    :param resid: Integer residue number.
    :type resid: int

    :param atoms: Atom name: ex., 'C1'
    :type name: :class:`~QMzyme.QMzymeAtom.QMzymeAtom`

    :param region:
    :type name: :class:`~QMzyme.QMzymeRegion.QMzymeResidue`
    """
    
    def __init__(self, resname, resid, atoms, region, chain=None):
        self.resname = resname
        self.resid = resid
        self.atoms = atoms
        self.region = region
        if chain is None:
            chain = self.atoms[0].get_chain()
        self.chain = chain

    def __repr__(self):
        rep =  f"<QMzymeResidue resname: {self.resname}, resid: {self.resid}, chain: "
        if self.chain is None:
            rep += "Not Specified>"
        else:
            rep += f"{self.chain}>"
        return rep

    def get_atom(self, atom_name):
        """
        Selects QmzymeAtom with specific atom name within the region.

        :param atom_name: The name of the atom (e.g. "CA")
        :type atom: str
        """
        for atom in self.atoms:
            if atom.name == atom_name:
                return atom

    def set_chain(self, value: str):
        """
        Sets chain attribute to QMzymeAtom objects within QMzymeRgion.

        :param segid: Chain name of the QMzymeAtom
        :type segid: str
        """
        self.chain = value

    def guess_charge(self, verbose=True):
        """
        Guesses charge based on the residue_charges information in QMzyme\configuration\__init__.py.
        QMzyme contains charge information of standard AMBER amino acid residues.
        If non-AMBER residues are present, it will raise an error. To update the charge of the
        unknown residue, the user can use QMzyme.data.residue_charges.update({'unknown residue name': int})

        :param verbose: Returns print statements including warning and estimated charge.
        :type verbose: bool
        """
        if hasattr(self.atoms[0], "charge"):
            self.read_charges(verbose)
            return
        txt = ''
        txt += f"\nEstimating total charge for QMzymeResidue {self.resname} based on residue naming conventions in QMzyme.data.residue_charges..."
        if self.resname not in residue_charges:
            if self.resname in ["WAT", "SOL"]:
                chrg = 0
            else:
                chrg = 'UNK'
        else: 
            chrg = residue_charges[self.resname.upper()]
        self.set_charge(chrg)
        if chrg != 'UNK':
            #self.set_charge(chrg)
            txt+=f"\nQMzymeResidue {self.resname} has an estimated charge of {chrg}."
        else:    
            txt+=f"\nQMzymeResidue {self.resname} has an unknown charge value based on conventional residue naming."
        if verbose == True:
            print(txt)
        #else:
        #    print(txt.split('\n')[-1])

    def read_charges(self, verbose=True):
        """
        Calculates total charge of QMzymeResidue by using tolopogy attribute 'charge'.

        :param verbose: Returns print statements including warning and estimated charge.
        :type verbose: bool
        """
        txt = ''
        txt+=f"\nCalculating total charge for QMzymeResidue {self.resname} based on charges read from topology attribute 'charge'..."
        chrg = 0
        for atom in self.atoms:
            chrg += atom.charge
        chrg = round(chrg)
        self.set_charge(chrg)
        txt+=f"\nQMzymeResidue {self.resname} has a total charge of {chrg}."
        if verbose == True:
            print(txt)

    def get_backbone_atoms(self, backbone_atoms=backbone_atoms):
        """
        Selects group of QMzymeAtom objects within QMzymeResidue that contains atom name within backbone_atoms dictionary.

        :param backbone_atoms: Sets of QMzymeAtom names that corresponds to the backbone atoms.
        :type backbone_atoms: dict, default {'C': 'C', 'CA': 'CA', 'H': 'H', 'HA': 'HA', 'N': 'N', 'O': 'O'}.
        """
        bb_atoms = []
        for atom_name, atom in backbone_atoms.items():
            if self.get_atom(atom) == None:
                if self.get_atom('HN') != None:
                    bb_atoms.append(self.get_atom('HN'))
                if self.get_atom('HC') != None:
                    bb_atoms.append(self.get_atom('HC'))
            else:
                bb_atoms.append(self.get_atom(atom))
        return bb_atoms

    def remove_atom(self, atom):
        """
        :param atom: The atom you want to remove from the QMzymeResidue. 
        :type atom: :class:`~QMzyme.QMzymeAtom.QMzymeAtom`, required

        .. warning:: Ths will modify the QMzymeResidue directly.
        """

        if atom in self.atoms:
            self.atoms.remove(atom)
            self.region.remove_atom(atom)

    @property
    def removed_atoms(self):
        """
        List of atoms removed from the QMzymeResidue. Often, this is done with
        remove_atom() method or truncate() method.
        """
        removed_atoms=[]
        u = self.region._universe
        sel = f'resid {self.resid} and resname {self.resname}'
        if self.chain is not None and self.chain != 'X' and self.chain != '':
            sel += f' and chainID {self.chain}'
        ag = u.select_atoms(sel)
        for atom in ag:
            if atom.name not in [a.name for a in self.atoms]:
                removed_atoms.append(atom.name)
        return removed_atoms
