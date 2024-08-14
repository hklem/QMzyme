###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
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

    def __repr__(self):
        return f"<QMzymeRegion {self.name} contains {self.n_atoms} atom(s) and {self.n_residues} residue(s)>"
    
    def __sub__(self, other, name=None):
        """
        Remove any atoms with IDs matching atom IDs in `other` reigon.
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
    
    def get_atom(self, id):
        for i in self.atoms:
            if i.id == id:
                return i
            
    def has_atom(self, id):
        if id in self.ids:
            return True
        return False
    
    def has_residue(self, resid):
        if resid in self.resids:
            return True
        return False
    
    def add_atom(self, atom: _QMzymeAtom, override_same_id=False):
        """
        :param atom: The atom you want to add to the QMzymeRegion. 
        :type atom: :class:`~QMzyme.QMzymeAtom.QMzymeAtom`, required

        .. warning:: Ths will modify the QMzymeRegion directly.
        """
        self.atoms.append(atom)
        self.atoms = self.sorted_atoms(override_same_id=override_same_id)

    def remove_atom(self, atom: _QMzymeAtom):
        """
        :param atom: The atom you want to remove from the QMzymeRegion. 
        :type atom: :class:`~QMzyme.QMzymeAtom.QMzymeAtom`, required

        .. warning:: Ths will modify the QMzymeRegion directly.
        """
        self.atoms.remove(atom)

    def sorted_atoms(self, override_same_id=False):
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
        for res in self.residues:
            if res.resid == resid:
                return res
            
    def rename(self, name):
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
        missing = []
        for atom in self.atoms:
            if not hasattr(atom, attr) or getattr(atom, attr) == None:
                missing.append(atom)
        if missing != []:
            raise UserWarning(f"The following atoms are missing {attr} information: {missing}")
        
    def set_method(self, method):
        """
        Used by the :module:`~QMzyme.CalculateModel` module when the region is passed to a Calculation Method class.
        """
        if type(method) != dict:
            method = method.__dict__
        self.method = method

    def set_charge(self, charge):
        self.charge = charge
        try:
            self.method["charge"] = charge
        except:
            pass

    def guess_charge(self, verbose=True):
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
            properties or attributs of self region not described from the QMzymeAtom level 
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
            properties or attributs of self region not described from the QMzymeAtom level 
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
        for atom in self.atoms:
            atom.segid = segid
    
    # def guess_bonds():
    #     """
    #     Method under development.
    #     """
    #     pass

    def summarize(self, filename=None):
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
        print(f"RMSD before alignment: {rmsd_before_alignment} r\AA")
        print(f"RMSD after alignment: {rmsd_after_alignment} r\AA")
        if update_region is True:
            self._atom_group.positions = aligned_positions
            for i, atom in enumerate(self.atoms):
                atom.position = aligned_positions[i]
        else:
            return aligned_positions




class QMzymeResidue(QMzymeRegion):
    """
    Subclass of QMzymeRegion.
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
        for atom in self.atoms:
            if atom.name == atom_name:
                return atom

    def set_chain(self, value: str):
        self.chain = value

    def guess_charge(self, verbose=True):
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
        if atom in self.atoms:
            self.atoms.remove(atom)
            self.region.remove_atom(atom)

    @property
    def removed_atoms(self):
        removed_atoms=[]
        u = self.region._universe
        sel = f'resid {self.resid} and resname {self.resname}'
        if self.chain is not None and self.chain != 'X':
            sel += f' and chainID {self.chain}'
        ag = u.select_atoms(sel)
        for atom in ag:
            if atom.name not in [a.name for a in self.atoms]:
                removed_atoms.append(atom.name)
        return removed_atoms
