###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
Product of the RegionBuilder class.
"""

from typing import TYPE_CHECKING, Any, Dict, Generic, List, Optional, TypeVar
from QMzyme.QMzymeAtom import QMzymeAtom
import warnings
from QMzyme import MDAnalysisWrapper as MDAwrapper

_QMzymeAtom = TypeVar("_QMzymeAtom", bound="QMzymeAtom")

class QMzymeRegion:
    def __init__(self, name, atoms: list, atom_group = None, layer = None):
        self.name = name
        self.atoms = atoms
        self.atom_group = atom_group
        self.layer = layer
        self.method = None

    def __repr__(self):
        return f"<QMzymeRegion {self.name} contains {self.n_atoms} atom(s) and {self.n_residues} residue(s)>"
    
        
    @property
    def ids(self):
        """
        Returns a list of the atom numbers/ids from the original starting structure. An atom id
        of an atom from the starting structure should not change. See idxs as an alternative.
        """
        return [atom.id for atom in self.atoms]
    
    @property
    def idxs(self):
        """
        Returns a list of the atom indices starting from 0. If the order of atoms changes the idx 
        assigned to an atom will change. See ids as an alternative.
        """
        return [idx for idx in range(self.n_atoms)]
    
    @property
    def resids(self):
        return sorted(list(set([atom.resid for atom in self.atoms])))
    
    @property
    def n_atoms(self):
        return len(self.atoms)
    
    @property
    def n_residues(self):
        #self.residues = list(set([atom.resid for atom in self.atoms]))
        #return len(self.resids)
        return len(self.residues)
    
    @property 
    def residues(self):
        residues = []
        for resid in self.resids:
            atoms = [atom for atom in self.atoms if atom.resid == resid]
            resname = atoms[0].resname
            res = QMzymeResidue(resname, resid, atoms)
            residues.append(res)
        return residues
    
    def set_layer(self, value: str):
        """
        Accepted layer values are "QM", "xTB" and "point_charges".
        """
        value = value.lower()
        if value not in ["qm", "xtb", "point_charges"]:
            raise UserWarning('Accepted layer values are QM, xTB and point_charges.')
        self.layer = value

    def set_atom_group(self, atom_group):
        self.atom_group = atom_group
        
    def get_atom_group(self):
        return self.atom_group
    
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
    
    def add_atom(self, atom: _QMzymeAtom):
        """
        :param atom: The atom you want to add to the QMzymeRegion. 
        :type atom: _QMzymeAtom. 
        """
        self.atoms.append(atom)
        self.atoms = self.sort_atoms()

    def sort_atoms(self):
        atoms = self.atoms
        ids = [atom.id for atom in self.atoms]
        return [x for _, x in sorted(zip(ids, atoms))]
    
    def get_residue(self, resid):
        for res in self.residues:
            if res.resid == resid:
                return res

    def write(self, filename=None):
        warnings.filterwarnings('ignore')
        # Housekeeping
        if filename is None:
            filename = f"{'_'.join(self.name.split(' '))}.pdb"
        ag = self.convert_to_AtomGroup()
        ag.write(filename)
        return filename

    def convert_to_AtomGroup(self):
        return MDAwrapper.build_universe_from_QMzymeRegion(self)
    
    def set_fixed_atoms(self, ids: list):
        for atom in self.atoms:
            if atom.id in ids:
                #atom.set_fixed(value = True)
                setattr(atom, "is_fixed", True)

    def get_ids(self, attribute: str, value):
        """
        Example: get_ids(attribute='type', value='CA')
        """
        ids = []
        for atom in self.atoms:
            if getattr(atom, attribute) == value:
                ids.append(atom.id)
        return ids
    
    def get_indices(self, attribute: str, value):
        ids = self.get_ids(attribute, value)
        return self.get_idxs_from_ids(ids)

    def get_idxs_from_ids(self, ids):
        """
        Example: get_ids(attribute='type', value='CA')
        """
        idxs = []
        for idx, atom in enumerate(self.atoms):
            if atom.id in ids:
                idxs.append(idx)
        return idxs
    
    def check_missing_attr(self, attr):
        missing = []
        for atom in self.atoms:
            if not hasattr(atom, attr) or getattr(atom, attr) == None:
                missing.append(atom)
        if missing != []:
            raise UserWarning(f"The following atoms are missing {attr} information: {missing}")
        
    def set_method(self, method, _type):
        if type(method) != dict:
            method = method.__dict__
        method["type"] = _type
        self.method = method

class QMzymeResidue(QMzymeRegion):
    def __init__(self, resname, resid, atoms, chain=None):
        self.resname = resname
        self.resid = resid
        self.atoms = atoms
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



