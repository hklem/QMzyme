###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
Product of the RegionBuilder class.
"""

from typing import TYPE_CHECKING, Any, Dict, Generic, List, Optional, TypeVar
from QMzyme.QMzymeAtom import QMzymeAtom
from QMzyme.converters import region_to_atom_group
import warnings
import numpy as np
import copy
from QMzyme import MDAnalysisWrapper as MDAwrapper
from QMzyme.data import protein_residues


_QMzymeAtom = TypeVar("_QMzymeAtom", bound="QMzymeAtom")

class QMzymeRegion:
    def __init__(self, name, atoms: list, atom_group = None):
        self.name = name
        self.atoms = atoms
        self._atom_group = atom_group
        self.method = None

    def __repr__(self):
        return f"<QMzymeRegion {self.name} contains {self.n_atoms} atom(s) and {self.n_residues} residue(s)>"
    
        
    @property
    def ids(self):
        """
        Returns a list of the atom numbers/ids from the original starting structure. An atom id
        of an atom from the starting structure should not change. See ix_array as an alternative.
        """
        return [atom.id for atom in self.atoms]
    
    @property
    def ix_array(self):
        """
        Returns a list of the atom indices starting from 0. If the order of atoms changes the ix 
        assigned to an atom will change. See ids as an alternative.
        """
        return [ix for ix in range(self.n_atoms)]
    
    @property
    def resids(self):
        """
        Returns a list of sorted residue IDs for residues within this region.
        """
        return sorted(list(set([atom.resid for atom in self.atoms])))
    
    @property
    def n_atoms(self):
        """
        Returns the number of atoms within this region.
        """
        return len(self.atoms)
    
    @property
    def n_residues(self):
        """
        Returns the number of residues within this region.
        """
        return len(self.residues)
    
    @property 
    def residues(self):
        """
        Returns a list of :class:`~QMzyme.QMzymeRegion.QMzymeResidue` instances within this region.
        """
        residues = []
        for resid in self.resids:
            atoms = [atom for atom in self.atoms if atom.resid == resid]
            resname = atoms[0].resname
            res = QMzymeResidue(resname, resid, atoms)
            residues.append(res)
        return residues
    
    @property
    def positions(self):
        """
        Returns a numpy array of the positions of all atoms in this region instance.
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
        if self._atom_group is None:
            return self.convert_to_atom_group()
        return self._atom_group
    
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
        :type atom: :class:`~QMzyme.QMzymeAtom.QMzymeAtom`. 
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
            
    def rename(self, name):
        self.name = name

    def write(self, filename=None, format='pdb'):
        from QMzyme.utils import check_filename
        warnings.filterwarnings('ignore')
        # Housekeeping
        if filename is None:
            filename = f"{'_'.join(self.name.split(' '))}.{format}"
        filename = check_filename(filename, format)
        #ag = self.convert_to_AtomGroup()
        #ag.write(filename)
        self.atom_group.write(filename)
        return filename

    # def convert_to_AtomGroup(self):
    #     return MDAwrapper.build_universe_from_QMzymeRegion(self)
    
    def set_fixed_atoms(self, ids: list= None, atoms=None):
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
        Example: get_ids(attribute='type', value='CA')

        Returns List[int]
        """
        ids = []
        for atom in self.atoms:
            if getattr(atom, attribute) == value:
                ids.append(atom.id)
        return ids
    
    def get_atoms(self, attribute: str, value):
        """
        Example: get_atoms(attribute='type', value='CA')

        Returns List[QMzymeAtom]
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
        Example: get_ids(attribute='type', value='CA')
        """
        ix_array = []
        for ix, atom in enumerate(self.atoms):
            if atom.id in ids:
                ix_array.append(ix)
        return ix_array
    
    def convert_to_atom_group(self):
        return region_to_atom_group(self)
    
    def check_missing_attr(self, attr):
        missing = []
        for atom in self.atoms:
            if not hasattr(atom, attr) or getattr(atom, attr) == None:
                missing.append(atom)
        if missing != []:
            raise UserWarning(f"The following atoms are missing {attr} information: {missing}")
        
    def set_method(self, method):
        if type(method) != dict:
            method = method.__dict__
        self.method = method

    def set_charge(self, charge):
        self.charge = charge
        try:
            self.method["charge"] = charge
        except:
            pass

    def guess_charge(self):
        if hasattr(self.atoms[0], "charge"):
            self.read_charges()
            return
        print(f"\nEstimating total charge for QMzymeRegion {self.name} based on protein residue naming conventions...")
        unk_res = []
        chrg = 0
        for res in self.residues:
            if res.resname not in protein_residues:
                if res.resname in ["WAT", "SOL"]:
                    continue
                unk_res.append(res)
                print(res, f"Charge: UNK")
            else: 
                q = protein_residues[res.resname.upper()]
                chrg += q
                print(res, f"Charge: {q}")
        self.set_charge(chrg)
        if unk_res == []:
            print(f"\nQMzymeRegion {self.name} has an estimated charge of {chrg}.")
        else:    
            print(f"\n!!!Charge estimation may be inaccurate due to presence of residue(s) with unknown charge: {unk_res}. Ignoring unknown residues in charge estimation!!!")
            print(f"QMzymeRegion {self.name} has an estimated total charge of {chrg}.")


    def read_charges(self):
        print(f"\nCalculating total charge for QMzymeRegion {self.name} based on charges read from topology attribute 'charge'...")
        chrg = 0
        for atom in self.atoms:
            chrg += atom.charge
        chrg = round(chrg)
        self.set_charge(chrg)
        print(f"QMzymeRegion {self.name} has a total charge of {chrg}.")


    def combine(self, other, name = ''):
        """
        Combine QMzymeRegion with another QMzymeRegion. 
        Duplicates are not retained.

        Parameters
        -----------
        other : QMzymeRegion
        name : Name of new QMzymeRegion.

        Returns
        ---------
        QMzymeRegion
            Combined QMzymeRegion
        """
        combined_atoms = copy.copy(self.atoms)
        for atom in other.atoms:
            if not atom.is_within(self):
                combined_atoms.append(atom)
        combined_region = QMzymeRegion(name=name, atoms=combined_atoms)
        return combined_region
    
    def subtract(self, other, name = ''):
        """
        Subtract components found in other QMzymeRegion. 

        Parameters
        -----------
        other : QMzymeRegion
        name : Name of new QMzymeRegion.

        Returns
        ---------
        QMzymeRegion
            Subtracted QMzymeRegion
        """
        atoms = []
        for atom in self.atoms:
            if not atom.is_within(other):
                atoms.append(atom)
        region = QMzymeRegion(name=name, atoms=atoms)
        return region
    
    def get_overlap(self, other):
        atoms = []
        for atom in self.atoms:
            if atom.is_within(other):
                atoms.append(atom)
        return atoms

    
    def guess_bonds():
        """
        Method under development.
        """
        pass


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



