"""
Product of the concrete builder class RegionBuilder.
"""

import copy
import numpy as np
from typing import TYPE_CHECKING, Any, Dict, Generic, List, Optional, TypeVar
from QMzyme.QMzymeAtom import QMzymeAtom
#from QMzyme.QMzymeResidue import QMzymeResidue
from MDAnalysis.core.groups import AtomGroup
from QMzyme import MDAnalysisWrapper as MDAwrapper

_QMzymeAtom = TypeVar("_QMzymeAtom", bound="QMzymeAtom")
_AtomGroup = TypeVar("_AtomGroup", bound="AtomGroup")

class QMzymeRegion:
    #def __init__(self, name, atoms: list[_QMzymeAtom], residues, atom_group: _AtomGroup):
    def __init__(self, name):
        self.name = name
        self.atoms = []

    def __repr__(self):
        return f"<QMzymeRegion {self.name} contains {self.n_atoms} atoms and {self.n_residues} residues>"
    
    @property
    def atom_group(self):
        return self.__atom_group
        
    @property
    def ids(self):
        return [atom.id for atom in self.atoms]
    
    @property
    def resids(self):
        return list(set([atom.resid for atom in self.atoms]))
    
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
            chain = atoms[0]._get_chain()
            res = QMzymeResidue(resname, resid, atoms, chain)
            residues.append(res)
        return residues
    
    
    def set_atom_group(self, atom_group):
        self.__atom_group = atom_group
    
    def sort(self, key='id', in_place=False):
        sorted_atoms = []
        original_atoms = copy.copy(self.atoms)
        for i in range(self.n_atoms):
            key_list = [getattr(atom, key) for atom in original_atoms]
            min = np.min(key_list)
            min_index = key_list.index(min)
            sorted_atoms.append(original_atoms[min_index])
            del original_atoms[min_index]
        if in_place is True:
            self.atoms = sorted_atoms
        else:
            return sorted_atoms

    def add_attr(self, attr_name, attr_value):
        setattr(self, attr_name, attr_value)
    
    def get_atom_group(self):
        return self.__atom_group
    
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
        atom = self.uniquify_atom(atom)
        self.atoms.append(atom)
        
    def uniquify_atom(self, atom):
        while atom.id in self.ids:
            atom.id += 1
        if atom.resid in self.resids:
            residue_atoms = self.get_residue(atom.resid).atoms
            i = 1
            while atom.name in [a.name for a in residue_atoms]:
                atom.name = f"{atom.element}{i}"
                i += 1
        return atom

    
    def get_residue(self, resid):
        for res in self.residues:
            if res.resid == resid:
                return res


    def write_file(self, filename=None):
        # Housekeeping
        if filename is None:
            filename = f"{'_'.join(self.name.split(' '))}.pdb"
        ag = self.convert_to_AtomGroup()
        ag.write(filename)

    def convert_to_AtomGroup(self):
        return MDAwrapper.build_universe_from_QMzymeRegion(self)


class QMzymeResidue(QMzymeRegion):
    def __init__(self, resname, resid, atoms, chain=None):
        self.resname = resname
        self.resid = resid
        self.atoms = atoms
        if chain is not None:
            self.chain = chain
        else:
            self.set_chain(self.atoms[0]._get_chain())

    def __repr__(self):
        return f"<QMzymeResidue resname: {self.resname}, resid: {self.resid}, chain: {self.chain}>"
    
    def get_atom(self, atom_name):
        for atom in self.atoms:
            if atom.name == atom_name:
                return atom

    def set_chain(self, value: str):
        self.set_chain(value)



