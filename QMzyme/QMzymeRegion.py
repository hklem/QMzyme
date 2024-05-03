"""
Product of the concrete builder class RegionBuilder.
"""

import copy
import numpy as np
from typing import TYPE_CHECKING, Any, Dict, Generic, List, Optional, TypeVar
from QMzyme.QMzymeAtom import QMzymeAtom
from MDAnalysis.core.groups import AtomGroup

_QMzymeAtom = TypeVar("_QMzymeAtom", bound="QMzymeAtom")
_AtomGroup = TypeVar("_AtomGroup", bound="AtomGroup")

class QMzymeRegion:
    #def __init__(self, name, atoms: list[_QMzymeAtom], residues, atom_group: _AtomGroup):
    def __init__(self, name):
        self.name = name
        self.atoms = []
        #self.residues = []

    def __repr__(self):
        return f"<QMzymeRegion {self.name} contains {self.n_atoms} atoms and {self.n_residues} residues>"
    
    def _is_unique(self, atom):
        if atom in self.atoms:
            return False
        else:
            return True
        
    # # Do I think an atoms property like this so it can be updated without re-instantiating QMzmeRegion?    
    # @property
    # def atoms(self):
    #     return self._atoms
    
    @property
    def atom_group(self):
        return self.__atom_group
        
    @property
    def n_atoms(self):
        return len(self.atoms)
    
    @property
    def n_residues(self):
        #self.residues = list(set([atom.resid for atom in self.atoms]))
        return len(self.resids)
    
    @property
    def ids(self):
        return [atom.id for atom in self.atoms]
    
    @property
    def resids(self):
        return list(set([atom.resid for atom in self.atoms]))
    
    def set_atom_group(self, atom_group):
        self.__atom_group = atom_group
    
    def sort(self, key='id'):
        sorted_atoms = []
        original_atoms = copy.copy(self.atoms)
        for i in range(self.n_atoms):
            key_list = [getattr(atom, key) for atom in original_atoms]
            min = np.min(key_list)
            min_index = key_list.index(min)
            sorted_atoms.append(original_atoms[min_index])
            del original_atoms[min_index]
        self.atoms = sorted_atoms

    def add_attr(self, attr_name, attr_value):
        setattr(self, attr_name, attr_value)
    
    def get_AtomGroup(self):
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
        self.atoms.append(atom)

    def get_residue_atoms(self, resid):
        return [atom for atom in self.atoms if atom.resid == resid]
    
    def get_resname(self, resid):
        return self.get_residue_atoms(resid)[0].resname

    def get_residue_atom(self, resid, atom_name):
        for atom in self.get_residue_atoms(resid):
            if atom.name == atom_name:
                return atom


