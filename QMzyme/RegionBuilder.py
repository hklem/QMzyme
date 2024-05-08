"""
Concrete builder class to construct QMzymeRegion.
"""

import warnings
import copy
from QMzyme.QMzymeRegion import QMzymeRegion
from QMzyme.QMzymeAtom import QMzymeAtom
from MDAnalysis.core.groups import Atom
from MDAnalysis.core.groups import AtomGroup

from typing import TYPE_CHECKING, Any, Dict, Generic, List, Optional, TypeVar, Union

_QMzymeAtom = TypeVar("_QMzymeAtom", bound="QMzymeAtom")
_QMzymeRegion = TypeVar("_QMzymeRegion", bound="QMzymeRegion")
_MDAtom = TypeVar("_MDAtom", bound="Atom")
_AtomGroup = TypeVar("_AtomGroup", bound="AtomGroup")
_Atom = TypeVar("_Atom", bound=Union["QMzymeAtom", "Atom"])

remove_mda_atom_props = ['chainID', 'level', 'universe', 'bfactor', 'altLoc', 'ix_array', 'segment', 'segindex']

class RegionBuilder:

    def __init__(self, name, atom_group = None):
        self.name = name
        #self.region = QMzymeRegion(name, atoms=[])
        #self.region = None
        self.atoms = []
        self.atom_group = atom_group
        if atom_group is not None:
            self.init_atom_group(atom_group)


    def __repr__(self):
        return f"<RegionBuilder: Current QMzymeRegion, {self.name}, "+\
               f"contains {self.n_atoms} atom(s) and {self.n_residues} residue(s)>"


    @property
    def n_residues(self):
        return len(list(set(atom.resid for atom in self.atoms)))
    

    @property
    def n_atoms(self):
        return len(self.atoms)


    def init_atom_group(self, atom_group):
        for atom in atom_group.atoms:
            self.init_atom(atom)
        self.atom_group = atom_group
        #self.region.set_atom_group(atom_group)
        

    def init_atom(self, atom):
        warnings.filterwarnings('ignore')
        atom_props = self.get_atom_properties(atom)
        atom = QMzymeAtom(**atom_props)
        self.uniquify_atom(atom)
        self.atoms.append(atom)
        #self.region.add_atom(atom)

    
    def uniquify_atom(self, atom):
        temp_region = QMzymeRegion(self.name, self.atoms)
        while atom.id in temp_region.ids:
            atom.id += 1
        if atom.resid in temp_region.resids:
            residue_atoms = temp_region.get_residue(atom.resid).atoms
            atom_names = [a.name for a in residue_atoms]
            name = atom.name
            if name in atom_names:
                i = 0
                while name in atom_names:
                    i += 1
                    name = f"{atom.element}{i}"   
                    atom.set_name(name)
        return atom


    def get_region(self):
        self.region = QMzymeRegion(self.name, self.atoms)
        self.region.set_atom_group(self.atom_group)
        return self.region


    def get_atom_properties(self, atom: _Atom):
        atom_attr_dict = {}
        for attr in dir(atom):
            if attr.startswith('_') or attr.startswith('get'):
                continue
            try:
                atom_attr_dict[attr] = getattr(atom, attr)
            except:
                pass
        if isinstance(atom, Atom):
            atom_attr_dict['chain'] = atom_attr_dict['chainID']
            d = copy.copy(atom_attr_dict)
            for attr in d:
                if attr in remove_mda_atom_props:
                    del atom_attr_dict[attr]
        try:
            del atom_attr_dict['region']
        except:
            pass
        return atom_attr_dict




