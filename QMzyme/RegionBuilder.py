###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
Class to build a QMzymeRegion object.
"""

import warnings
import copy
from QMzyme.QMzymeRegion import QMzymeRegion
from QMzyme.QMzymeAtom import QMzymeAtom
from QMzyme.converters import mda_atom_to_qmz_atom
from MDAnalysis.core.groups import Atom
from MDAnalysis.core.groups import AtomGroup

from typing import TYPE_CHECKING, Any, Dict, Generic, List, Optional, TypeVar, Union

_QMzymeAtom = TypeVar("_QMzymeAtom", bound="QMzymeAtom")
_QMzymeRegion = TypeVar("_QMzymeRegion", bound="QMzymeRegion")
_MDAtom = TypeVar("_MDAtom", bound="Atom")
_AtomGroup = TypeVar("_AtomGroup", bound="AtomGroup")
_Atom = TypeVar("_Atom", bound=Union["QMzymeAtom", "Atom"])

#remove_mda_atom_props = ['chainID', 'level', 'universe', 'bfactor', 'altLoc', 'ix_array', 'segment', 'segindex']
remove_mda_atom_props = ['level', 'universe', 'bfactor', 'altLoc', 'ix_array', 'segment', 'segindex']

class RegionBuilder:

    def __init__(self, name, atom_group = None, universe = None):
        self.name = name
        self.atoms = []
        self.region = None
        self.atom_group = atom_group
        self.universe = universe
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
        """
        It is assumed that the atoms in the selection are already unique. 
        """
        for atom in atom_group.atoms:
            self.init_atom(atom, uniquify=True)
        self.atom_group = atom_group
        self.universe = atom_group.universe
        region = self.get_region()
        #return region

    def init_atom(self, atom, uniquify=True):
        warnings.filterwarnings('ignore')
        if isinstance(atom, Atom):
            atom = mda_atom_to_qmz_atom(atom)
        else:
            atom_props = self.get_atom_properties(atom)
            if uniquify is True:
                self.uniquify_atom(atom_props)
            atom = QMzymeAtom(**atom_props)
        # if self.atoms != []:
        #     print(atom, self.atoms[-1])
        self.atoms.append(atom)

    def uniquify_atom(self, atom_props):
        temp_region = QMzymeRegion(self.name, self.atoms)
        resid = atom_props['resid']
        element = atom_props['element']
        if resid in temp_region.resids:
            residue_atoms = temp_region.get_residue(resid).atoms
            names = [a.name for a in residue_atoms]
            i = 0
            while atom_props['name'] in names:
                i += 1
                atom_props['name'] = f"{element}{i}" 
        while atom_props['id'] in temp_region.ids:
            atom_props['id'] += 1

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
            try:
                atom_attr_dict['chain'] = atom_attr_dict['chainID']
            except:
                atom_attr_dict['chain'] = 'X'
            d = copy.copy(atom_attr_dict)
            for attr in d:
                if attr in remove_mda_atom_props:
                    del atom_attr_dict[attr]
        try:
            del atom_attr_dict['region']
        except:
            pass
        return atom_attr_dict

    def get_region(self):
        self.region = QMzymeRegion(self.name, self.atoms, self.atom_group, self.universe)
        return self.region
