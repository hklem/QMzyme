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

    def __init__(self, name, atom_group: _AtomGroup=None):
        self.atom = None
        #self.residue = None
        self.region = QMzymeRegion(name)
        if atom_group is not None:
            self.region.set_atom_group(atom_group)
            for atom in atom_group.atoms:
                self.init_atom(atom)

    def init_atom(self, atom):
        warnings.filterwarnings('ignore')
        atom_props = self.get_atom_properties(atom)
        #print(atom_props)
        self.atom = QMzymeAtom(region=self.region, **atom_props)
        self.region.add_atom(self.atom)

    def get_region(self):
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
            #del atom_attr_dict['altloc']
        try:
            del atom_attr_dict['region']
        except:
            pass
        return atom_attr_dict




