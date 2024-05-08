"""
Product of the builder class RegionBuilder.
"""

import copy
import numpy as np
from typing import TYPE_CHECKING, Any, Dict, Generic, List, Optional, TypeVar
from QMzyme.QMzymeAtom import QMzymeAtom
from MDAnalysis.core.groups import AtomGroup
from QMzyme import MDAnalysisWrapper as MDAwrapper

_QMzymeAtom = TypeVar("_QMzymeAtom", bound="QMzymeAtom")
_AtomGroup = TypeVar("_AtomGroup", bound="AtomGroup")

class QMzymeRegion:
    #def __init__(self, name, atoms: list[_QMzymeAtom], residues, atom_group: _AtomGroup):
    # def __init__(self, name):
    def __init__(self, name, atoms: list, atom_group= None):
        self.name = name
        self.atoms = atoms
        self.atom_group = atom_group

    def __repr__(self):
        return f"<QMzymeRegion {self.name} contains {self.n_atoms} atom(s) and {self.n_residues} residue(s)>"
    
    # @property
    # def atoms(self):
    #     return self._atoms
    
    # @atoms.setter
    # def atoms(self, value):
    #     self.atoms = value
        
    # @property
    # def atom_group(self):
    #     return self.atom_group
        
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
            res = QMzymeResidue(resname, resid, atoms)
            residues.append(res)
        return residues

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
        
    # def uniquify_atom(self, atom):
    #     atom = copy.copy(atom)
    #     if self.atoms == None:
    #         return atom
    #     ids = self.ids
    #     while atom.id in ids:
    #         atom.id += 1
    #     if atom.resid in self.resids:
    #         residue_atoms = self.get_residue(atom.resid).atoms
    #         atom_names = [a.name for a in residue_atoms]
    #         name = atom.name
    #         if name in atom_names:
    #             i = 0
    #             while name in atom_names:
    #                 i += 1
    #                 name = f"{atom.element}{i}"
    #             atom.set_name(name)
    #     return atom

    
    def get_residue(self, resid):
        for res in self.residues:
            if res.resid == resid:
                return res


    def write(self, filename=None):
        # Housekeeping
        if filename is None:
            filename = f"{'_'.join(self.name.split(' '))}.pdb"
        ag = self.convert_to_AtomGroup()
        ag.write(filename)

    def convert_to_AtomGroup(self):
        return MDAwrapper.build_universe_from_QMzymeRegion(self)
    
    def set_fixed_atoms(self, ids: list):
        for id in ids:
            atom = self.get_atom(id)
            atom.set_fixed()


class QMzymeResidue(QMzymeRegion):
    def __init__(self, resname, resid, atoms, chain=None):
        self.resname = resname
        self.resid = resid
        self.atoms = atoms
        if chain is None:
            chain = self.atoms[0].get_chain()
        self.chain = chain

    # def set_atoms(self, atoms):
    #     return [atom for atom in atoms]

    # @property
    # def atoms(self):
    #     return self.__atoms

    # @atoms.setter
    # def atoms(self, value):
    #     self.__atoms = value

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



