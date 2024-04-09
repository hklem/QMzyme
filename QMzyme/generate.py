###############################################################################
# Code written by Heidi Klem while at
# Colorado State University as a graduate student
# in the Paton and McCullagh groups and at the
# National Institute of Standards and Technology
# as an NRC Postdoc (Fed).
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

'''Generate QM-based enzyme model.'''

import numpy as np
import json
import inspect
from typing import Optional
from datetime import datetime
import os
from QMzyme.calculate import *
from QMzyme.Biopython.Structure import Structure
from QMzyme.BiopythonWrapper import BiopythonWrapper
from QMzyme.BiopythonWrapper import res_dispatch
from QMzyme import utils
from QMzyme.utils import(
    download_pdb,
    get_coords,
    get_atoms,
    get_outlines,
    to_dict,
    filename_format,
    coords_from_pdb,
    res_charges,
    )

protein_residues = ['ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'CYM', 'CYS', 'CYX',
                    'GLH', 'GLN', 'GLU', 'GLY', 'HIS', 'HID', 'HIE', 'HIP',
                    'HYP', 'ILE', 'LEU', 'LYN', 'LYS', 'MET', 'PHE', 'PRO',
                    'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HSE', 'HSD', 'HSP',
                    'SEC', 'PYL']
positive_residues = ['HIP', 'LYS', 'ARG']
negative_residues = ['ASP', 'GLU']
solvent_list = ['HOH', 'WAT', 'T3P', 'SOL']

elements = ['H','He','Li','Be','B','C','N','O','F','Ne',
           'Na','Mg','Al','Si','P','S','Cl','Ar','K', 'Ca',
           'Sc', 'Ti', 'V','Cr', 'Mn', 'Fe', 'Co', 'Ni',
           'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
           'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
           'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
           'I', 'Xe','Cs', 'Ba','La', 'Ce', 'Pr', 'Nd', 'Pm',
           'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm',
           'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir',
           'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
           'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am',
           'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
           'Rf', 'Db', 'Sg', 'Bh','Hs', 'Mt', 'Ds', 'Rg', 'Cn',
           'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

class GenerateModel(Structure):

    def __init__(self, structure, id=None, child=0):
        if type(structure) == str:
            filename = structure
            structure = BiopythonWrapper.load_structure(structure, id)
            setattr(structure, 'pdb_file', filename)
        self.__dict__ = structure.__dict__.copy()

        if not hasattr(self, 'base'):
            self.base = self.child_list[child]
        self.models = []
        self.QMzyme_details = {}
        func = inspect.currentframe().f_code.co_name
        self.QMzyme_details = utils.record_execution(self.QMzyme_details, func)
        setattr(self, 'id', id)
        delattr(self, 'xtra')


    def __repr__(self):
        return f"<QMzyme Structure id={self._id}>"
    

    def set_catalytic_center(
        self, resname: Optional[str] = None, resnumber: Optional[int] = None, 
        chain: Optional[str] = None, overwrite=True
    ):        
        '''
        Function to define the center of the QMzyme model. This is
        typically the ligand/substrate. Currently only supports the
        definition of full residues. To create a catalytic center with multiple
        residues run set_catalytic_center again with the new residue info
        and set 'overwrite=False'. This will append the new residue information
        to the existing catalytic_center attribute. 

        Parameters
        ----------
        resname : str, optional
            Three letter code matching residue name. 
        resnumber : int, optional
            Number matching residue id/seq number.
        chain : str, optional
            Single letter matching residue chain. May be necessary to 
            uniquely identify a residue.  
        overwrite : bool, required, default=True
            To clear any existing catalytic_center definition. Set to 
            False if you would like to append current catalytic_center.

        Notes
        -----

        '''
        residues = []
        _list, _dict = ['resname', 'resnumber', 'chain'], {}
        for key, val in locals().items():
            if val != None and key in _list:
                _dict[key] = val

        for res in self.base.get_residues(): 
            if False not in [val == res_dispatch[key](res) for key, val in _dict.items()]:
                residues.append(res)
                
        if overwrite is True:
            self.catalytic_center = []
            self.catalytic_center = residues
        if overwrite is False:
            self.catalytic_center += residues


    def within_distance(self, distance_cutoff, store_model=True):
            residues = []
            center_coords = []
            for res in self.catalytic_center:
                for atom in res.get_atoms():
                    center_coords.append(atom.get_coord())
            neighbors = BiopythonWrapper.get_neighbors(self.base, center_coords, distance_cutoff)
            residues = [atom.get_parent() for atom in neighbors]
            residues = list(set(residues))
            residues = BiopythonWrapper.order_residues(residues)
            method = {}
            method['type'] = 'within_distance'
            method['cutoff'] = distance_cutoff
            method['catalytic_center'] = self.catalytic_center
            if store_model is True:
                self.store_model(residues, method)
            else:
                return BiopythonWrapper.init_model(self, residues, method)


    def truncate(self, scheme = 'CA_terminal'):
            for res in self.child_list[-1].get_residues():
                if res.resname not in protein_residues:
                    continue
                if BiopythonWrapper.has_Nterm_neighbor(res) is False and res.resname != 'PRO':
                    BiopythonWrapper.cap_terminus(res, 'N')
                if BiopythonWrapper.has_Cterm_neighbor(res) is False:
                    BiopythonWrapper.cap_terminus(res, 'C')
            self.models[-1] = self.child_list[-1]


    def store_model(self, residues, method):
        m = BiopythonWrapper.init_model(self, residues, method)
        self.add(m)
        if self.child_list[0] == self.base:
            del self.child_list[0]
            del self.child_dict[0]
        self.models.append(self.child_list[-1])
        print(f"Model {self.child_list[-1].id} created by the {method['type']} method has been stored to {self.__repr__()}.")


    def write_pdb(self, entity, filename=""):
        entity.pdb_file = BiopythonWrapper.write_pdb(entity, filename)


    def calculateQM(self, model=None, functional=None, 
                basis_set=None, opt=True, freq=True, freeze_atoms = [], 
                charge = 0, mult = 1, mem='32GB', 
                nprocs=16, program='gaussian', suffix=''):
        if model == None:
            model = self.child_list[-1]
        if charge != None:
            model.charge = charge
        if mult != None:
            model.multiplicity = mult
        CalculateQM(model, functional, basis_set, opt, freq, freeze_atoms, charge, mult, mem, nprocs, program, suffix)


    def write_json(self, filename=None):
        _dict = {}
        for model in self.models:
            _dict['Model '+str(model.id)] = BiopythonWrapper.make_model_dict(model)
        QMzyme_dict = {self.__repr__()[1:-1]: _dict}
        if filename == None:
            filename = self.id+'.json'
        with open(filename, 'w') as f:
            json.dump(QMzyme_dict, f, indent=4, sort_keys=True)
