'''
Module in charge of create input files for QM-only or QM/MM calculations. This 
module integrates the `AQME QPREP <https://aqme.readthedocs.io/en/latest/API/aqme.qprep.html>`_ 
workflow.

Prerequisites
...............

    *   

Notes
...............

    *   Currently optimized to generate QM-only Gaussian input files.
'''

import copy
import os
from QMzyme.aqme.qprep import qprep
import QMzyme.MDAnalysisWrapper as MDAWrapper
from QMzyme.RegionBuilder import RegionBuilder


class CalculateModel:
    def __init__(self, model):
        self.set_layers()
        for layer, region in self.layers.items():
            if 'QM' in layer:
                # prepare QM calc
                self.write_qm_input(region)


    def set_layers(self):
        self.layers = {}
        for region in self.regions:
            while region.layer in self.layers:
                count = int(region.layer[-1])
                region.layer = f"{region.layer[:2]}{count+1}"
            self.layers[region.layer] = region



class QM_Region(object):
    def __init__(self, region):
        region.layer = 'QM1'
        self.atoms = copy.copy(region.atoms)
        self.region = copy.copy(region)
        self.qm_input = ''

    def set_basis_set(self, basis_set: str, selection='all', overwrite=False):
        # self.basis_set.append(basis_set)
        # ids = self.region.ids
        # if selection != 'all':
        #     atom_group = self.region.convert_to_AtomGroup().select_atoms(selection)
        #     ids = atom_group.ids
        #     if len(self.basis_set) > 1:
        #         self.set_gen_atoms([])
        #         self.set_gen_basis_set(basis_set)
        # for id in ids:
        #     atom = self.region.get_atom(id)
        #     atom.basis_set = basis_set
        if hasattr(self, 'basis_set') and overwrite is False:
            raise UserWarning(f"QM_Region already has basis_set defined as "+
                              f"{self.basis_set}. To overwrite, set overwrite=False. If "+
                              f"you want to use two basis sets, see set_gen_basis_set().")
        else:
            self.basis_set = basis_set
            for atom in self.atoms:
                self.assign_bs_to_atom(basis_set, atom)
    
    def assign_bs_to_atom(self, basis_set, atom):
        atom.basis_set = basis_set

    def set_gen_atoms(self, elements: list, basis_set: str):
        """
        Accepts atom elements as input.
        """
        elements = ' '.join(elements)
        sel = 'element ' + elements
        atom_group = self.region.convert_to_AtomGroup().select_atoms(sel)
        for atom in atom_group.atoms:
            qmz_atom = self.region.get_atom(atom.id)
            self.assign_bs_to_atom(basis_set, atom)
        self.bs_nogen = self.basis_set
        
    def set_functional(self, functional: str):
        # self.functional.append(functional)
        # ids = self.region.ids
        # if selection != 'all':
        #     atom_group = self.region.convert_to_AtomGroup().select_atoms(selection)
        #     ids = atom_group.ids
        # for id in ids:
        #     atom = self.region.get_atom(id)
        #     atom.functional = functional
        self.functional = functional
        for atom in self.atoms:
            atom.functional = self.functional
    
    def set_charge(self, charge: int):
        self.charge = charge

    def set_multiplicity(self, mult: int):
        self.mult = mult

    def set_qm_input(self, input: str):
        self.qm_input = input

    def set_qm_end(self, input: str):
        self.qm_end = input

    def set_chk(self, value: bool):
        self.chk = value

    def set_chk_path(self, value: str):
        self.chk_path = value

    def set_memory(self, value: str):
        self.mem = value

    def set_n_processors(self, value: int):
        self.n_procs = value
            
    # def check(self):
    #     """
    #     Method to ensure every atom has been assigned a basis_set and functional, and other properties have been assigned.
    #     """
    #     self.region.check_missing_attr('basis_set')
    #     self.region.check_missing_attr('functional')
    #     self.check_missing_attr('charge', 'multiplicity')
    #     print("All checks passed.")
    
    def check_missing_attr(self, *args):
        for attr in args:
            if not hasattr(self, attr):
                raise UserWarning(f"The current region is missing {attr} information. Run set_{attr}() to fulfill.")
            
    def write_qm_input(self, program, filename=None):
        # if self.basis_set[0] not in qm_input:
        #     qm_input = f"{self.basis_set[0]} {qm_input}"
        # if self.functional[0] not in qm_input:
        #     qm_input = f"{self.functional[0]} {qm_input}"
        if self.basis_set not in self.qm_input:
            self.qm_input = f"{self.basis_set} {self.qm_input}"
        if self.functional not in self.qm_input:
            self.qm_input = f"{self.functional} {self.qm_input}"
        
        self.freeze = self.region.get_indices('is_fixed', True)
        file = self.region.write(filename)
        self.files = file

        exclude = ['region', 'atoms', 'basis_set', 'functional']
        calc_info = self.__dict__
        for info in exclude:
            del calc_info[info]
        qprep(program='gaussian', **calc_info)

        # qprep(
        #     files = file,
        #     charge = self.charge,
        #     mult = self.multiplicity,
        #     qm_input = qm_input,
        #     freeze = freeze,
        #     program = program,
        #     mem = mem,
        #     nprocs = nprocs,
        #     suffix = None
        # )

        # clean up file name because AQME QPREP adds _conf ending.
        calc_file = f"./QCALC/{file.split('.pdb')[0]}.com"
        try:
            os.rename(f"./QCALC/{file.split('.pdb')[0]}_conf_1.com", calc_file)
        except:
            pass
            
        return calc_info


class MM_Region:
    def __init__(self, region):
        self.region = region
        self.atoms = region.atoms
        self.layer = 'MM1'
        self.force_field = []

    def set_program(self, name: str):
        self.program = name
    
    def set_charge(self, charge: int):
        self.charge = charge

    def set_multiplicity(self, mult: int):
        self.multiplicity = mult

    def set_force_field(self, force_field: str, selection='all'):
        self.force_field.append(force_field)
        ids = self.region.ids
        if selection != 'all':
            atom_group = self.region.convert_to_AtomGroup().select_atoms(selection)
            ids = atom_group.ids
        for id in ids:
            atom = self.region.get_atom(id)
            atom.force_field = force_field

