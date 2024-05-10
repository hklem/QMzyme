'''
Module in charge of creating input files for QM-only or QM/MM calculations. This 
module integrates the `AQME QPREP <https://aqme.readthedocs.io/en/latest/API/aqme.qprep.html>`_ 
workflow.

Notes
...............
    *   Currently optimized to generate QM-only Gaussian input files.
'''

import copy
import os
from QMzyme.aqme.qprep import qprep


# class MultiscaleCalculation:
#     """
#     Class to prepare input for a multiscale calculateion (i.e., QM1/QM2 or QM/MM).
#     """
#     def __init__(self, low_layer, high_layer):
#         self.set_layers()
#         for layer, region in self.layers.items():
#             if 'QM' in layer:
#                 # prepare QM calc
#                 self.write_qm_input(region)

#     def set_layers(self):
#         self.layers = {}
#         for region in self.regions:
#             while region.layer in self.layers:
#                 count = int(region.layer[-1])
#                 region.layer = f"{region.layer[:2]}{count+1}"
#             self.layers[region.layer] = region


class QM_Region(object):
    """
    Class to prepare a QMzymeRegion for a QM calculation.
    """
    def __init__(self, region):
        """
        :param region: QMzyme region to treat at the QM level.
        :type region: QMzymeRegion
        """
        region.layer = 'QM1'
        self.atoms = copy.copy(region.atoms)
        self.region = copy.copy(region)
        self.qm_input = ''

    def set_basis_set(self, basis_set: str, overwrite=False):
        """
        Method to assign basis set to region.
        
        :param basis_set: Basis set to assign
        :type basis_set: str, required
        :param overwrite: Option to assign a different basis_set once one has already been assigned.
        :type overwrite: bool, default=False
        """
        if hasattr(self, 'basis_set') and overwrite is False:
            raise UserWarning(f"QM_Region already has basis_set defined as "+
                              f"{self.basis_set}. To overwrite, set overwrite=False. If "+
                              f"you want to use two basis sets, see set_gen_basis_set().")
        else:
            self.basis_set = basis_set
            for atom in self.atoms:
                self._assign_bs_to_atom(basis_set, atom)
    
    def _assign_bs_to_atom(self, basis_set, atom):
        atom.basis_set = basis_set

    def set_gen_atoms(self, elements: list, basis_set: str):
        """
        Method to define what elements to assign a 2nd basis set to.

        :param elements: A list of element names for atoms to be treated with the gen(ECP).
        :type elements: list[str], required
        :param basis_set: Basis set used for gen(ECP) atoms (i.e. 'def2svp')
        :type basis_set: str, required
        """
        elements = ' '.join(elements)
        sel = 'element ' + elements
        atom_group = self.region.convert_to_AtomGroup().select_atoms(sel)
        for atom in atom_group.atoms:
            qmz_atom = self.region.get_atom(atom.id)
            self.assign_bs_to_atom(basis_set, atom)
        self.get_atoms = elements
        self.bs_nogen = self.basis_set
        
    def set_functional(self, functional: str):
        """
        Method to assign functional.
        :param functional: The QM level of theory to use.
        :param type: str, required
        """
        self.functional = functional
        for atom in self.atoms:
            atom.functional = self.functional
    
    def set_charge(self, charge: int):
        """
        Method to assign charge to region.
        :param charge: Cumulative atomic charge of the region.
        :type charge: int, required
        """
        self.charge = charge

    def set_multiplicity(self, mult: int):
        """
        Method to assign multiplicity to region.
        :param mult: Multiplicity of the region.
        :type mult: int, required
        """
        self.mult = mult

    def set_qm_input(self, input: str):
        """
        Keywords line for new input files (i.e. 'opt freq').
        The functional and basis set will be included automatically after 
        being set. 
        """
        self.qm_input = input

    def set_qm_end(self, input: str):
        """
        Final line(s) in the new input files.
        """
        self.qm_end = input

    def set_chk(self, value: bool):
        """
        Include the chk input line in input file for Gaussian calculations.
        """
        self.chk = value

    def set_chk_path(self, value: str):
        """
        PATH to store CHK files. For example, if chk_path='root/user, the 
        chk line of the input file would be %chk=root/user/FILENAME.chk.
        """
        self.chk_path = value

    def set_memory(self, value: str):
        """
        Memory for the QM calculations (i) Gaussian: total memory; 
        (ii) ORCA: memory per processor. Must include units.
        """
        self.mem = value

    def set_n_processors(self, value: int):
        """
        Number of processors used in the QM calculations.
        """
        self.n_procs = value

    def check_missing_attr(self, *args):
        for attr in args:
            if not hasattr(self, attr):
                raise UserWarning(f"The current region is missing {attr} information. Run set_{attr}() to fulfill.")
            
    def write_qm_input(self, program, filename=None):
        # check necessary args have been set
        self.check_missing_attr('basis_set', 'functional', 'charge', 'mult')

        # combine basis_set and functional into qm_input
        if self.basis_set not in self.qm_input:
            self.qm_input = f"{self.basis_set} {self.qm_input}"
        if self.functional not in self.qm_input:
            self.qm_input = f"{self.functional} {self.qm_input}"
        
        # check if any atoms have is_fixed=True
        self.freeze = self.region.get_indices('is_fixed', True)

        # create pdb file that will be read by aqme qprep
        #we could also just input the atoms and coords separately
        file = self.region.write(filename)
        self.files = file

        # object attributes that don't correspeond to aqme qprep keywords.
        exclude = ['region', 'atoms', 'basis_set', 'functional']
        calc_info = self.__dict__
        for info in exclude:
            del calc_info[info]
        
        # use qprep to write input file
        qprep(program='gaussian', **calc_info)

        # clean up file name because AQME QPREP adds _conf ending.
        calc_file = f"./QCALC/{file.split('.pdb')[0]}.com"
        try:
            os.rename(f"./QCALC/{file.split('.pdb')[0]}_conf_1.com", calc_file)
        except:
            pass
            
        return calc_info


# class MM_Region:
#     def __init__(self, region):
#         self.region = region
#         self.atoms = region.atoms
#         self.layer = 'MM1'
#         self.force_field = []

#     def set_program(self, name: str):
#         self.program = name
    
#     def set_charge(self, charge: int):
#         self.charge = charge

#     def set_multiplicity(self, mult: int):
#         self.multiplicity = mult

#     def set_force_field(self, force_field: str, selection='all'):
#         self.force_field.append(force_field)
#         ids = self.region.ids
#         if selection != 'all':
#             atom_group = self.region.convert_to_AtomGroup().select_atoms(selection)
#             ids = atom_group.ids
#         for id in ids:
#             atom = self.region.get_atom(id)
#             atom.force_field = force_field

