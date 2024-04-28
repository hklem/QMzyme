###############################################################################
# Code written by Heidi Klem while at
# Colorado State University as a graduate student
# in the Paton and McCullagh groups and at the
# National Institute of Standards and Technology
# as an NRC Postdoc (Fed).
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

'''
Module in charge of creating the QM enzyme models given a starting structure. 
This module can be used to call CalculateModel to generate QM calculation 
input files.

Required input
...............

    *   Prepared starting structure in PDB format.
    *   Residue selection for catalytic center definition.
    *   Distance cutoff for region selection.

Output options
...............

    *   QMzyme structure object containing all model information.
    *   PDB files for catalytic center or truncated models.
    *   QM calculation input files for Gaussian or Orca packages.
    *   JSON file containing all model information.
'''

import os
import numpy as np
import json
import copy
import inspect
import warnings
from QMzyme.calculate import *
from QMzyme import MDAnalysisWrapper
from QMzyme.utils import record_execution

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

class GenerateModel():

    def __init__(self, pdb_file, id=None):
        """
        Constructor method for GenerateModel class.

        :param pdb_file: PDB file name.  
        :type pdb_file: str, required
        :param id: Name to associate structure with, defaults to file name. Will be used in naming of any files produced.
        :type id: str
        """

        self.models = []
        self.QMzyme_calls = {}
        self.catalytic_center = None
        func = inspect.currentframe().f_code.co_name
        u = MDAnalysisWrapper.init_universe(pdb_file)
        if id is None:
            id = os.path.basename(pdb_file).split('.')[0]

        setattr(self, 'pdb_file', pdb_file)
        setattr(self, 'universe', u)
        record_execution(self.QMzyme_calls, func)
        setattr(self, 'id', id)

        Hs_present = False
        for atom in u.atoms:
            if atom.element == 'H':
                Hs_present = True
                break
        if Hs_present == False:
            raise Exception("PDB structure does not contain hydrogens. Please pre-process this structure. We recommend using tLeap. See Documentation.")


    def __repr__(self):
        return f"<QMzyme id={self.id}>"


    def set_catalytic_center(
        self,
        selection, 
        overwrite=True,
        save_pdb=False,
        filename=None
    ):        
        '''
        Function to define the center of the QMzyme model. This is
        typically the ligand/substrate. To add onto an already defined catalytic center 
        run set_catalytic_center again with the new selection and set 'overwrite=False'. 

        :param selection: Selection of atoms to be made- based on `MDAnalysis selection command language <https://docs.mdanalysis.org/stable/documentation_pages/selections.html>`_.
        :type selection: str, required`

        :param overwrite: To clear any existing catalytic_center definition. Set to False if you would like to append current catalytic_center. Defaults to True.
        :type overwrite: bool, default: overwrite=True

        '''
        self.catalytic_center_definition = selection
        catalytic_center = MDAnalysisWrapper.select_atoms(self.universe, selection) 

        if overwrite is True:
            setattr(self, 'catalytic_center', catalytic_center)
        if overwrite is False:
            self.catalytic_center = self.catalytic_center.concatenate(catalytic_center)
        if save_pdb is True:
            pass
        self.atoms = self.catalytic_center

        func = inspect.currentframe().f_code.co_name
        record_execution(self.QMzyme_calls, func)

        if save_pdb is True:
            if filename is None:
                filename=f'{self.id}_catalytic_center.pdb'
            MDAnalysisWrapper.write_pdb(self.atoms, filename)


    def size_scan(self, minimum_size, save_last_pdb=True, save_all_pdb=False):
        n_atoms = 0
        cutoff = 0
        print(minimum_size)
        while n_atoms < minimum_size:
            print(n_atoms, cutoff)
            cutoff +=1
            self.within_distance(self, distance_cutoff=cutoff, save_pdb=save_all_pdb)
            n_atoms = self.atoms.n_atoms
        if save_last_pdb is True:
            MDAnalysisWrapper.write_pdb(self.atoms)


    def within_distance(self, distance_cutoff, save_pdb=False, filename=None):
        '''
        Function to select all residues that have at least one atom within a
        specified distance to any atom in self.catalytic_center. By default,
        the resulting model is appended to self.models. 

        :param distance_cutoff: Numerical value specifying selection cutoff.
        :type distance_cutoff: int or float, required

        :param store_model: To specify if model should be recorded in QMzyme object. Set to False if you are just messing around with this function and do not want everything recorded. Defaults to True.
        :type store_model: bool

        .. code-block:: text

            Notes:

            The region will not be spherically selected unless your catalytic_center definition 
            only contains one atom. Therefore, the shape of the selection region depends
            on the shape of the catalytic_center.
        '''

        neighbors = MDAnalysisWrapper.get_neighbors(
            self.universe.select_atoms('all'), 
            self.catalytic_center, 
            distance_cutoff
        )

        self.neighbors = neighbors

        residues = [atom.residue for atom in neighbors]
        residues = MDAnalysisWrapper.order_residues(residues)
        atoms = []

        method = {'type': inspect.currentframe().f_code.co_name,
                  'cutoff': distance_cutoff, 
                  'catalytic_center': self.catalytic_center_definition,
                  'neighbors': neighbors,
                  'residues': residues}
        
        self.residues = residues
        for res in residues:
            atoms+=[res.atoms]
        self.atoms = sum(atoms)

        func = inspect.currentframe().f_code.co_name
        record_execution(self.QMzyme_calls, func)

        if save_pdb is True:
            if filename is None:
                filename=f'{self.id}_cutoff_{distance_cutoff}.pdb'
            MDAnalysisWrapper.write_pdb(self.atoms, filename)


    def truncate(self, residues=None, scheme='CA_terminal', filename='truncated_model.pdb', store_result=True):
        '''
        Function to remove extraneous atoms in preparation for model calculation. The added hydrogens 
        will have a bond length of 1.09 Angstroms (equilibrium CH bond lenght), along the bond vector of the original atom the H is replacing.

        :param scheme: See `QMzyme Documentation <https://hklem-qmzyme-documentation.readthedocs.io>`_ for explanations of each truncation scheme. Defaults to 'CA_terminal'.
        :type scheme: str

        .. code-block:: text

            Notes:

            Currently, the only available scheme is CA_terminal, which will remove 
            backbone atoms only on residues without their sequence neighbor present in
            the model, and cap the C-alpha atom with hydrogen(s). I.e., if the following
            resnumbers are in the model: [15, 23, 24, 25], then the N-terminal and C-terminal 
            backbone atoms will be removed from residue 15 and the C-alpha will be converted 
            to a methyl group, the N-terminal backbone atoms of residue 23 will be removed, but
            the C-terminal backbone atoms will remain, as will the N-terminal backbone atoms of
            residues 24 and 25, but the C-terminal backbone atoms of 25 will be removed. 
            
            Additional schemes will be created in the future.

        '''
        self.removed_sidechain = []
        self.removed_Nterm = []
        self.removed_Cterm = []
        self.truncated_model = []

        #Need to init new universe so original one isn't overwritten
        u = MDAnalysisWrapper.init_universe(self.pdb_file)
        residues = [MDAnalysisWrapper.get_parallel_residue(res, u) for res in self.residues]
        neighbors = [MDAnalysisWrapper.get_parallel_atom(atom, u) for atom in self.neighbors]

        for i, res in enumerate(residues):
            if res.resname not in protein_residues:
                self.truncated_model.append(res.atoms)
                continue
            side_chain_atoms = []
            keep_sc = False
            #Loop over side chain atoms. If none are within cutoff distance don't include them
            for atom in res.atoms:
                if atom.name not in ['C', 'O', 'N', 'H']:
                    side_chain_atoms.append(atom)
                    if atom in neighbors or res.resname == 'GLY':
                        keep_sc = True
            if keep_sc is True:
                self.truncated_model += side_chain_atoms
            if keep_sc is False:
                self.removed_sidechain.append(res)
                CA_atom = MDAnalysisWrapper.get_atom(res, 'CA')
                CB_atom = MDAnalysisWrapper.get_atom(res, 'CB')
                HA_atom = MDAnalysisWrapper.get_atom(res, 'HA')
                self.truncated_model.append(CA_atom)
                self.truncated_model.append(HA_atom)
                cap_atom = MDAnalysisWrapper.cap_backbone_CA(CB_atom)
                self.truncated_model.append(cap_atom)
            C_atom = MDAnalysisWrapper.get_atom(res, 'C')
            O_atom = MDAnalysisWrapper.get_atom(res, 'O')
            try:
                next_res = residues[i+1]
                next_N_atom = MDAnalysisWrapper.get_atom(next_res, 'N')
            except:
                # last residue, so cap CA at C_atom
                #MDAnalysisWrapper.cap_backbone_CA(C_atom)
                C_atom = MDAnalysisWrapper.cap_backbone_CA(C_atom)
                self.truncated_model.append(C_atom)
                continue
            # If very first residue, cap CA at N_atom... unless it's proline
            if res == self.residues[0]:
                self.removed_Nterm.append(res)
                N_atom = MDAnalysisWrapper.get_atom(res, 'N')
                if res.resname == 'PRO':
                    cap_atom = MDAnalysisWrapper.cap_backbone_N(N_atom)
                    self.truncated_model.append(cap_atom)
                else:
                    N_cap= MDAnalysisWrapper.cap_backbone_CA(N_atom)
                    self.truncated_model.append(N_cap)
            #Look at next residue in list to decide how to treat current C term and next N term
            # if next res is not in sequence, truncate current res CA at C_atom and
            # truncate next res CA at N_atom... unless it's proline
            if res.resid+1 != next_res.resid:
                self.removed_Cterm.append(res)
                next_N_atom = MDAnalysisWrapper.get_atom(next_res, 'N')
                if next_res.resname == 'PRO':
                    cap_atom = MDAnalysisWrapper.cap_backbone_N(next_N_atom)
                    self.truncated_model.append(cap_atom)
                else:
                    self.removed_Nterm.append(next_res)
                    C_cap = MDAnalysisWrapper.cap_backbone_CA(C_atom)
                    self.truncated_model.append(C_cap)
                    N_cap = MDAnalysisWrapper.cap_backbone_CA(next_N_atom)
                    self.truncated_model.append(N_cap)
            #If next res is in sequence keep all C term and next res N term backbone atoms
            else:
                self.truncated_model.append(C_atom)
                self.truncated_model.append(O_atom)
                self.truncated_model.append(next_N_atom)
                if next_res.resname != 'PRO':
                    next_H_atom = MDAnalysisWrapper.get_atom(next_res, 'H')
                    self.truncated_model.append(next_H_atom)

        self.truncated_model = sum(self.truncated_model)
        MDAnalysisWrapper.write_pdb(self.truncated_model, filename)
        func = inspect.currentframe().f_code.co_name
        record_execution(self.QMzyme_calls, func)


    def calculateQM(self, model=None, functional=None, 
                basis_set=None, opt=True, freq=True, freeze_atoms = [], 
                charge = 0, mult = 1, mem='32GB', 
                nprocs=16, program='gaussian', suffix=''):
        """
        Function to generate QM only calculation input file.

        :param model: Defaults to last model in the list.
        :type model: Model

        :param functional: QM theory functional to use, i.e., 'b3lyp'.
        :type functional: str, required

        :param basis_set: QM theory basis set to use, i.e., '6-31g(d)'. 
        :type basis_set: str, required

        :param opt: Set up input to perform optimization if True. Default is True.
        :type opt: bool

        :param freq: Set up input to perform frequency analysis if True. Default is True.
        :type freq: bool

        :param freeze_atoms: List containing the PDB style names of any atoms to be frozen during optimization. E.x., ['CA'] will assign frozen atom flags to all C-alpha atoms.
        :type freeze_atoms: list

        :param charge: Electronic charge of model. 
        :type charge: int, required
        
        :param mult: Multiplicity of model.
        :type mult: int, required
        
        :param mem: Amount of memory to specify in input file. Default is '32GB'.
        :type mem: str

        :param nprocs: Amount of processors to specify for job. Default is 16. 
        :type nprocs: int

        :param program: Name of QM program to be used. Options are 'gaussian' or 'orca'. 
        :type program: str

        :param suffix: To be appended to end of calculation input file name. Default to 'calc_1'.
        :type suffix: str

        """

        if model == None:
            model = self.truncated_model
        if charge != None:
            model.charge = charge
        if mult != None:
            model.multiplicity = mult
        CalculateQM(model, functional, basis_set, opt, freq, freeze_atoms, charge, mult, mem, nprocs, program, suffix)

        func = inspect.currentframe().f_code.co_name
        record_execution(self.QMzyme_calls, func)

    # def write_json(self, filename=None):
    #     """
    #     Function to write JSON file containing all information regarding QMzyme run.

    #     :param filename: File name. Defaults to structure id. 
    #     :type filename: str

    #     """
    #     _dict = {}
    #     _dict['Original structure'] = BiopythonWrapper.make_model_dict(self.base)
    #     for model in self.models:
    #         _dict['Model '+str(model.id)] = BiopythonWrapper.make_model_dict(model)
    
    #     QMzyme_dict = {'id': self._id}
    #     QMzyme_dict['Structure file'] = self.pdb_file
    #     QMzyme_dict['QMzyme_calls'] = self.QMzyme_calls
    #     QMzyme_dict = {**QMzyme_dict,**_dict}

    #     if filename == None:
    #         filename = self.id+'.json'
    #     with open(filename, 'w') as f:
    #         json.dump(QMzyme_dict, f, indent=4, sort_keys=True)

    #     func = inspect.currentframe().f_code.co_name
    #     record_execution(self.QMzyme_calls, func)


    ########################
    # DEPRECATED FUNCTIONS #
    ########################

    def QMXTB_input(self, file=None, suffix='', substrate_charge=0, mult=1, qm_atoms='', mem='32GB', nprocs=16, program='orca', qm_input=None, verbose=True):
        warnings.warn("This function is no longer available. Revert back to QMzyme 0.9.34 to use this function.", DeprecationWarning)

    def catalytic_center(self, sel=None, res_name=None, res_number=None, chain=None, output_file=None, save_file=True, save_json=None, verbose=None):
        warnings.warn("This function is no longer available and has been replaced with 'set_catalytic_center'. Revert back to QMzyme 0.9.34 to use the original function.", DeprecationWarning)

    def subsystem(self, distance_cutoff=0, output_file=None, save_file=True, starting_pdb=None, include_residues={}, save_json=None, verbose=None):
        warnings.warn("This function is no longer available and has been replaced with 'within_distance'. Revert back to QMzyme 0.9.34 to use the original function.", DeprecationWarning)

    def size_scan(self, threshold=1000, starting_cutoff=6, output_file=None, verbose=True):
        warnings.warn("This function is no longer available. Revert back to QMzyme 0.9.34 to use this function.", DeprecationWarning)
