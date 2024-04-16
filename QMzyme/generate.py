###############################################################################
# Code written by Heidi Klem while at
# Colorado State University as a graduate student
# in the Paton and McCullagh groups and at the
# National Institute of Standards and Technology
# as an NRC Postdoc (Fed).
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

'''Generate QM-based enzyme models.'''

import numpy as np
import json
import inspect
import warnings
from QMzyme.calculate import *
from QMzyme.Biopython.Structure import Structure
from QMzyme import BiopythonWrapper
from QMzyme.BiopythonWrapper import res_dispatch
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

class GenerateModel(Structure):

    def __init__(self, structure, id=None, model_id=0):
        """
        Constructor method for GenerateModel class.

        :param structure: PDB or mmCIF file.  
        :type structure: str, required
        :param id: Name to associate structure with, defaults to file name.
        :type id: str
        :param model_id: Model ID to base original structure on, defaults to 0. 
        :type model_id: int
        """

        if type(structure) == str:
            filename = structure
            structure = BiopythonWrapper.load_structure(structure, id)
            setattr(structure, 'pdb_file', filename)
        self.__dict__ = structure.__dict__.copy()

        Hs_present = False
        for atom in self.get_atoms():
            if atom.element == 'H':
                Hs_present = True
        if Hs_present == False:
            raise Exception("PDB structure does not contain hydrogens. Please pre-process this structure.")

        if not hasattr(self, 'base'):
            #self.base = self.child_dict[model_id]
            self.base = self.child_list[model_id]
        self.models = []
        self.QMzyme_calls = {}
        func = inspect.currentframe().f_code.co_name
        record_execution(self.QMzyme_calls, func)
        setattr(self, 'id', id)
        delattr(self, 'xtra')


    def __repr__(self):
        return f"<QMzyme Structure id={self._id}>"


    def set_catalytic_center(
        self, 
        overwrite=True,
        **kwargs
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
        overwrite : bool, default=True
            To clear any existing catalytic_center definition. Set to 
            False if you would like to append current catalytic_center.

        Notes
        -----
        '''
        _dict=(kwargs)
        residues = []
        for res in self.base.get_residues(): 
            if False not in [val == res_dispatch[key](res) for key, val in _dict.items()]:
                residues.append(res)   
        if overwrite is True:
            self.catalytic_center = []
            self.catalytic_center = residues
        if overwrite is False:
            self.catalytic_center += residues

        func = inspect.currentframe().f_code.co_name
        record_execution(self.QMzyme_calls, func)


    def within_distance(self, distance_cutoff, store_model=True):
        '''
        Function to select all residues that have at least one atom within a
        specified distance to any atom in self.catalytic_center. By default,
        the resulting model is appended to self.models. 

        Parameters
        ----------
        distance_cutoff : int or float, required
            Numerical value specifying selection cutoff.

        Notes
        -----
        The cutoff is not done radially, unless your catalytic_center definition 
        only contains one atom. Therefore, the shape of the selection region depends
        on the shape of the catalytic_center.
        '''

        residues = []
        center_coords = []
        for res in self.catalytic_center:
            for atom in res.get_atoms():
                center_coords.append(atom.get_coord())
        neighbors = BiopythonWrapper.get_neighbors(self.base, center_coords, distance_cutoff)
        residues = [atom.get_parent() for atom in neighbors]
        residues = list(set(residues))
        residues = BiopythonWrapper.order_residues(residues)
        method = {'type': inspect.currentframe().f_code.co_name,
                  'cutoff': distance_cutoff, 
                  'catalytic_center': self.catalytic_center}
        if store_model is True:
            self.store_model(residues, method)
        else:
            return BiopythonWrapper.init_model(self, residues, method)
        
        func = inspect.currentframe().f_code.co_name
        record_execution(self.QMzyme_calls, func)


    def truncate(self, scheme='CA_terminal'):
        '''
        Function to remove extraneous atoms in preparation for model calculation.
        This will be performed only on the most recently created model.

        Parameters
        ----------
        scheme : str, default='CA_terminal'
            See documentation for explanations of each truncation scheme.

        Notes
        -----
        Currently, the only available scheme is CA_terminal, which will remove 
        backbone atoms only on residues without their sequence neighbor present in
        the model, and cap the C-alpha atom with hydrogen(s). I.e., if the following
        resnumbers are in the model: [15, 23, 24, 25], then the N-terminal and C-terminal 
        backbone atoms will be removed from residue 15 and the C-alpha will be converted 
        to a methyl group, the N-terminal backbone atoms of residue 23 will be removed, but
        the C-terminal backbone atoms will remain, as will the N-terminal backbone atoms of
        residues 24 and 25, but the C-terminal backbone atoms of 25 will be removed. 
        Additional schemes will be created in the future.

        The added Hydrogens will have a bond length of 1.00, along the bond vector of the
        original atom the H is replacing.
        '''
        for res in self.child_list[-1].get_residues():
            if res.resname not in protein_residues:
                continue
            if BiopythonWrapper.has_Nterm_neighbor(res) is False and res.resname != 'PRO':
                BiopythonWrapper.cap_terminus(res, 'N')
            if BiopythonWrapper.has_Cterm_neighbor(res) is False:
                BiopythonWrapper.cap_terminus(res, 'C')
        self.models[-1] = self.child_list[-1]

        func = inspect.currentframe().f_code.co_name
        record_execution(self.QMzyme_calls, func)


    def store_model(self, residues, method):
        """
        Function to store new model to QMzyme structure object.

        Parameters
        ----------
        residues : list
            List of residue objects that comprise the model.
        method : dict
            Dictionary containing details of how that model was generated.

        Notes
        -----
        """
        m = BiopythonWrapper.init_model(self, residues, method)
        self.add(m)
        if self.child_list[0] == self.base:
            del self.child_list[0]
            del self.child_dict[0]
        self.models.append(self.child_list[-1])
        print(f"Model {self.child_list[-1].id} created by the {method['type']} method has been stored to {self.__repr__()}.")

        func = inspect.currentframe().f_code.co_name
        record_execution(self.QMzyme_calls, func)

    def write_pdb(self, entity=None, filename=""):
        """
        Function to write PDB file.

        Parameters
        ----------
        entity : Biopython structure, model, or chain object, defaults to last generated model
        filename : str
            File name (should have '.pdb' suffix). Defaults to object id. 

        Notes
        -----
        """
        if entity == None:
            entity = self.models[-1]
        entity.pdb_file = BiopythonWrapper.write_pdb(entity, filename)

        func = inspect.currentframe().f_code.co_name
        record_execution(self.QMzyme_calls, func)


    def calculateQM(self, model=None, functional=None, 
                basis_set=None, opt=True, freq=True, freeze_atoms = [], 
                charge = 0, mult = 1, mem='32GB', 
                nprocs=16, program='gaussian', suffix=''):
        """
        Function to generate QM only calculation input file.

        Parameters
        ----------
        model : QMzyme model object, defaults to last model in the list.
        functional : str, QM theory functional to use.
        basis_set : str, QM theory basis set to use. 
        opt : bool, set up input to perform optimization if true. Default is True.
        freq : bool, set up input to perform frequency analysis if true. Default is True.
        freeze_atoms : list containing the PDB style names of any atoms to be frozen during optimization.
            E.x., ['CA'] will assign frozen atom flags to all C-alpha atoms.
        charge : int, electronic charge of model. 
        mult : int, multiplicity of model.
        mem : str, amount of memory to specify in input file. Default is '32GB'.
        nprocs: int, amount of processors to specify for job. Default is 16. 
        program: str, name of QM program to be used. Options are 'gaussian' or 'orca'. 
        suffix: str, to be appended to end of calculation input file name. Default is ''.

        Notes
        -----
        """
        if model == None:
            model = self.child_list[-1]
        if charge != None:
            model.charge = charge
        if mult != None:
            model.multiplicity = mult
        CalculateQM(model, functional, basis_set, opt, freq, freeze_atoms, charge, mult, mem, nprocs, program, suffix)

        func = inspect.currentframe().f_code.co_name
        record_execution(self.QMzyme_calls, func)

    def write_json(self, filename=None):
        """
        Function to write JSON file containing all information regarding QMzyme run.

        Parameters
        ----------
        filename : str
            File name. Defaults to structure id. 

        Notes
        -----
        """
        _dict = {}
        _dict['Original structure'] = BiopythonWrapper.make_model_dict(self.base)
        for model in self.models:
            _dict['Model '+str(model.id)] = BiopythonWrapper.make_model_dict(model)
    
        QMzyme_dict = {'id': self._id}
        QMzyme_dict['Structure file'] = self.pdb_file
        QMzyme_dict['QMzyme_calls'] = self.QMzyme_calls
        QMzyme_dict = {**QMzyme_dict,**_dict}

        if filename == None:
            filename = self.id+'.json'
        with open(filename, 'w') as f:
            json.dump(QMzyme_dict, f, indent=4, sort_keys=True)

        func = inspect.currentframe().f_code.co_name
        record_execution(self.QMzyme_calls, func)


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
