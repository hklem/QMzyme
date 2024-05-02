###############################################################################
# Code written by Heidi Klem while at
# Colorado State University as a graduate student
# in the Paton and McCullagh groups and at the
# National Institute of Standards and Technology
# as an NRC Postdoc (Fed).
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
Module in charge of creating the QM enzyme models given a starting structure. 
This module can be used to call CalculateModel to generate QM calculation 
input files.

Required input
...............

    *   Prepared starting structure in PDB format.
"""

from QMzyme import MDAnalysisWrapper
from QMzyme.QMzymeBuilder import QMzymeModel
from QMzyme.QMzymeBuilder import QMzymeRegion
import inspect
import os

protein_residues = ['ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'CYM', 'CYS', 'CYX',
                    'GLH', 'GLN', 'GLU', 'GLY', 'HIS', 'HID', 'HIE', 'HIP',
                    'HYP', 'ILE', 'LEU', 'LYN', 'LYS', 'MET', 'PHE', 'PRO',
                    'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HSE', 'HSD', 'HSP',
                    'SEC', 'PYL']

class GenerateModel():

    def __init__(self, pdb_file, id=None):
        """
        Constructor method for GenerateModel class.

        :param pdb_file: PDB file name.  
        :type pdb_file: str, required
        :param id: Name to associate structure with, defaults to file name. Will be used in naming of any files produced.
        :type id: str
        """
        u = MDAnalysisWrapper.init_universe(pdb_file)
        if id is None:
            id = os.path.basename(pdb_file).split('.')[0]

        setattr(self, 'pdb_file', pdb_file)
        setattr(self, 'universe', u)
        setattr(self, 'id', id)
        # Raise warning if hydrogens are not present- can cause issues with truncate methods
        Hs_present = False
        for atom in u.atoms:
            if atom.element == 'H':
                Hs_present = True
                break
        if Hs_present == False:
            raise UserWarning("PDB structure does not contain hydrogens. Please pre-process this structure. We recommend using tLeap. See Documentation.")

    def __repr__(self):
        return f"<QMzyme id={self.id}>"
    
    # def set_region(self, region_name, AtomGroup):
    #     QMzymeModel.set_region(self, region_name, AtomGroup)

    def set_catalytic_center(
        self,
        selection):        
        '''
        Function to define the center of the QMzyme model. This is
        typically the ligand/substrate. To add onto an already defined catalytic center 
        run set_catalytic_center again with the new selection and set 'overwrite=False'. 

        :param selection: Selection of atoms to be made- based on `MDAnalysis selection command language <https://docs.mdanalysis.org/stable/documentation_pages/selections.html>`_.
        :type selection: str, required`

        :param overwrite: To clear any existing catalytic_center definition. Set to False if you would like to append current catalytic_center. Defaults to True.
        :type overwrite: bool, default: overwrite=True
        '''
        catalytic_center = self.universe.select_atoms(selection) 
        QMzymeModel.set_region(self, 'catalytic_center', catalytic_center)


    def within_distance(self, distance_cutoff, store_region=False, region_name=None):
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
        if not hasattr(self, 'catalytic_center'):
            raise UserWarning("You must first define a catalytic_center. See method `set_catalytic_center()`.")

        neighbors = MDAnalysisWrapper.get_neighbors(
            self.universe.select_atoms('all'),
            self.catalytic_center.get_AtomGroup(),
            #self.catalytic_center._AtomGroup, 
            distance_cutoff)

        residues = neighbors.residues.sorted_unique
        atoms = residues.atoms
        
        if store_region is True:
            if region_name is None:
                region_name = f'cutoff_{distance_cutoff}'
            QMzymeModel.set_region(self, region_name, atoms)
            region = QMzymeModel.get_region(self, region_name)
            setattr(region, 'neighbors', QMzymeRegion(f'neighbors_{region_name}', neighbors).atoms)
        
        if store_region is False:
            if region_name is None:
                region_name = f'cutoff_{distance_cutoff}'
            return (QMzymeRegion(f'neighbors_{region_name}', neighbors), QMzymeRegion(f'neighbors_byres_{region_name}', atoms))
    
    def distance_scan(self, cutoffs=None, minimum_size=None, store_region=False):
        regions = []
        if minimum_size is not None:
            cutoff = 0
            n_atoms = 0
            while n_atoms < minimum_size:
                cutoff +=1
                regions.append(self.within_distance(distance_cutoff=cutoff, store_region=store_region))
                n_atoms = regions[-1].n_atoms
            if store_region == False:
                return regions[-1]
        else:
            for cutoff in cutoffs:
                regions.append(self.within_distance(distance_cutoff=cutoff, store_region=store_region))
            if store_region == False:
                return regions
            

    def truncate(self, region, scheme='CA_terminal', 
                 store_region=False, region_name=None):
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

        ######## can i use self._neighbors to make a private attribute of a QMzymeRegion???
        removed_Nterm = []
        removed_Cterm = []
        truncated_model = []
        new_atoms = []

        # Need to init new universe so original one isn't overwritten
        u = MDAnalysisWrapper.init_universe(self.pdb_file)
        residues = [MDAnalysisWrapper.get_parallel_residue(res, u) for res in region.get_AtomGroup().residues]
        print(residues)
        for i, res in enumerate(residues):
            if res.resname not in protein_residues:
                truncated_model += list(res.atoms)
                continue
            side_chain_atoms = []
            keep_sc = False
            #Loop over side chain atoms. If none are within cutoff distance don't include them
            for atom in res.atoms:
                if atom.name not in ['C', 'O', 'N', 'H']:
                    truncated_model.append(atom)
            C_atom = MDAnalysisWrapper.get_atom(res, 'C')
            O_atom = MDAnalysisWrapper.get_atom(res, 'O')
            try:
                next_res = residues[i+1]
                next_N_atom = MDAnalysisWrapper.get_atom(next_res, 'N')
            except:
                # last residue, so cap CA at C_atom
                C_cap = MDAnalysisWrapper.cap_backbone_CA(C_atom)
                new_atoms.append(C_cap)
                continue
            # If very first residue, cap CA at N_atom... unless it's proline
            if res == residues[0]:
                removed_Nterm.append(res)
                N_atom = MDAnalysisWrapper.get_atom(res, 'N')
                if res.resname == 'PRO':
                    cap_atom = MDAnalysisWrapper.cap_backbone_N(N_atom)
                    new_atoms.append(cap_atom)
                else:
                    N_cap= MDAnalysisWrapper.cap_backbone_CA(N_atom)
                    new_atoms.append(N_cap)
            # Look at next residue in list to decide how to treat current C term and next N term
            #if next res is not in sequence, truncate current res CA at C_atom and
            #truncate next res CA at N_atom... unless it's proline
            if res.resid+1 != next_res.resid:
                removed_Cterm.append(res)
                next_N_atom = MDAnalysisWrapper.get_atom(next_res, 'N')
                if next_res.resname == 'PRO':
                    cap_atom = MDAnalysisWrapper.cap_backbone_N(next_N_atom)
                    new_atoms.append(cap_atom)
                else:
                    removed_Nterm.append(next_res)
                    C_cap = MDAnalysisWrapper.cap_backbone_CA(C_atom)
                    new_atoms.append(C_cap)
                    N_cap = MDAnalysisWrapper.cap_backbone_CA(next_N_atom)
                    new_atoms.append(N_cap)
            #If next res is in sequence keep all C term and next res N term backbone atoms
            else:
                truncated_model.append(C_atom)
                truncated_model.append(O_atom)
                truncated_model.append(next_N_atom)
                if next_res.resname != 'PRO':
                    next_H_atom = MDAnalysisWrapper.get_atom(next_res, 'H')
                    truncated_model.append(next_H_atom)

        combined = truncated_model+new_atoms
        u = MDAnalysisWrapper.build_universe(combined)

        if store_region is False:
            if region_name is None:
                region_name = region.name+'_truncated'
            return QMzymeRegion(region_name, u.select_atoms('all'))
        if store_region is True:
            if region_name is None:
                region_name = region.name+'_truncated'
            QMzymeModel.set_region(self, region_name, u.select_atoms('all'))
    
        
    





    
