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
from datetime import datetime
import os
from aqme.qprep import qprep
from QMzyme.rdkit_wrapper import(
    rdkit_info,
    h_cap,
    atom_coords,
    remove_atoms,
    fix_h_bond,
    centroid_coords,
    define_residue,
    store_mol_pdb,
    mol_from_pdb,
    )
from QMzyme.utils import(
    download,
    get_coords,
    get_atoms,
    get_outlines,
    to_dict,
    name_output,
    coords_from_pdb,
    res_charges,
    )
try:
    from rdkit import Chem
    from rdkit.Chem import rdDistGeom
except ModuleNotFoundError:
    print("RDKit is not installed! You can install the program with 'pip install rdkit' or according to https://anaconda.org/conda-forge/rdkit.")

try:
    import MDAnalysis
    from QMzyme.mdanalysis_wrapper import res_selection, atom_selection
    mdanalysis = True
except:
    mdanalysis = False

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

section_spacer = '#############################################################'

class GenerateModel:

    def __init__(self, calculation='QM-only',
                 protein_file=None,
                 pdb_code=None,
                 save_json=True,
                 verbose=True):
        '''
        Initialize QMzyme model.
        calculation -- string defining the type of calculation this
                    model will be used for. Only option currently is
                    'QM-only'. 'QM/MM' will be supported in the near future.
        protein_file -- string defining the protein or system .pdb
                    file. This can be preprocessed, or a fresh download
                    from rcsb.org.
        pdb_code -- string defining the PDB ID code. This option is
                    required if you do not currently have a protein pdb
                    file for your system of interest. Defining the
                    4-letter PDB code here will automatically download
                    it from the Protein Data Bank. It is highly recommended
                    that you double check the downloaded file and become
                    acquainted with it.
        '''

        if protein_file is None:
            if pdb_code is None:
                raise ValueError("Either protein_file, or pdb_code must"+
                                 " be defined when instantiating model."+
                                 " If you currently do not have a pdb file"+
                                 " use the pdb_code argument to automatically"+
                                 " download the file from the PDB server.")
            else:
                protein_file = download(pdb_code)

        self.protein_file = protein_file
        self.protein_prefix = os.path.basename(protein_file).split('.pdb')[0]
        self.calculation = calculation
        self.verbose = verbose
        self.save_json = save_json

        try:
            x = open(protein_file)
        except FileNotFoundError as e:
            print(e)
        else:
            try:
                self.protein_mol = Chem.MolFromPDBFile(protein_file,
                                                   removeHs=False,
                                                   sanitize=False)
                self.protein_mol.GetAtoms()[0]
            except AttributeError as e:
                print(f"{e}: Protein file {protein_file} was found but"+
                      " rdkit was unable to create mol object."+
                      " This may mean the file does not follow proper PDB"+
                      " formatting.")

        file = self.protein_prefix+'_QMzyme'
        a = 0
        while file+'.json' in os.listdir():
            a += 1
            file = self.protein_prefix + f'_QMzyme_{a+1}'
        self.log_file = file + '.log'
        self.json_file = file + '.json'

        string = "INITIALIZING... QMZYME OBJECT: "
        string += "{}\n".format(self.protein_prefix)
        timestamp = str(datetime.now()).split('.')[0]
        string += "TIMESTAMP: {}\n".format(timestamp)
        if self.verbose == True:
            print(string)

        # initialize dict for json file
        self.dict = {'Starting structure': self.protein_file}
        self.dict['Timestamp'] = str(datetime.now()).split('.')[0]

###############################################################################
    def catalytic_center(self, sel='', res_name=None, res_number=None, chain=None,
                         output_file=None, save_file=True, save_json=None, verbose=None):
        '''
        Function to define the center of the QMzyme model. This is
        typically the ligand/substrate. Currently only supports the
        definition of 1 residue (i.e., lengths of res_name,
        res_number, and chain inputs must not exceed 1. Future
        versions will adapt for this functionality. If you wish to
        use multiple residues right now you can bypass this by
        manipulating the PDB file: for example, change the residue
        names to LIG for all residues you want in the catalytic center,
        and set res_name='LIG'.)
        res_name        - string of the three-letter residue name.
        res_number        - integer of the residue number.
        chain            - string of the chain the residue of interest is in.
        output_file        - string to define the name of the output .pdb file.
                    Generates the following attributes of the class object: a
                    catalytic_center_mol rdkit mol object, catalytic_center_definition,
                    and creates a .pdb file containing onlt the catalytic center atoms.
        '''
        if verbose is None:
            verbose = self.verbose
        if save_json is None:
            save_json = self.save_json
        verbose_str = "INITIALIZING... CATALYTIC CENTER\n"
        timestamp = str(datetime.now())
        verbose_str += "TIMESTAMP: {}\n".format(timestamp)

        count = 0
        definition, catalytic_center=[],[]
        cat_center_def = {}
        file_suffix = '_'

        if mdanalysis is True and sel != ('' or None):
            cat_center_def = res_selection(self.protein_file,sel=sel)
        if mdanalysis is False and sel != ('' or None):
            raise ValueError("The sel argument is only compatible with the MDAnalysis "+ 
                             "package which was not found. Either install MDAnalysis or "+
                             "define the selection using any combination of res_name=str, "+
                             "res_number=int or chain=str that will result in a unique definition.")

        else:
            if chain is not None:
                cat_center_def['Chain'] = chain
                definition.append('Chain')
                catalytic_center.append(chain)
            if res_name is not None:
                cat_center_def['Residue name'] = res_name
                definition.append('res_name')
                catalytic_center.append(res_name)
            if res_number is not None:
                cat_center_def['Residue number'] = res_number
                definition.append('res_number')
                catalytic_center.append(res_number)

        catalytic_center_mol = Chem.RWMol(self.protein_mol)
        previous_res = None
        catalytic_center_residues = []
        for atom in reversed(self.protein_mol.GetAtoms()):
            remove=False
            for key,value in cat_center_def.items():
                if rdkit_info(atom,key) != value:
                    remove=True
            if remove is True:
                catalytic_center_mol.RemoveAtom(atom.GetIdx())
                continue
            count += 1
            current_res = define_residue(atom)
            if current_res != previous_res:
                catalytic_center_residues.append(current_res)
            previous_res=current_res

        if count==0:
            raise ValueError("WARNING: No atoms found matching"+
                    " catalytic center definition.")

        if len(catalytic_center_residues)>1:
            raise Warning("Catalytic center definition is not unique"+
                      " and multiple residues were therefore included: {}."\
                      .format(catalytic_center_residues)+" Please ensure this"+
                      " is what is intended!")

        verbose_str += "DEFINITION: {}\n".format(file_suffix[1:-1])
        n_atoms = catalytic_center_mol.GetNumAtoms()
        verbose_str += "N_ATOMS: {}\n".format(n_atoms)

        if save_file is True:
            if output_file is None:
                output_file = self.protein_prefix+'_catalytic_center.pdb'
            outfile = name_output(name=output_file,suffix='.pdb')
            Chem.MolToPDBFile(catalytic_center_mol,outfile)
            verbose_str += "OUTPUT_FILE: {}\n".format(outfile)
            self.catalytic_center_pdb = get_outlines(outfile)
        if save_file is False:
            self.catalytic_center_pdb = store_mol_pdb(catalytic_center_mol)

        verbose_str += "The following object attributes are now available:\n"
        verbose_str += "\tself.catalytic_center_definition\n"
        verbose_str += "\tself.catalytic_center_mol\n"
        verbose_str += "\tself.catalytic_center_pdb"

        if verbose is True:
            print(verbose_str)

        self.catalytic_center_definition = current_res
        self.catalytic_center_mol = catalytic_center_mol

        if save_json is True:
            info = cat_center_def
            info['Number of atoms'] = n_atoms
            try:
                info['Output file'] = outfile
            except:
                pass
            self.dict = to_dict(key='Catalytic center', data=info, 
                                dict=self.dict, json_file=self.json_file)

###############################################################################
    def subsystem(self, distance_cutoff=0, output_file=None, save_file=True, 
                    starting_pdb=None, include_residues={}, 
                    save_json=None, verbose=None):
        '''
        Function that selects all residues that have at least
        one atom within the cutoff distance from the predefined catalytic
        center atoms.
        distance_cutoff        - integer
        output_file        - string defining the .pdb output file name.
        starting_pdb    - PDB file to be used instead of the file from model 
            initialization. This is useful if, for example, you are creating 
            multiple subsystems from the same initial PDB with varying cutoff distances, in
            which case the starting_pdb should have been generated
            with a larger cutoff distance than what is currently specified.

        Generates subsystem_mol attribute, and creates .pdb file.
        '''
        if verbose is None:
            verbose = self.verbose
        if save_json is None:
            save_json = self.save_json
        verbose_str = "INITIALIZING... SUBSYSTEM SELECTION\n"
        timestamp = str(datetime.now())
        verbose_str += "TIMESTAMP: {}\n".format(timestamp)
        verbose_str += "CUTOFF: {}\n".format(distance_cutoff)
        if distance_cutoff==0:
            raise ValueError("Please specify a distance cutoff by"
                             " 'distance_cutoff={int}'.")

        keep_residue=[]
        if starting_pdb is not None:
            mol = mol_from_pdb(starting_pdb)
        else:
            mol = self.protein_mol
        new_mol = Chem.RWMol(mol)
        
        all_coords = get_coords(self.protein_file)
        lig_coords = coords_from_pdb(self.catalytic_center_pdb)

        for atom in reversed(mol.GetAtoms()):
            current_res = define_residue(atom)
            if current_res in (include_residues or keep_residue):
                continue
            else:
                distances = [np.linalg.norm(all_coords[atom.GetIdx()]-coord) for coord in lig_coords]
                if True in (d < distance_cutoff for d in distances):
                    keep_residue.append(current_res)

        for atom in reversed(mol.GetAtoms()):
            current_res=define_residue(atom)
            if current_res not in (include_residues or keep_residue):
                new_mol.RemoveAtom(atom.GetIdx())

        self.distance_cutoff=distance_cutoff
        self.subsystem_mol = new_mol

        if save_file is True:
            if output_file is None:
                output_file = '{}_'.format(self.protein_prefix)
                output_file += 'subsystem_distance_cutoff'
                output_file += '{}.pdb'.format(self.distance_cutoff)
            outfile = name_output(name=output_file,suffix='.pdb')
            Chem.MolToPDBFile(new_mol,outfile)
            verbose_str += "OUTPUT_FILE: {}\n".format(outfile)
            self.subsystem_pdb = get_outlines(outfile)    
        if save_file is False:
            self.subsystem_pdb = store_mol_pdb(new_mol)

        n_atoms = new_mol.GetNumAtoms()
        verbose_str += "N_ATOMS: {}\n".format(n_atoms)
        verbose_str += "The following object attributes have been generated:\n"
        verbose_str += "\tself.distance_cutoff\n"
        verbose_str += "\tself.subsystem_mol\n"
        verbose_str += "\tself.subsystem_pdb\n"

        if verbose is True:
            print(verbose_str)

        # write json
        if save_json == True:
            info = {
                'Number of atoms': n_atoms,
                'Distance cutoff': self.distance_cutoff,
                }
            try:
                info['Output file'] = outfile
            except:
                info['Output file'] = 'Not saved'
            self.dict = to_dict(key='Subsystem selection', data=info, 
                                dict=self.dict, json_file=self.json_file)

###############################################################################
    def truncate(self, scheme='CA_terminal', output_file=None,
                 skip_res_names=solvent_list, skip_res_numbers=[],
                 remove_res_numbers=[], remove_atom_ids=[],
                 remove_sidechains=[], keep_backbones=[],
                 constrain_atoms=['CA'], add_hydrogens=False,
                 exclude_solvent=False, save_file=True, save_json=None, 
                 verbose=None):
        '''
        Function to prepare truncated QMzyme model.
        scheme            - string to define what truncation scheme to use.
                    Currently only option is 'CA_terminal'. This will remove
                    backbone atoms of residues that are not bondedto other
                    residues. I.e., if you have the following residues in
                    the subsystem (ALA120, GLY121, ASP200), the N-terminus
                    of ALA120 will be removed, the C-terminus of GLY121 will
                    be removed, and both termini will be removed from ASP
                    200 like a typical methyl-capping scheme.
        output_file         - string to define the output .pdb file name.
        skip_res_names        - list containing strings of residue names you
                    do not want to truncate at all.
        skip_res_numbers     - list containing integers of residue numbers you
                    do not want to truncate.
        remove_res_numbers    - list containing integers of residue numbers you
                    want completely removed from the model.
        remove_atom_ids        - list containing integers of atom IDs you want
                    completed removed from the model.
        remove_sidechains    - list containing strings of residue names you
                    want the sidechains removed of. This can be useful if its
                    a bulky residue pointing away from the catalytic center, and
                    only its backbone is important to consider in the model.
        keep_backbones        - list containing strings of residue names that
                    you want to keep backbone atoms of, regardless of the
                    capping scheme.
        constrain_atoms        - list containing strings of the atom names
                    (according to PDB format) that you want to constrain in the
                    calculation.
        add_hydrogens        - boolean, default=False. If True, hydrogens will
                    be added to your model with the reduce package.
                    Recommdended if you started with a .pdb file that did not
                    contain hydrogens, but you should still evaluate the
                    final structure.
        exclude_solvent        - boolean, default=False. If True, solvent
                    molecules will be removed from the model.

        Generates the truncated_subsystem_mol and constrain_atom_list
        attributes, and saves new .pdb file.
        '''
        if verbose is None:
            verbose = self.verbose
        if save_json is None:
            save_json = self.save_json
        verbose_str = "INITIALIZING... SUBSYSTEM TRUNCATION\n"
        timestamp = str(datetime.now())
        verbose_str += "TIMESTAMP: {}\n".format(timestamp)
        verbose_str += "SCHEME: {}\n".format(scheme)
        verbose_str += "CUTOFF: {}\n".format(self.distance_cutoff)

        ### Heidi's to do: add CA capping scheme, create capping summary as a
        # returned item, allow capping summary to be a function input so users
        # can specify the cap they want for each residue.

        if add_hydrogens is True:
            Chem.MolToPDBFile(self.subsystem_mol,'temp0.pdb')
            self.subsystem_mol = utils.add_H(pdb_file='temp0.pdb',remove_file=True)
            os.remove('temp0.pdb')
            self.subsystem_mol = self.protonate_solvent(mol=self.subsystem_mol)

        new_mol = Chem.RWMol(self.subsystem_mol)
        proline_count,bb_atom_count = 0,0
        constrain_list, N_terminus, C_terminus=[],[],[]
        previous_res = None
        remove_ids = []
        residues = []
        fix_N = []

        for atom in reversed(self.subsystem_mol.GetAtoms()):
            current_res = define_residue(atom)
            if current_res==self.catalytic_center_definition:
                continue
            if atom.GetIdx() in remove_atom_ids:
                remove_ids.append(atom.GetIdx())
                continue
            if current_res['Residue name'] in skip_res_names:
                continue
            if current_res['Residue number'] in skip_res_numbers:
                continue
            if current_res['Residue number'] in remove_res_numbers:
                remove_ids.append(atom.GetIdx())
                continue
            if exclude_solvent is True:
                if current_res['Residue name'] in solvent_list:
                    remove_ids.append(atom.GetIdx())

            name = rdkit_info(atom)
            if current_res != previous_res:
                if name == 'C':
                    C_id = atom.GetIdx()
                    C_atom = atom
                    bb_atom_count += 1
                elif name == 'O':
                    O_id = atom.GetIdx()
                    bb_atom_count += 1
                elif name == 'N':
                    N_atom = atom
                    N_id = atom.GetIdx()
                    bb_atom_count +=1
                elif name == 'CA':
                    CA_atom = atom
                    CA_id = atom.GetIdx()
                    bb_atom_count +=1
                if bb_atom_count != 4:
                    continue
                C_bonds = C_atom.GetNeighbors()
                N_bonds = N_atom.GetNeighbors()
                CA_bonds = CA_atom.GetNeighbors()
                previous_res = current_res
                bb_atom_count = 0
            if scheme == 'CA_terminal':
                if 'C' not in rdkit_info(list(N_bonds)):
                    if current_res['Residue name'] == 'PRO':
                        fix_N.append(N_id)
                    else:
                        new_mol = h_cap(new_mol, N_id)
                        for i in N_bonds:
                            if rdkit_info(i) != 'CA':
                                remove_ids.append(i.GetIdx())
                if 'N' not in rdkit_info(list(C_bonds)):
                    remove_ids.append(O_id)
                    new_mol = h_cap(new_mol, C_id)
                    for i in C_bonds:
                        if rdkit_info(i) != 'CA':
                            remove_ids.append(i.GetIdx())
            if scheme == 'CA_all':
                remove_ids.append(C_id)
                remove_ids.append(O_id)
                new_mol = h_cap(new_mol,C_id)
                if current_res['Residue name'] != 'PRO':
                    remove_ids.append(N_id)
                    new_mol = h_cap(new_mol,N_id)
                elif 'C' not in rdkit_info(list(N_bonds)):
                    fix_N.append(N_id)

            # How should free GLY be treated? floating methyl, or keep backbone?
            if current_res['Residue name'] == 'GLY':
                if (C_id and N_id) in remove_ids:
                    remove_ids.append(CA_id)
                    for i in CA_bonds:
                        remove_ids.append(i.GetIdx())

#TO DO: automate code to remove unecessary side chains
# still need to test
        if len(remove_sidechains) > 0:
            print("NEED TO TEST THIS FUNCTION.")
            #for atom in reversed(self.subsystem_mol.GetAtoms()):
            #    current_res = define_residue(atom)
            #    if current_res['Residue name'] in remove_sidechains:
            #        name = rdkit_info(atom)
            #        #if name in ('CA','C','O','N'):
            #        if name in ('CA','C','O','N','H'):
            #            continue
            #        elif name == 'CB':
            #            #cap_atom = new_mol.GetAtomWithIdx(atom.GetIdx())
            #            #cap_atom.SetAtomicNum(1)
            #            #cap_atom.GetPDBResidueInfo().SetName(' H* ')
            #            new_mol = h_cap(new_mol, atom.GetIdx())
            #            continue
            #        #if name == ' H  ':
            #        #    continue
            #        elif 'CA' not in [rdkit_info(x) for x in atom.GetNeighbors()]:
            #            remove_ids.append(atom.GetIdx())

        new_mol = remove_atoms(new_mol,remove_ids)

        #Fix cap atom bond lengths
        new_mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(new_mol,Chem.SanitizeFlags.SANITIZE_FINDRADICALS\
                         |Chem.SanitizeFlags.SANITIZE_KEKULIZE\
                         |Chem.SanitizeFlags.SANITIZE_SETAROMATICITY\
                         |Chem.SanitizeFlags.SANITIZE_SETCONJUGATION\
                         |Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION\
                         |Chem.SanitizeFlags.SANITIZE_SYMMRINGS,\
                         catchErrors=True)
        for atom in new_mol.GetAtoms():
            current_res = define_residue(atom)
            name = rdkit_info(atom)
            if current_res not in residues:
                residues.append(current_res)
            if name == 'H*':
                new_mol = fix_h_bond(new_mol, atom)
            if current_res['Residue name'] == 'PRO':
                if name == 'N':
                    if 'C' not in rdkit_info(list(atom.GetNeighbors())):
                        atom.SetHybridization(Chem.HybridizationType.SP2)
                        fix_N.append(atom.GetIdx())

            # Record atom ids to add to constrain list
            if name in constrain_atoms:
                constrain_list.append(atom.GetIdx()+1)
        if len(fix_N) > 0:
            new_mol = Chem.rdmolops.AddHs(new_mol,addCoords=True,onlyOnAtoms=fix_N)

        if save_file is True:
            if output_file is None:
                output_file = '{}_'.format(self.protein_prefix)
                output_file += 'truncated_subsystem_distance_cutoff'
                output_file += '{}.pdb'.format(self.distance_cutoff)
            outfile = name_output(name=output_file,suffix='.pdb')
            self.filename = outfile
            Chem.MolToPDBFile(new_mol,outfile)
            verbose_str += "OUTPUT_FILE: {}\n".format(outfile)
            #with open(outfile) as f:
                #self.truncated_subsystem_pdb = f.readlines()
            self.truncated_subsystem_pdb = get_outlines(outfile)
        if save_file is False:
            Chem.MolToPDBFile(new_mol,'temp1.pdb')
            self.truncated_subsystem_pdb = get_outlines('temp1.pdb')
            os.remove('temp1.pdb')

        if proline_count > 0:
            print("WARNING: Subsystem model contains {}".\
                  format(proline_count)+" proline(s). The N backbone"
                  " atom was automatically kept.")
        subsystem_charge = 0
        res_list = ''
        for res in residues:
            res_list += res['Residue name']+str(res['Residue number'])+','
            if res['Residue name'] in positive_residues:
                subsystem_charge+=1
            if res['Residue name'] in negative_residues:
                subsystem_charge-=1

        verbose_str += "N_ATOMS: {}\n".format(new_mol.GetNumAtoms())
        verbose_str += "CHARGE: {}\n".format(subsystem_charge)
        verbose_str += "NOTE: charge does NOT include the catalytic center "
        verbose_str += "and is based on AMBER amino acid naming conventions.\n "
        verbose_str += "MODEL_COMPONENTS: {}\n".format(res_list[:-1])
        verbose_str += "The following object attributes are now available:\n"
        verbose_str += "\tself.subsystem_charge\n"
        verbose_str += "\tself.model_atom_count\n"
        verbose_str += "\tself.truncated_subsystem_mol\n"
        verbose_str += "\tself.truncated_subsystem_pdb\n"
        verbose_str += "\tself.constrain_atom_list\n"
        verbose_str += "\tself.residues\n"
        if verbose is True:
            print(verbose_str)

        self.subsystem_charge = subsystem_charge
        self.model_atom_count = new_mol.GetNumAtoms()
        self.truncated_subsystem_mol = new_mol
        self.constrain_atom_list = constrain_list
        self.model_residues = residues
        self.residues = residues
        
        # write json
        if save_json == True:
            info = {
                'Number of atoms': new_mol.GetNumAtoms(),
                'Distance cutoff': self.distance_cutoff,
                'Residues': residues,
                'C-alpha atom indices': constrain_list,
                'Subsystem charge': subsystem_charge
                }
            try:
                info['Output file'] = outfile
            except:
                pass
            self.dict = to_dict(key='Truncated subsystem', data=info, 
                                dict=self.dict, json_file=self.json_file)

###############################################################################
    def size_scan(self, threshold=1000, starting_cutoff=6,
                     output_file=None, verbose=True):
        cutoff = starting_cutoff-1
        pass_verbose = verbose
        verbose_str = "PERFORMING... SIZE SCAN\n"
        verbose_str += "Scanning distance cutoffs starting from "
        verbose_str += "{} Angstroms until minimal ".format(starting_cutoff)
        verbose_str += "model size of {} atoms is met.".format(threshold)
        if verbose is True:
            print(verbose_str)
        n_atoms = 0
        while n_atoms<threshold:
            cutoff+=1
            self.subsystem(distance_cutoff=cutoff,
                             save_file=False,verbose=pass_verbose)
            self.truncate(save_file=False,verbose=pass_verbose)
            n_atoms = self.model_atom_count
        if output_file is None:
            output_file = '{}_'.format(self.protein_prefix)
            output_file += 'truncated_subsystem_distance_cutoff'
            output_file += '{}_'.format(self.distance_cutoff)
            output_file += '{}atoms.pdb'.format(self.model_atom_count)
        self.truncate(save_file=True,
                      output_file=output_file,
                      verbose=pass_verbose)

###############################################################################
    def QMXTB_input(self,file=None,suffix='',substrate_charge=0,mult=1,
                    qm_atoms='',mem='32GB',nprocs=16,program='orca',
                    qm_input=None,verbose=True):
        '''
        CURRENTLY UNDER DEVELOPMENT. 
        qm_atoms is a string that gets passed to MDAnalysis to make the selection.'''
        chrg = self.subsystem_charge+substrate_charge
        if '!' in qm_input:
            qm_input = qm_input.split('!')[-1]
        if 'qm/xtb' not in qm_input.lower():
            qm_input = 'QM/XTB '+qm_input
        if file is None:
            try:
                file = self.filename
            except:
                file = self.protein_prefix+'truncated_subsystem_distance_'
                file += 'cutoff'+self.distance_cutoff+'.pdb'
                with open(file, "a") as f:
                    f.writelines(self.truncated_subsystem_pdb)
        # calculate charge of qm region
        qm_chrg = res_charges(res_selection(pdb=file, sel=qm_atoms))

            
###############################################################################
    def QM_input(self,file=None,suffix='', non_protein_charge=0,mult=1,mem='32GB',
                 nprocs=16,program='gaussian',level_of_theory=None,verbose=True, 
                 save_json=False):
        if suffix.startswith('_') is True:
            suffix = suffix.split('_')[-1]
        if file is None:
            try:
                file = self.filename
            except:
                file = self.protein_prefix+'truncated_subsystem_distance_'
                file += 'cutoff'+self.distance_cutoff+'.pdb'
                with open(file, "a") as f:
                    f.writelines(self.truncated_subsystem_pdb)
        qm_input = {"starting structure": file,
                    "atom types": get_atoms(file),
                    "coords": get_coords(file),
                    "charge": self.subsystem_charge+non_protein_charge,
                    "multiplicity": mult,
                    "memory": mem,
                    "# processors": nprocs,
                    "level of theory": level_of_theory}
        if save_json is True:
            outfile = file.split('.')[0]+suffix+'.json'
            with open(outfile, 'w+') as f:
                json.dump(qm_input, f, indent=1)
        try: qm_input["frozen atoms"] = self.constrain_atom_list
        except: pass
        verbose_str, json_info = qm_only(qm_input, program, suffix=suffix, verbose=verbose)
        if verbose is True:
            print(verbose_str)
        # write json
        info = {'Calculation file': file.split('.pdb')[0]+suffix+'.com',
                'Calculation': level_of_theory,
                'Charge': self.subsystem_charge+non_protein_charge,
                'Multiplicity': mult}
        self.to_dict(section='QM preparation', dict=info)


