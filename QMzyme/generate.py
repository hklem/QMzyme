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
from datetime import datetime
from rdkit import Chem
import os
from rdkit.Chem import rdDistGeom
from aqme.qprep import qprep
from QMzyme import utils
from QMzyme.mdanalysis_wrapper import res_selection
from QMzyme.rdkit_wrapper import(
    rdkit_info,
    h_cap,
    atom_coords,
    remove_atoms,
    fix_h_bond,
    centroid_coords,
    define_residue,
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

section_spacer = '#############################################################'

class GenerateModel:

    def __init__(self, calculation='QM-only',
                 protein_file=None,
                 pdb_code=None):
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
                protein_file = utils.download(pdb_code)

        self.protein_file=protein_file
        self.protein_prefix = protein_file.split('.pdb')[0]
        self.calculation = calculation

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
                      " formatting. Run self.check_pdb() function to clean.")
        # begin storing information
        file = self.protein_prefix+'_QMzyme_1'
        while file+'.log' in os.listdir():
            a = int(file.split('QMzyme_')[1].split('.')[0])
            file = self.protein_prefix
            file += '_QMzyme_{}'.format(a+1)
        self.log_file = file+'.log'
        self.json_file = file+'.json'
        string = "INITIALIZING... QMZYME OBJECT: "
        string += "{}\n".format(self.protein_prefix)
        timestamp = str(datetime.now())
        string += "TIMESTAMP: {}\n".format(timestamp)
        self.log(string)

        # initialize dict for json file
        self.to_dict(dict={'Starting structure': self.protein_file})

        available_pks = []
        try:
            import MDAnalysis
            self.mdanalysis = True
        except:
            self.mdanalysis = False

###############################################################################
    def name_output(self, name=None, suffix='.out'):
        if name.endswith(suffix):
            return name
        else:
            return name+suffix

###############################################################################
    def store_mol_pdb(self, mol):
        Chem.MolToPDBFile(mol,'temp2.pdb')
        with open('temp2.pdb') as f:
            data = f.readlines()
        os.remove('temp2.pdb')
        return data

###############################################################################
    def log(self, string):
        with open(self.log_file, 'a') as f:
            string += section_spacer
            string += '\n'
            f.writelines(string)

###############################################################################
    def to_dict(self, section=None, dict={}):
        import json
        if section is None:
            self.data = dict
        elif section in self.data.keys():
            for key in self.data[section].keys():
                if type(dict[key]) is not list:
                    if type(self.data[section][key]) is not list:
                        self.data[section][key] = [self.data[section][key]]
                else:
                    if type(self.data[section][key][0]) is not list:
                        self.data[section][key] = [self.data[section][key]]
                self.data[section][key].append(dict[key])
        else:
            self.data[section] = dict
        with open(self.json_file, "w") as f:
            json.dump(self.data, f)

###############################################################################
    def catalytic_center(self, sel='', res_name=None, res_number=None, chain=None,
                         output_file=None, save_file=True, verbose=True):
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
        verbose_str = "INITIALIZING... CATALYTIC CENTER\n"
        timestamp = str(datetime.now())
        verbose_str += "TIMESTAMP: {}\n".format(timestamp)

        count = 0
        definition, catalytic_center=[],[]
        cat_center_def = {}
        file_suffix = '_'

        if self.mdanalysis is True and sel != '':
            cat_center_def = res_selection(self.protein_file,sel=sel)

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
            outfile = self.name_output(name=output_file,suffix='.pdb')
            Chem.MolToPDBFile(catalytic_center_mol,outfile)
            verbose_str += "OUTPUT_FILE: {}\n".format(outfile)
            with open(outfile) as f:
                self.catalytic_center_pdb = f.readlines()
        if save_file is False:
            self.catalytic_center_pdb = self.store_mol_pdb(catalytic_center_mol)
        self.log(verbose_str)

        verbose_str += "The following object attributes are now available:\n"
        verbose_str += "\tself.catalytic_center_definition\n"
        verbose_str += "\tself.catalytic_center_mol\n"
        verbose_str += "\tself.catalytic_center_pdb"

        if verbose is True:
            print(verbose_str)

        self.catalytic_center_definition = current_res
        self.catalytic_center_mol = catalytic_center_mol

        #write json
        print(cat_center_def)
        #for key in cat_center_def.keys():
        #    if len(cat_center_def[key]) == 1:
        #        if type(cat_center_def[key]) is list:
        #            cat_center_def[key] = cat_center_def[key][0]
        info = cat_center_def
        info['Number of atoms'] = n_atoms
        try:
            info['Output file'] = outfile
        except:
            pass
        self.to_dict(section='Catalytic center', dict=info)

###############################################################################
    def active_site(self, distance_cutoff=0, include_residues=[],
                    output_file=None, save_file=True,
                    solvent=solvent_list, intermediate_mol=None,
                    center=None, verbose=True):
        '''
        Function that selects all residues that have at least
        one atom within the cutoff distance from the predefined catalytic
        center atoms. If there are residues you want to include that may
        not be within the distance cutoff you can specify that in
        include_residues=[].
        distance_cutoff        - integer
        include_residues    - list of a tuple of strings:
                    (str(chain),str(res_name),int(res_number)).
                    I.e., include_residues=[('A','GLY',101),('A','ASP',20)].
                    It follows the define_residue() format.
        output_file        - string defining the .pdb output file name.
        intermediate_mol    - rdkit mol object to be used instead of the
            protein mol that was created during model initialization. This is
            useful if, for example, you are creating multiple active sites
            from the same initial PDB with varying cutoff distances, in
            which case the intermediate_mol object must have been generated
            with a larger cutoff distance than what is currently specified.
        center        - (Work in progress) string defining the .pdb file
            of your catalytic center. Default is None. This is helpful to
            use if you want to alter the ligand in any way, for example,
            add hydrogens to certain atoms. Please not the coordinates of
            the molecule need to be consistent with the coordinates of
            the protein.

        Generates active_site_mol attribute, and creates .pdb file.
        '''
        verbose_str = "INITIALIZING... ACTIVE SITE SELECTION\n"
        timestamp = str(datetime.now())
        verbose_str += "TIMESTAMP: {}\n".format(timestamp)
        verbose_str += "CUTOFF: {}\n".format(distance_cutoff)
        if distance_cutoff==0:
            raise ValueError("Please specify a distance cutoff by"
                             " 'distance_cutoff={int}'.")

        res_name, res_number, res_chain = [], [], []
        add_residue=[]
        if intermediate_mol is not None:
            mol = intermediate_mol
        else:
            mol = self.protein_mol

        temp_mol = Chem.RWMol(mol)
        new_mol = Chem.RWMol(mol)
        if center is not None:
            self.catalytic_center_mol = Chem.MolFromPDBFile(center)

        centroid = centroid_coords(self.catalytic_center_mol)
        distances = [np.linalg.norm(atom_coords(mol, atom)-centroid)\
            for atom in self.catalytic_center_mol.GetAtoms()]
        distance_buffer = np.max(distances)

        # First pass: goes over all atoms in protein and cat. centroid
        for atom in reversed(mol.GetAtoms()):
            if rdkit_info(atom).startswith('H'):
                continue
            current_res=define_residue(atom)
            if current_res in (include_residues or add_residue):
                continue
            else:
                coords = atom_coords(mol,atom)
                atomic_distance = np.linalg.norm(coords-centroid)
                if atomic_distance < distance_cutoff+distance_buffer:
                    add_residue.append(current_res)
                else:
                    temp_mol.RemoveAtom(atom.GetIdx())

        # Second pass: goes over all atoms in catalytic center
        keep_residue=[]
        res_dict = {'Chain':[],
                    'Residue Name':[],
                    'Residue Number':[]}

        for atom1 in reversed(temp_mol.GetAtoms()):
            current_res=define_residue(atom1)
            if current_res in include_residues:
                keep_residue.append(current_res)
                continue
            if current_res in keep_residue:
                continue
            else:
                coords1 = atom_coords(temp_mol,atom1)
            for atom2 in self.catalytic_center_mol.GetAtoms():
                if current_res in keep_residue:
                    continue
                coords2 = atom_coords(self.catalytic_center_mol,atom2)
                atomic_distance = np.linalg.norm(coords1-coords2)
                if atomic_distance < distance_cutoff:
                    res_dict['Chain'].append(current_res['Chain'])
                    res_dict['Residue Name'].append(current_res['Residue name'])
                    res_dict['Residue Number'].append(current_res['Residue number'])
                    keep_residue.append(current_res)

        for atom in reversed(mol.GetAtoms()):
            current_res=define_residue(atom)
            if current_res not in keep_residue:
                new_mol.RemoveAtom(atom.GetIdx())

        if save_file is True:
            if output_file is None:
                output_file = '{}_'.format(self.protein_prefix)
                output_file += 'active_site_distance_cutoff'
                output_file += '{}.pdb'.format(self.distance_cutoff)
            outfile = self.name_output(name=output_file,suffix='.pdb')
            Chem.MolToPDBFile(new_mol,outfile)
            verbose_str += "OUTPUT_FILE: {}\n".format(outfile)
            with open(outfile) as f:
                self.active_site_pdb = f.readlines()
        if save_file is False:
            self.active_site_pdb = self.store_mol_pdb(new_mol)

        n_atoms = new_mol.GetNumAtoms()
        verbose_str += "N_ATOMS: {}\n".format(n_atoms)
        self.log(verbose_str)

        verbose_str += "The following object attributes have been generated:\n"
        verbose_str += "\tself.distance_cutoff\n"
        verbose_str += "\tself.active_site_mol\n"
        verbose_str += "\tself.active_site_pdb\n"
        verbose_str += "\tself.active_site_residues\n"

        if verbose is True:
            print(verbose_str)
        self.distance_cutoff=distance_cutoff
        self.active_site_mol = new_mol
        self.active_site_residues = res_dict

        # write json
        info = {
            'Number of atoms': n_atoms,
            'Distance cutoff': self.distance_cutoff,
            }
        try:
            info['Output file'] = outfile
        except:
            pass
        self.to_dict(section='Active site selection',dict=info)

###############################################################################
    def protonate_solvent(self,mol=None,solvent=solvent_list):
        '''
        Function to add hydrogens to water molecules because reduce does not
        recognize all water residue names. I found that for some weird
        reason you cannot do rdkit's Chem.rdmolops.AddHs and include
        PDB residue information when Hs are already present, so we can't
        simply add Hs to the relevant water oxygen atoms because then
        the PDB information cannot be included if any other Hs are present
        in the structure. To avoid that situation, this function separates
        water molecules from the mol object, protonates them, adds coordinate
        information back, and then combines with final mol object. Note: if
        any solvent molecules already have bound H, they will be skipped.

        mol : rdkit mol object, must be defined.
        solvent : List of strings containing names of solvent residues in PDB.
            default is solvent=['HOH','WAT','T3P','SOL']

        Returns mol object.
        '''
        solvent_mol = Chem.RWMol()
        no_solvent_mol = Chem.RWMol(mol)
        pos = []
        for atom in reversed(mol.GetAtoms()):
            if rdkit_info(atom,'Residue name') in solvent:
                if rdkit_info(atom) == 'O':
                    if len(atom.GetNeighbors()) != 2:
                        pos.append(atom_coords(mol,atom))
                        no_solvent_mol.RemoveAtom(atom.GetIdx())
                        solvent_mol.AddAtom(atom)
        if solvent_mol.GetNumAtoms() > 0:
            rdDistGeom.EmbedMolecule(solvent_mol)
            for i,atom in enumerate(reversed(solvent_mol.GetAtoms())):
                solvent_mol.GetConformer().SetAtomPosition(atom.GetIdx(),pos[i])
                atom.SetHybridization(Chem.HybridizationType.SP3)

            solvent_mol = Chem.rdmolops.AddHs(solvent_mol,
                                              addCoords=True,
                                              addResidueInfo=True)

            combined_mol = Chem.CombineMols(no_solvent_mol,solvent_mol)
            return combined_mol
        else:
            return mol


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
        self.log(verbose_str)
        n_atoms = 0
        while n_atoms<threshold:
            cutoff+=1
            self.active_site(distance_cutoff=cutoff,
                             save_file=False,verbose=pass_verbose)
            self.truncate(save_file=False,verbose=pass_verbose)
            n_atoms = self.model_atom_count
        if output_file is None:
            output_file = '{}_'.format(self.protein_prefix)
            output_file += 'truncated_active_site_distance_cutoff'
            output_file += '{}_'.format(self.distance_cutoff)
            output_file += '{}atoms.pdb'.format(self.model_atom_count)
        self.truncate(save_file=True,
                      output_file=output_file,
                      verbose=pass_verbose)

###############################################################################
    def truncate(self, scheme='CA_terminal', output_file=None,
                 skip_res_names=solvent_list, skip_res_numbers=[],
                 remove_res_numbers=[], remove_atom_ids=[],
                 remove_sidechains=[], keep_backbones=[],
                 constrain_atoms=['CA'], add_hydrogens=False,
                 exclude_solvent=False, save_file=True, verbose=True):
        '''
        Function to prepare truncated QMzyme model.
        scheme            - string to define what truncation scheme to use.
                    Currently only option is 'CA_terminal'. This will remove
                    backbone atoms of residues that are not bondedto other
                    residues. I.e., if you have the following residues in
                    the active site (ALA120, GLY121, ASP200), the N-terminus
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
                    a bulky residue pointing away from the active site, and
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

        Generates the truncated_active_site_mol and constrain_atom_list
        attributes, and saves new .pdb file.
        '''
        verbose_str = "INITIALIZING... ACTIVE SITE TRUNCATION\n"
        timestamp = str(datetime.now())
        verbose_str += "TIMESTAMP: {}\n".format(timestamp)
        verbose_str += "SCHEME: {}\n".format(scheme)
        verbose_str += "CUTOFF: {}\n".format(self.distance_cutoff)

        ### Heidi's to do: add CA capping scheme, create capping summary as a
        # returned item, allow capping summary to be a function input so users
        # can specify the cap they want for each residue.

        if add_hydrogens is True:
            Chem.MolToPDBFile(self.active_site_mol,'temp0.pdb')
            self.active_site_mol = utils.add_H(pdb_file='temp0.pdb',remove_file=True)
            os.remove('temp0.pdb')
            self.active_site_mol = self.protonate_solvent(mol=self.active_site_mol)

        new_mol = Chem.RWMol(self.active_site_mol)
        proline_count,bb_atom_count = 0,0
        constrain_list, N_terminus, C_terminus=[],[],[]
        previous_res = None
        remove_ids = []
        residues = []
        fix_N = []

        for atom in reversed(self.active_site_mol.GetAtoms()):
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
            #for atom in reversed(self.active_site_mol.GetAtoms()):
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
                output_file += 'truncated_active_site_distance_cutoff'
                output_file += '{}.pdb'.format(self.distance_cutoff)
            outfile = self.name_output(name=output_file,suffix='.pdb')
            self.filename = outfile
            Chem.MolToPDBFile(new_mol,outfile)
            verbose_str += "OUTPUT_FILE: {}\n".format(outfile)
            with open(outfile) as f:
                self.truncated_active_site_pdb = f.readlines()
        if save_file is False:
            Chem.MolToPDBFile(new_mol,'temp1.pdb')
            with open('temp1.pdb') as f:
                self.truncated_active_site_pdb = f.readlines()
            os.remove('temp1.pdb')

        if proline_count > 0:
            print("WARNING: Active site model contains {}".\
                  format(proline_count)+" proline(s). The N backbone"
                  " atom was automatically kept.")
        active_site_charge = 0
        res_list = ''
        for res in residues:
            res_list += res['Residue name']+str(res['Residue number'])+','
            if res['Residue name'] in positive_residues:
                active_site_charge+=1
            if res['Residue name'] in negative_residues:
                active_site_charge-=1

        verbose_str += "N_ATOMS: {}\n".format(new_mol.GetNumAtoms())
        verbose_str += "CHARGE: {}\n".format(active_site_charge)
        verbose_str += "NOTE: charge does NOT include the catalytic center "
        verbose_str += "and is based on AMBER amino acid naming conventions.\n "
        verbose_str += "MODEL_COMPONENTS: {}\n".format(res_list[:-1])
        self.log(verbose_str)
        verbose_str += "The following object attributes are now available:\n"
        verbose_str += "\tself.active_site_charge\n"
        verbose_str += "\tself.model_atom_count\n"
        verbose_str += "\tself.truncated_active_site_mol\n"
        verbose_str += "\tself.truncated_active_site_pdb\n"
        verbose_str += "\tself.constrain_atom_list\n"
        verbose_str += "\tself.residues\n"
        if verbose is True:
            print(verbose_str)

        self.active_site_charge = active_site_charge
        self.model_atom_count = new_mol.GetNumAtoms()
        self.truncated_active_site_mol = new_mol
        self.constrain_atom_list = constrain_list
        self.model_residues = residues
        self.residues = residues

        # write json
        info = {
            'Number of atoms': new_mol.GetNumAtoms(),
            'Distance cutoff': self.distance_cutoff,
            'Residues': residues,
            'C-alpha atom indices': constrain_list,
            'Active site charge': active_site_charge
            }
        try:
            info['Output file'] = outfile
        except:
            pass
        self.to_dict(section='Truncated active site', dict=info)

###############################################################################
    def QM_input(self,file=None,suffix='',substrate_charge=0,mult=1,mem='32GB',nprocs=16,qm_input=None,verbose=True):
        chrg = self.active_site_charge+substrate_charge
        if file is None:
            try:
                file = self.filename
            except:
                file = self.protein_prefix+'truncated_active_site_distance_'
                file += 'cutoff'+self.distance_cutoff+'.pdb'
                with open(file, "a") as f:
                    f.writelines(self.truncated_active_site_pdb)
        qprep(files=file,
              charge=chrg,
              mult=mult,
              freeze=self.constrain_atom_list,
              qm_input=qm_input,
              program='gaussian',
              mem=mem,
              nprocs=nprocs,
              suffix=suffix)
        if suffix!='':
            suffix = '_'+suffix+'_'
        os.rename('./QCALC/'+file.split('.pdb')[0]+suffix+'_conf_1.com',
                  './QCALC/'+file.split('.pdb')[0]+suffix+'.com')
        verbose_str = "INITIALIZING... AQME.QPREP QM INPUT FILE GENERATION\n"
        verbose_str += "STARTING_STRUCTURE: {}\n".format(file)
        verbose_str += "CALCULATION: {}\n".format(qm_input)
        verbose_str += "CHARGE: {}\n".format(chrg)
        verbose_str += "MULTIPLICITY: {}\n".format(mult)
        verbose_str += "MEM: {}\n".format(mem)
        verbose_str += "NPROCS: {}\n".format(nprocs)
        verbose_str += "FROZEN_ATOMS: {}\n".format(self.constrain_atom_list)
        verbose_str += "INPUT_FILE: {}\n".format(file.split('.pdb')[0]+suffix+'.com')
        if verbose is True:
            print(verbose_str)
        self.log(verbose_str)

        # write json
        info = {
            'Calculation file': file.split('.pdb')[0]+suffix+'.com',
            'Calculation': qm_input,
            'Charge': chrg,
            'Multiplicity': mult,
            'Memory': mem,
            'Number of processors': nprocs,
            }
        self.to_dict(section='QM preparation', dict=info)


