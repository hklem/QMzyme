import os
import json
from datetime import datetime
import numpy as np
from urllib.request import urlopen
from rdkit import Chem
#from QMzyme.rdkit_wrapper import check_pdb_rdkit
from QMzyme.protein_parser import collect_pdb_data, pdb_deposit_info


##### GENERAL #####

delimeter='---------------------------------------------------------\n'
solvent_list=['HOH','WAT','T3P','SOL'] 
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

protein_residues = ['ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'CYM', 'CYS', 'CYX',
                    'GLH', 'GLN', 'GLU', 'GLY', 'HIS', 'HID', 'HIE', 'HIP',
                    'HYP', 'ILE', 'LEU', 'LYN', 'LYS', 'MET', 'PHE', 'PRO',
                    'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HSE', 'HSD', 'HSP',
                    'SEC', 'PYL']

positive_residues = ['HIP', 'LYS', 'ARG']

negative_residues = ['ASP', 'GLU']
###############################################################################
def name_output(name=None, suffix='.out'):
    if name.endswith(suffix):
        return name
    else:
        return name+suffix

###############################################################################
def write_log(file, text):
    log = open(file,'w')
    log.write("{}\n".format(text))
    log.write(delimeter)
    log.close()

def log(string, file):
    with open(self.log_file, 'a') as f:
        string += section_spacer
        string += '\n'
        f.writelines(string)

def json_type_encoder(dict):
    for key in dict.keys():
        if isinstance(dict[key], np.integer):
            dict[key] = int(dict[key])
        if isinstance(dict[key], np.floating):
            dict[key] = float(dict[key])
        if isinstance(dict[key], np.ndarray):
            dict[key] = dict[key].tolist()
    return dict

def to_dict(key=None, data=None, dict={}, json_file=''):
    if type(data) is dict:
        data = json_type_encoder(data)
    if key == 'Catalytic center':
        if key in dict.keys():
            raise Exception("A catalytic center has been previous " +
                            "defined in this object. Please initialize a " +
                            "new object via QMzyme.GenerateModel() to " +
                            "continue with a new catalytic center definiton."
                            )
        dict[key] = data
    elif 'QMzyme 1' not in dict.keys():
        dict['QMzyme 1'] = {key: data}
    else:
        number = 0
        for model in dict.keys():
            if 'QMzyme' in model:
                number += 1
        if f'QMzyme {number}' not in dict.keys():
            dict[f'QMzyme {number}'] = {key: data}
        elif key not in dict[f'QMzyme {number}'].keys():
            dict[f'QMzyme {number}'][key] = data
        else:
            dict[f'QMzyme {number+1}'] = {key: data}
    with open(json_file, "w") as f:
        json.dump(dict, f, indent=4)

    return dict


###############################################################################
def download(pdb_code):
    base_url = 'http://www.pdb.org/pdb/download/downloadFile.do?fileFormat'+\
     '=pdb&compression=NO&structureId='
    for structure in pdb_code.split():
         pdb_url = base_url + structure[:4]
         output_file = structure[:4] + '.pdb'

         with urlopen(pdb_url) as response, open(output_file, 'wb') as outfile:
             data = response.read()
             outfile.write(data)
             print("Downloading {} as {}.".format(structure[:4], output_file))
             print(pdb_url)
                
    return output_file

###############################################################################
def check_pdb(file,clean=False):
    '''
    Function to assess PDB format, fix any issues that might break 
    the QMzyme code, and gather useful basic information.
    clean            - boolean, default=True. If False, the PDB format will 
    not be checked. Usually not a problem when downloading 
    straight from rcsb.org.

    Prints out basic information of the pdb file, and  if 
    issues are found the PDB file will be fixed, and the original 
    will be copied to a new file with the suffix '_original.pdb'.
    '''

    wat_count, protein_res_count, non_protein_res_count = 0,0,0
    data, protein_seq, non_protein_seq = [], [], []
    h_present, no_chain_info = False, False
    previous_res = None
    edits_made = False
    residue_count = 0
    res_num_previous = None
    non_protein_residues = {}
    non_protein_residues['Chain'] = []
    non_protein_residues['Name'] = []
    non_protein_residues['Number'] = []
    non_protein_chemical_name = {}

    pdb_data = collect_pdb_data(file=file)
    deposit_info = pdb_deposit_info(file=file)

    with open(file,'r') as f:
        data = f.readlines()
    for line in data:
        if 'HET   ' in line:
            non_protein_res_count += 1
            chain, name, number = line[12], line[7:10], line[13:17]
            non_protein_residues['Chain'].append(chain)
            non_protein_residues['Name'].append(name)
            non_protein_residues['Number'].append(number)
        if 'HETNAM' in line:
            non_protein_chemical_name[line[11:14]]=line[15:].split('  ')[0]

    residues_reordered=False
    if clean==True:
        for i,line in enumerate(data):
            if pdb_info('record_type',line) in ['ATOM','HETATM']:
                atom_type = pdb_info('atom_name',line)
                res_num = pdb_info('res_number',line)
                res_name = pdb_info('res_name',line)
                if res_num != res_num_previous:
                    if res_name not in solvent_list:
                        residue_count += 1
                        if int(res_num.split()[0]) != residue_count:
                            residues_reordered = True
                            edits_made=True
                            if residue_count>999:
                                res_num = str(residue_count)
                            elif residue_count>99:
                                res_num = ' '+str(residue_count)
                            elif residue_count>9:
                                res_num = '  '+str(residue_count)
                            else:
                                res_num = '   '+str(residue_count)
                            #data[i]=line[0:22]+res_num+line[26:]
                            #line=line[0:22]+res_num+line[26:]
                            if pdb_info('res_name', data[i+1]) == res_name:
                                residue_count -= 1
                res_num_previous = res_num
                try:
                    element = pdb_info('element_symbol', line)
                except ValueError:
                    element=None

                if atom_type[0]!=' ':
                    if atom_type[0]=='H':
                        if len(atom_type)==4:
                            if element is None:
                                edits_made=True
                                #data[i]=line[0:76]+' H'+line[78:]
                            continue
                if atom_type[:2] in elements:
                    if element is None:
                        edits_made=True
                        #data[i]=line[0:76]+atom_type[:2]+line[78:]
                    elif atom_type[:2] == element:
                        continue
                if atom_type[1] in elements:
                    if element is None:
                        edits_made=True
                        #data[i]=line[0:76]+atom_type[:2]+line[78:]
                    elif atom_type[:2] == element:
                        continue
                    #else:
                else:
                    if atom_type[1]+atom_type[2].lower() in elements:
                        edits_made=True
                        #data[i]=line[0:12]+atom_type[1:]+' '+line[16:76]+\
                            #atom_type[1:3]+line[78:]
                            
                        print("SUSPECTED ISSUE FOUND ON ATOM{}: Element symbol"\
                              .format(pdb_info('atom_number',line))+
                              " with two letters ({}) must begin in the 13th"\
                              .format(atom_type[1:3])+" columnspace according"+
                              " to PDB format.")
                    elif atom_type[2] in elements:
                        edits_made=True
                        unk_sym = atom_type[:2].split()[0]
                        #if len(unk_sym) == 2:
                        #    data[i] = line[0:12]+' '+atom_type[2:]+' '+\
                        #              line[16:76]+' '+atom_type[2]+line[78:]
                        #if len(unk_sym) == 1:
                        #    data[i] = line[0:12]+' '+atom_type[2:]+' '+\
                        #              line[16:76]+' '+atom_type[2]+line[78:]
                        
                        print("SUSPECTED ISSUE FOUND ON ATOM{}"\
                              .format(pdb_info('atom_number',line)) +
                              ": Atom type ({}) ".format(atom_type) +
                              "is preceeded by unknown symbols ({})."\
                              .format(atom_type[:2]))
        #if residues_reordered==True:
        #    print("WARNING: Some residues were renumbered.")
        #if edits_made is True:
        #    new_file = self.protein_prefix+'_original.pdb'
        #    cmd = 'cp {} {}'.format(self.protein_file,new_file)
        #    os.system(cmd)
        #    print("Saved original pdb file as {}.".format(new_file)+
        #          " Suspected issue(s)/warning(s) have been resolved in {}."\
        #          .format(self.protein_file))
        #    with open(self.protein_file, 'w+') as g:
        #        for i,line in enumerate(data):
        #            g.writelines(data[i])
    protein_res_count,h_present,non_protein_residues,non_protein_residue_count,non_protein_chemical_names = check_pdb_rdkit(file)
    return protein_res_count,h_present,non_protein_residues,non_protein_residue_count,non_protein_chemical_names



###############################################################################
def add_H(pdb_file=None, output_file=None, remove_file=False):
    '''
    Function that calls the reduce function from AmberTools 
    to add hydrogens.
    pdb_file        - string defining the .pdb file to be reduced.
    output_file        - string defining the reduced .pdb file name.

    Returns new rdkit mol object that has been reduced and saves 
    the reduced .pdb. 
    '''
        
    if pdb_file is None:
        raise ValueError("PDB file must be defined.")
    if output_file is None:
        output_file=pdb_file.split('.pdb')[0]+'_reduced.pdb'
    cmd = "reduce -NOFLIP -Quiet {} > {}".format(pdb_file,output_file)
    os.system(cmd)
    new_mol = Chem.MolFromPDBFile(output_file,
                                  removeHs=False,
                                  sanitize=False)
    if remove_file is True:
        os.remove(output_file)
        
    return new_mol

##### PDB SPECIFIC #####
pdb_format = {'record_type':[0,6],
              'atom_number':[6,11], 
              'atom_name':[12,16],
              'alt_loc':[16],
              'res_name':[17,20],
              'chain_id':[21],
              'res_number':[22,26],
              'insertion_code':[26],
              'x':[30,38],
              'y':[38,46],
              'z':[46,54],
              'occupancy':[54,60],
              'temperature_factor':[60,66],
              'seg_id':[72,76],
              'element_symbol':[76,78],
              'charge':[78,80]}

def pdb_info(info,line):
    return line[pdb_format[info][0]:pdb_format[info][1]]

def collect_pdb_data(file=None,data=None):
    if file is not None:
        try: 
            with open(file, 'r') as f:
                data=f.readlines()
        except:
            raise FileNotFoundError("test file {} not found.".format(file))
    else:
        if data is None:
            raise Exception("Must define either file or lines from a pdb file.")

    pdb_data = {'res_name_sequence':[],'atom_sequence':[],'number_of_atoms':0,'atom_coords':[],'res_num_sequence':[]}
    for line in data:
        if pdb_info('record_type',line) in ['ATOM','HETATM']:
            pdb_data['res_name_sequence'].append(pdb_info('res_name',line))
            pdb_data['res_num_sequence'].append(pdb_info('res_number',line))
            pdb_data['number_of_atoms'] += 1
            pdb_data['atom_sequence'].append(pdb_info('atom_name',line))
            pdb_data['atom_coords'].append((pdb_info('x',line),
                                          pdb_info('y',line),
                                          pdb_info('z',line)))

    return pdb_data

def rmsd(xyz1, xyz2):
    delta = xyz1 - xyz2
    rmsd = (delta ** 2.0).sum(1).mean() ** 0.5
    return rmsd

def compute_translation_and_rotation(mobile, target):
    #meta data
    n_atoms = mobile.shape[0]
    n_dim = mobile.shape[1]
    mu1 = np.zeros(n_dim)
    mu2 = np.zeros(n_dim)
    for i in range(n_atoms):
        for j in range(n_dim):
            mu1[j] += mobile[i,j]
            mu2[j] += target[i,j]
    mu1 /= n_atoms
    mu2 /= n_atoms
    mobile = mobile - mu1
    target = target - mu2

    correlation_matrix = np.dot(np.transpose(mobile), target)
    V, S, W_tr = np.linalg.svd(correlation_matrix)
    if np.linalg.det(V) * np.linalg.det(W_tr) < 0.0:
        V[:, -1] = -V[:, -1]
    rotation = np.dot(V, W_tr)

    translation = mu2 - np.dot(mu1,rotation)

    return translation, rotation

def kabsch_transform(mobile, translation, rotation):
    mobile_prime = np.dot(mobile,rotation) + translation
    return mobile_prime

def coords_from_json(file):
    coords = []
    data = get_outlines(file)
    for atom in range(len(np.array(data["atoms"]["coords"]["3d"]))//3):
        coords.append([data["atoms"]["coords"]["3d"][atom*3+x] for x in range(3)])
    return np.array(coords)

def coords_from_pdb(file):
    coords = []
    if type(file) is str:
        data = get_outlines(file)
    else:
        data = file
    count = 0
    for i, line in enumerate(data):
        if line.split()[0] in ['ATOM', 'HETATM']:
            coords.append([float(line.split()[x+6]) for x in range(3)])
    return np.array(coords)

def coords_from_gout(file):
    coords = []
    data = get_outlines(file)
    done = False
    for i, line in enumerate(data):
        if 'Standard orientation' in line:
            for j in range(i+5, i+20000):
                if '--------' in data[j]:
                    break
                else:
                    coords.append([float(data[j].split()[x+3]) for x in range(3)])
            break
    return coords

def get_frozen_coords(coords, frozen_atoms):
    f_coords = np.array([coords[x] for x in frozen_atoms])
    return f_coords
        
def get_coords(file):
    dispatch = {
                'pdb': lambda file: coords_from_pdb(file),
                'log': lambda file: coords_from_gout(file),
                'out': lambda file: coords_from_gout(file),
                'json': lambda file: coords_from_json(file),
                }
    coords = dispatch[file.split('.')[-1]](file)
    return coords

def frozen_atoms_from_json(file):
    data = get_outlines(file)
    return data["metadata"]["frozenatoms"] # zero indexed

def frozen_atoms_from_gout(file):
    frozen_atoms = []
    outlines = get_outlines(file)
    done = False
    atom_idx = -1 # to make zero indexed
    for i,line in enumerate(outlines):
        if done is True:
            break
        if 'Symbolic Z-matrix' in line:
            done = True
            for j in range(i + 2, i + 20000):
                atom_idx += 1
                l = outlines[j].split()
                if len(l) == 0:
                    break
                elif l[1] == '-1':
                    frozen_atoms.append(atom_idx)
    return frozen_atoms

def frozen_atoms_from_ginp(file):
    frozen_atoms = []
    outlines = get_outlines(file)
    atom_idx = -1 # to make zero indexed
    for line in outlines:
        l = line.split()
        try:
            if l[0] in elements:
                atom_idx += 1
                if l[1] == '-1':
                    frozen_atoms.append(atom_idx)
        except:
            pass
    return frozen_atoms

def get_frozen_atoms(file):
    dispatch = {
                'log': lambda file: frozen_atoms_from_gout(file),
                'out': lambda file: frozen_atoms_from_gout(file),
                'com': lambda file: frozen_atoms_from_ginp(file),
                'gjf': lambda file: frozen_atoms_from_ginp(file),
                'json': lambda file: frozen_atoms_from_json(file),
                }
    coords = dispatch[file.split('.')[-1]](file)
    return coords

def get_outlines(file):
    with open(file, 'r') as f:
        if file.endswith('json'):
            outlines = json.load(f)
        else:
            outlines = f.readlines()
    return outlines

def align(mobile, reference, atom_indices=[]):
    '''
    Mobile and reference are numpy arrays with dimensions (N, 3) where N is the number of atoms.
    Atom by default is an empty list and alignment will be performed for all atoms.
    '''
    if type(mobile) is str:
        mobile = get_coords(mobile)
    if type(reference) is str:
        reference = get_coords(reference)
        
    original = mobile
    if len(atom_indices) > 0:
        ref = []
        mob = []
        for atom in atom_indices:
            #ref.append([float(reference[atom+pad].split()[6+x]) for x in range(3)])
            ref.append(reference[atom])
            mob.append(mobile[atom])
        reference = np.array(ref)
        mobile = np.array(mob)
    
    t, r = compute_translation_and_rotation(mobile, reference)
    aligned_coords = kabsch_transform(original, t, r)
    return aligned_coords

def write_pdb_coords(original_pdb, new_coords, suffix='optimized', dec=4):
    file_name = os.path.basename(original_pdb).split('.')[0]
    if type(new_coords) is str:
        new_coords = get_coords(new_coords)
    pdb = get_outlines(original_pdb)
    atom = -1
    with open(file_name+'_'+suffix+'.pdb', 'w+') as f:
        for i,line in enumerate(pdb):
            if line.startswith('ATOM'):
                atom += 1
                for j,coord in enumerate(['x', 'y', 'z']):
                    xyz = pdb_info(coord, line)
                    new_coord = format(new_coords[atom][j], '.{}f'.format(dec))
                    if len(new_coord) < len(xyz):
                        for s in range(len(xyz)-len(new_coord)):
                            new_coord = ' '+new_coord
                    line = line.replace(xyz, new_coord)
                f.write(line)
            elif line.startswith('CONECT'):
                continue
            else:
                f.write(line)

def get_atoms(file):
    atoms = []
    with open(file, 'r') as f:
        for line in f.readlines():
            if line.startswith('ATOM'):
                atoms.append(pdb_info('element_symbol',line).strip())
    return atoms

def res_charges(residues):
        charge=0 
        if type(residues) is dict:
            for i,res in enumerate(residues):
                residues[i] = res['Residue Name']
        for res in residues:
            if res in positive_residues:
                charge+=1
            elif res in negative_residues:
                charge-=1

        return charge