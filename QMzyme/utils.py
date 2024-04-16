###############################################################################
# Code written by Heidi Klem while at
# Colorado State University as a graduate student
# in the Paton and McCullagh groups and at the
# National Institute of Standards and Technology
# as an NRC Postdoc (Fed).
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

'''Utility functions for QMzyme package.'''

import os
import json
import numpy as np
import datetime
from QMzyme.Biopython.Data.PDBData import protein_letters_3to1_extended

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

protein_residues = [r for r in protein_letters_3to1_extended.keys()]

positive_residues = ['HIP', 'LYS', 'ARG']

negative_residues = ['ASP', 'GLU']


def record_execution(QMzyme_details, function_name):
    now = datetime.datetime.now()
    date = now.strftime("%m/%d/%Y")
    time = now.strftime("%H:%M:%S")
    count = 1
    for key in QMzyme_details.keys():
        if function_name in key:
            count += 1
    if count > 1:
        if function_name in list(QMzyme_details.keys()):
            QMzyme_details[f'{function_name}_run_1'] = QMzyme_details[function_name]
            del QMzyme_details[function_name]
        function_name += f'_run_{count}'

    QMzyme_details[function_name] = {}
    QMzyme_details[function_name]['Date'] = date
    QMzyme_details[function_name]['Time'] = time


def filename_format(name=None, suffix='.out'):
    if name.endswith(suffix):
        return name
    else:
        if suffix.startswith('.'):
            return f"{name}{suffix}"
        else: return f"{name}.{suffix}"

def json_type_encoder(dict):
    for key in dict.keys():
        if isinstance(dict[key], np.integer):
            dict[key] = int(dict[key])
        if isinstance(dict[key], np.floating):
            dict[key] = float(dict[key])
        if isinstance(dict[key], np.ndarray):
            dict[key] = dict[key].tolist()
    return dict

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

def rmsd(xyz1, xyz2, align=False):
    if align == True:
        t, r = compute_translation_and_rotation(xyz1, xyz2)
        xyz1 = kabsch_transform(xyz1, t, r)
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
        if pdb_info('record_type',line).strip() in ['ATOM','HETATM']:
            coords.append([float(pdb_info(coord, line)) for coord in ['x', 'y', 'z']])
            #coords.append([float(line.split()[x+6]) for x in range(3)])
    return np.array(coords)

def coords_from_gout(file):
    coords = []
    data = get_outlines(file)
    for i, line in enumerate(data):
        if 'Standard orientation' in line:
            for j in range(i+5, i+20000):
                if '--------' in data[j]:
                    break
                else:
                    coords.append([float(data[j].split()[x+3]) for x in range(3)])
            break
    return np.array(coords)

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
    Atom indices are 0 indexed.
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
            if pdb_info('record_type',line).strip() in ['ATOM','HETATM']:
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


class AddArgs:
    pass

#def set_args(kwargs, var_dict):
def set_args(**kwargs):
    # set default options and options provided
    options = AddArgs()
    # dictionary containing default values for options
    #for key in var_dict:
    #    vars(options)[key] = var_dict[key]
    for key in kwargs:
        vars(options)[key] = kwargs[key]
        # if key in var_dict:
        #     vars(options)[key] = kwargs[key]
        # elif key.lower() in var_dict:
        #     vars(options)[key.lower()] = kwargs[key.lower()]
        # else:
        #     print("Warning! Option: [", key,":",kwargs[key],"] provided but no option exists, try the online documentation to see available options for each module.",)
    return vars(options)


    ########################
    # DEPRECATED FUNCTIONS #
    ########################

    def add_H(pdb_file=None, output_file=None, remove_file=False):
        warnings.warn("This function is no longer available. Revert back to QMzyme 0.9.34 to use this function.", DeprecationWarning)

    def check_pdb(file, clean=False):
        warnings.warn("This function is no longer available. Revert back to QMzyme 0.9.34 to use this function.", DeprecationWarning)

    def download(pdb_code):
        warnings.warn("This function is no longer available. Revert back to QMzyme 0.9.34 to use the original function.", DeprecationWarning)

    def to_dict(key=None, data=None, dict={}, json_file=''):
        warnings.warn("This function is no longer available. Revert back to QMzyme 0.9.34 to use this function.", DeprecationWarning)

    def collect_pdb_data(file=None, data=None):
        warnings.warn("This function is no longer available. Revert back to QMzyme 0.9.34 to use this function.", DeprecationWarning)

    def log(string, file):
        warnings.warn("This function is no longer available. Revert back to QMzyme 0.9.34 to use this function.", DeprecationWarning)

    def log(string, file):
        warnings.warn("This function is no longer available. Revert back to QMzyme 0.9.34 to use this function.", DeprecationWarning)

    def write_log(file, text):
        warnings.warn("This function is no longer available. Revert back to QMzyme 0.9.34 to use this function.", DeprecationWarning)

