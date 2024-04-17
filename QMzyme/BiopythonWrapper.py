import os
import numpy as np
import  QMzyme
from QMzyme.Biopython.StructureBuilder import StructureBuilder
from QMzyme.Biopython import NeighborSearch
from QMzyme.Biopython.PDBParser import PDBParser
from QMzyme.Biopython.MMCIFParser import MMCIFParser
from QMzyme.Biopython.PDBIO import PDBIO
from QMzyme.utils import filename_format

model_dispatch = {'id': lambda model: model.id,
                  'full_id': lambda model: model.full_id,
                  'serial_number': lambda model: model.serial_number,
                  'chains_list': lambda model: model.child_list,
                  'chains_dict': lambda model: model.child_dict,
                  'xtra': lambda model: model.extra,
                  'structure': lambda model: model.parent}

chain_dispatch = {'id': lambda chain: chain.id,
                  'full_id': lambda chain: chain.full_id,
                  'residues_list': lambda chain: chain.child_list,
                  'residues_dict': lambda chain: chain.child_dict,
                  'xtra': lambda chain: chain.extra,
                  'model': lambda chain: chain.parent,
                  'structure': lambda chain: model_dispatch['structure'](chain.parent)}

res_dispatch = {'resname': lambda res: res.resname,
                'resnumber': lambda res: res.id[1],
                'chain': lambda res: res.parent.id,
                'segid': lambda res: res.segid,
                'id': lambda res: res.id,
                'full_id': lambda res: res.full_id,
                'atoms_list': lambda res: res.child_list,
                'atoms_dict': lambda res: res.child_dict,
                'xtra': lambda res: res.extra,
                'model': lambda res: chain_dispatch['model'](res.parent),
                'structure': lambda res: chain_dispatch['structure'](res.parent)}

atom_dispatch = {'name': lambda atom: atom.name,
                 'fullname': lambda atom: atom.fullname,
                 'coord': lambda atom: atom.coord,
                 'b_factor': lambda atom: atom.bfactor,
                 'bfactor': lambda atom: atom.bfactor,
                 'occupancy': lambda atom: atom.occupancy,
                 'altloc': lambda atom: atom.altloc,
                 'full_id': lambda atom: atom.full_id,
                 'id': lambda atom: atom.id,
                 'serial_number': lambda atom: atom.serial_number,
                 'element': lambda atom: atom.element,
                 'mass': lambda atom: atom.mass,
                 'pqr_charge': lambda atom: atom.pqr_charge,
                 'radius': lambda atom: atom.radius,
                 'xtra': lambda atom: atom.extra,
                 'residue': lambda atom: atom.parent,
                 'resname': lambda atom: res_dispatch['resname'](atom.parent),
                 'resnumber': lambda atom: res_dispatch['resnumber'](atom.parent),
                 'chain': lambda atom: res_dispatch['chain'](atom.parent),
                 'segid': lambda atom: res_dispatch['segid'](atom.parent),
                 'model': lambda atom: res_dispatch['model'](atom.parent),
                 'structure': lambda atom: res_dispatch['structure'](atom.parent)}


def write_pdb(pdb_object, filename="", verbose=True):
    io = PDBIO()
    io.set_structure(pdb_object)
    if filename == "":
        for i in pdb_object.full_id:
            filename += f"{i}_"
        filename = filename[:-1]

    filename = filename_format(filename,'.pdb')
    io.save(f"{filename}", preserve_atom_numbering=True)
    if verbose == True:
        print(f"File {filename} created.")
    return filename

def order_residues(residues_list):
    new_list = []
    chains = list(set([res.parent.id for res in residues_list]))
    for chain in chains:
        res_list = [res for res in residues_list if res.parent.id == chain]
        max = np.max([res.id[1] for res in res_list])
        min = np.min([res.id[1] for res in res_list])
        for i in range(min, max+1):
            for res in res_list:
                if res.id[1] == i:
                    new_list.append(res)
    return new_list

def count_atoms(entity):
    return len(list(entity.get_atoms()))

def count_residues(entity):
    return len(list(entity.get_residues()))

def count_nonprotein_residues(entity):
    count = 0
    for res in entity.get_residues():
        if res.resname not in QMzyme.protein_residues:
            count += 1
    return count

def count_chains(entity):
    return len(list(entity.get_chains()))

def load_structure(file, id):
    if file.endswith('pdb'):
        p = PDBParser(QUIET=True)
    elif file.endswith('cif'):
        p = MMCIFParser(QUIET=True)
    else:
        print(f"File format for {file} not recognized. Please use filename ending in '.pdb' or '.cif'.")
    if id is None:
        id = os.path.basename(file).split('.')[0]
    structure = p.get_structure(f'{id}', f'{file}')
    return structure

def get_neighbors(entity, coords: [list], cutoff):
    #return will be one list of all atoms that were within cutoff
    neighbors = []
    neighbor_search = NeighborSearch(list(entity.get_atoms()))
    if len(np.shape(coords)) == 1:
        coords = [coords]
    for coord in coords:
        neighbors += neighbor_search.search(coord,cutoff)
    return list(set(neighbors)) 

def init_model(structure, model_residues, method = {}):
    """
    xtra is a dictional where you can add whatever you want.
    """ 
    s = StructureBuilder()
    s.init_structure(structure.id)
    #s.structure = structure # with this option, things get messed up because of duplicates
    id = structure.child_list[-1].id + 1
    s.init_model(model_id = id)
    #structure.init_model(model_id = id)
    #for c in structure.parent.get_chains():
    for c in structure.get_chains():
        for res in model_residues:
            if res.parent.id != c.id:
                continue
            s.init_chain(chain_id = c.id)        
            s.init_seg(segid = ' ')
            resname = res.resname
            field = res.id[0]
            resseq = res.id[1]
            icode = res.id[2]
            s.init_residue(resname, field, resseq, icode)
            for atom in res.get_atoms():
                name = atom.name
                coord = atom.coord
                b_factor = atom.bfactor
                occupancy = atom.occupancy
                altloc = atom.altloc
                fullname = atom.fullname
                serial_number = atom.serial_number
                element = atom.element
                pqr_charge = atom.pqr_charge
                radius = atom.radius
                is_pqr = False
                if pqr_charge is not None:
                    is_pqr = True
                s.init_atom(name, coord, b_factor, occupancy, altloc, fullname, serial_number, element, pqr_charge, radius, is_pqr)
    s.model.method = method

    return s.model

# def check_terminal_neighbors(model):
#     for chain in model.get_chains():
#         prev_res = chain.child_list[0]
#         prev_res.truncation = {'N_terminus': 'cut', 'C_terminus' : 'cut'}
#         for res in chain.child_list[1:]:
#             if res.resname not in QMzyme.protein_residues:
#                 continue
#             res.truncation = {'N_terminus': 'cut', 'C_terminus' : 'cut'}
#             if res.id[1]-1 == prev_res.id[1]:
#                 res.truncation['N_terminus'] = 'keep'
#                 prev_res.truncation['C_terminus'] = 'keep'
#             elif res.resname == 'PRO':
#                 res.truncation['N_terminus'] = 'keep'     
#             prev_res = res

def check_terminal_neighbor(residue, inc):
    resid = residue.id[1]
    for res in residue.parent.get_residues():
        if res.id[1] == resid+inc:
            return True
    return False
        
def has_Nterm_neighbor(residue):
    return check_terminal_neighbor(residue,-1)

def has_Cterm_neighbor(residue):
    return check_terminal_neighbor(residue,1)

def h_cap(atom, bonded_atom, Hbond_length = 1.00):
    coords = change_bond_length(bonded_atom.coord, atom.coord, 
                                                    Hbond_length)
    if atom.id == 'N':
        id = 'HN'
    if atom.id == 'C':
        id = 'HC'
    else:
        id = 'H'
    new_atom = {'element': 'H', 'name': 'Hcap', 
                'coord': coords, 'fullname': 'Hcap', 'mass': 1.00794,
                'id': id}
    alter_atom(atom, new_atom)

def change_bond_length(fixed_coords, mobile_coords, new_length):
    M = new_length/np.linalg.norm(fixed_coords-mobile_coords)
    q = fixed_coords-(M*(fixed_coords-mobile_coords))
    return q

def alter_atom(atom, new_atom_dict):
    #setattr(atom, 'original_atom', (atom.parent.full_id, (atom.id, atom.altloc)))
    for key, val in new_atom_dict.items():
        if hasattr(atom, key):
            setattr(atom, key, val)
    setattr(atom, 'full_id', (atom.parent.full_id, (atom.id, atom.altloc)))
    
def remove_atom(atom):
    atom.parent.detach_child(atom.id)

def cap_terminus(residue, terminus):
    """
    terminus options are 'N' or 'C'.
    """
    model = residue.parent.parent
    remove_list = []
    if terminus == 'N':
        remove_name = ['H']
    if terminus == 'C':
        remove_name = ['O']
    for atom in residue.get_atoms():
        if atom.id == terminus:
            replace_atom = atom
        if atom.id == 'CA':
            CAatom = atom
        if atom.id in remove_name:
            remove_list.append(atom)
    # Necessary to avoid messing with very first or last res in full sequence if they were capped for simulation
    if remove_list != []:
        h_cap(replace_atom, CAatom)
        remove_atom(remove_list[0])

def get_atom_idx(entity, atom_name):
    """
    Parameters
    ----------
    :param entity:
    :type entity: Entity, required
    :param atom_name: Value corresponding to the atom attribute 'name'.
    :type atom_name: str or list, required

    Returns
    -------
    - List of indeces corresponding to atoms in entity with atom_name. 0 indexed.
    """
    idx = []
    if type(atom_name) == str:
        atom_name = [atom_name]

    for i, atom in enumerate(entity.get_atoms()):
        for name in atom_name:
            if atom.name == name:
                idx.append(i)
    return idx

def get_element_idx(entity, atom_element):
    """
    Parameters
    ----------
    :param entity:
    :type entity: Entity, required
    :param atom_element: Value corresponding to the atom attribute 'element'.
    :type atom_element: str, required

    Returns
    -------
    - List of indeces corresponding to atoms in entity of that element type. 0 indexed.
    """
    idx = []
    for i, atom in enumerate(entity.get_atoms()):
        if atom.element == atom_element:
            idx.append(i)
    return idx

def make_atom_dict(atom, idx):
    skip_keys = ['altloc' 'full_name', 'parent', '_sorting_keys', 'level', 'disordered_flag', 'anisou_array', 'siguij_array', 'sigatm_array', 'xtra']
    atom_dict = {}
    for key, val in atom.__dict__.items():
        if key in skip_keys:
            continue
        if key == 'coord':
            val = [float(c) for c in val]
        atom_dict['idx'] = idx
        atom_dict[key] = val
    return atom_dict


def make_res_dict(res, atom_count):
    _dict = {}
    skip_keys = ['altloc', 'child_list', 'child_dict', 'parent', 'level', 'disordered_flag', 'xtra']
    for key, val in res.__dict__.items():
        if key == 'child_list':
            _dict['Atoms'] = {}
            for atom in val:
                _dict['Atoms'][atom.id] = make_atom_dict(atom, atom_count)
                atom_count += 1
        if key in skip_keys:
            continue
        _dict[key] = val
    return _dict, atom_count

def make_model_dict(model):
    model_dict = {'Residues': {}}
    skip_keys = ['level', 'serial_num', 'parent', 'child_list', 'child_dict', 'xtra', 'calculations', 'method']
    atom_count = 0
    for res in model.list_residues():
        resname = res.resname
        resnumber = res.id[1]
        model_dict['Residues'][f'{resname}{resnumber}'], atom_count = make_res_dict(res, atom_count)
    for key, val in model.__dict__.items():
        if key not in skip_keys:
            model_dict[key] = val
        elif key == 'method':
            cc = []
            for res in val['catalytic_center']:
                if type(res) is str:
                    cc.append(res)
                else:
                    cc.append(f'{res.resname}{res.id[1]}')
            val['catalytic_center'] = cc
            model_dict[key] = val
        elif key == 'calculations':
            model_dict['calculations'] = []
            for c in val:
                model_dict['calculations'].append(c.__dict__)
    return model_dict