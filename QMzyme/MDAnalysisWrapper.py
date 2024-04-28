"""
Code to integrate MDAnalysis utilities in QMzyme. 
"""
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.pkdtree import *
from QMzyme import utils

def init_universe(topology, traj=None):
    if traj is None:
        return mda.Universe(topology)
    else:
        return mda.Universe(topology, traj)

def select_atoms(universe, selection):
    """
    :param universe: MDAnalysis Universe object.
    :param selection: Selection of atoms to be made- based on `MDAnalysis selection command language <https://docs.mdanalysis.org/stable/documentation_pages/selections.html>`_.
    :type selection: str, required
    """
    u = universe.select_atoms(selection)
    return u

def concatenate_atoms(ag1, ag2):
    """
    To create an atom group with unique atoms (no duplicates) from two atom groups.
    """
    return

def get_neighbors(ag1, ag2, radius, remove_duplicates=True):
    """
    Returns list of atoms in distance based atom group.
    """
    # Using the fast C based code
    tree = PeriodicKDTree()
    atoms = []
    full_system = ag1
    sub_system = ag2

    # To ensure bigger atom group selection is used to set_coords
    if len(ag2) > len(ag1):
        full_system = ag2
        sub_system = ag1
    tree.set_coords(full_system.positions)
    pairs = tree.search_tree(sub_system.positions, radius)

    for pair in pairs:
        atom = full_system[pair[1]]
        if remove_duplicates is True:
            if atom not in atoms:
                atoms.append(atom)
        else:
            atoms.append(atom)
    
    return atoms

def order_residues(residues, remove_duplicates=True):
    """
    :param residues: 
    :type residues: list of <Residue> objects, required

    Returns list of <Residue> objects ordered by increasing resid.
    """
    if remove_duplicates is True:
        return sorted(set(residues))
    else:
        return sorted(residues)

def CA_cap(residue, cap_atom_name, Hbond_length=1.09):
    """
    :param atom: The Atom to be replaced with H.
    Hbond_length defaults to equilibrium C-H bond length of 1.09.
    """
    for atom in residue.atoms:
        if atom.name == 'CA':
            CA_atom = atom
        elif atom.name == cap_atom_name:
            cap_atom = atom
    # Try statement in case the cap_atom doesn't exist
    try:       
        coords = utils.set_bond_length(mobile_coords=cap_atom.position, fixed_coords=CA_atom.position, new_length=Hbond_length)
        new_atom = {
            'element': 'H', 
            'name': 'Hcap', 
            'position': coords, 
            'mass': 1.00794,
        }
        alter_atom(cap_atom, new_atom)
    except:
        pass

def get_atom(residue, atom):
    if type(atom) is str:
        return residue.atoms[list(residue.atoms.names).index(atom)]
    if type(atom) is int:
        return residue.atoms[list(residue.atoms.ids).index(atom)]

def name_cap_H(residue):
    if 'H1' not in residue.atoms.names:
        name = 'H1'
    elif 'H2' not in residue.atoms.names:
        name = 'H2'
    else:
        name = 'H3'
    return name

def cap_backbone_CA(target_atom, bond_length=1.09, cap='H'):
    """
    :param target_atom: The atom bound to CA that will be converted.
    """
    CA_atom = get_atom(target_atom.residue, 'CA')
    if cap == 'H':
        return cap_H(target_atom, CA_atom, bond_length)
    
def cap_backbone_N(N_atom, bond_length=1.01, cap='H'):
    """
    :param N_atom: The backbone N Atom bound to the preceeding residue backbone C Atom that will be converted.
    """
    C_atom = get_atom(get_preceding_residue(N_atom.residue), 'C')
    if cap == 'H':
        return cap_H(C_atom, N_atom, bond_length)

def cap_backbone_C(C_atom, bond_length=1.09, cap='H'):
    """
    :param C_atom: The backbone C Atom bound to the next residue backbone N Atom that will be converted.
    :type C_atom: Universe.Atom, required
    """
    print(list(get_next_residue(C_atom.residue).atoms))
    N_atom = get_atom(get_next_residue(C_atom.residue), 'N')
    if cap == 'H':
        return cap_H(N_atom, C_atom, bond_length)

def cap_H(atom, fixed_atom, bond_length):
    """
    :param atom: The Atom to be converted to hydrogen.
    :param fixed_atom: The Atom bound to atom that serves as anchor point.
    """
    res = atom.residue
    new_position = utils.set_bond_length(atom.position, fixed_atom.position, bond_length)
    new_name = name_cap_H(res)
    new_atom_dict = {
            'element': 'H', 
            'name': new_name, 
            'position': new_position, 
            'mass': 1.00794,
        }
    new_atom = alter_atom(atom, new_atom_dict)
    return res.atoms.select_atoms(f'name {new_atom.name}')

def cap_proline_H(N_atom, bond_length=1.01):
    u = N_atom.universe
    resindex = N_atom.residue.resindex
    prev_res = u.residues[resindex-1]
    prev_C_atom = prev_res.atoms.select_atoms('name C')[0]
    new_position = utils.set_bond_length(prev_C_atom.position, N_atom.position, bond_length)
    new_name = name_cap_H(N_atom.residue)
    new_atom_dict = {
            'element': 'H', 
            'name': new_name, 
            'position': new_position, 
            'mass': 1.00794,
            'resid': N_atom.residue.resid,
            'residue': N_atom.residue,
            #'resindex': resindex,
            'resname': N_atom.residue.resname,
            'resnum': N_atom.residue.resnum
        }
    new_atom = alter_atom(prev_C_atom, new_atom_dict)
    return N_atom.residue.atoms.select_atoms(f'name {new_atom.name}')

def alter_atom(atom, new_atom_dict):
    for key, val in new_atom_dict.items():
        if 'res' in key:
            setattr(atom.residue, key, val)
        elif hasattr(atom, key):
            setattr(atom, key, val)
    return atom

def check_seq_neighbor(residue_list, residue, increment):
    resid = residue.resid
    resids = [res.resid for res in residue_list]
    if resid + increment in resids:
        return True
    else:
        return False
    
def has_Cterm_neighbor(residue_list, residue):
    return check_seq_neighbor(residue_list, residue, 1)

def has_Nterm_neighbor(residue_list, residue):
    return check_seq_neighbor(residue_list, residue, -1)

def write_pdb(selection, filename):
    selection.write(filename)

def get_atom_idx(atom_group, name_list):
    ids = []
    for i, atom in enumerate(atom_group):
        if atom.name in name_list:
            ids.append(i)
    return ids

def closest_waters(full_system, sub_system, Nwaters):
    pass

def get_preceding_residue(residue):
    u = residue.universe
    ids = [res.resid for res in u.residues]
    r = u.residues[ids.index(residue.resid-1)]
    return r
    
def get_next_residue(residue):
    u = residue.universe
    ids = [res.resid for res in u.residues]
    r = u.residues[ids.index(residue.resid+1)]
    return r


def select_sidechain(residue):
    return residue.atoms.select_atoms("not backbone and not name H")


def get_residue(residue, other_universe):
    resid = residue.resid
    for res in other_universe.residues:
        if res.resid == resid:
            return res
    
def get_parallel_residue(residue, other_universe):
    resid = residue.resid
    for res in other_universe.residues:
        if res.resid == resid:
            return res

def get_parallel_atom(atom, other_universe):
    atom_id = atom.id
    for atom in other_universe.atoms:
        if atom.id == atom_id:
            return atom



def fix_pdb(original_universe, pdb_file):

    with open(pdb_file) as f:
        data = f.readlines()
    u = mda.Universe(pdb_file)
    for res in u.residues:
        if 'N' not in res.atoms.names:
            og_res = get_residue(res, original_universe)
            a = get_atom(og_res, 'N')
            new_name = name_cap_H(N_atom.residue)
            new_atom_dict = {
                'element': 'H', 
                'name': new_name, 
                'position': new_position, 
                'mass': 1.00794,
            }


            alter_atom()

    for line in data:
        if line.startswith('ATOM') is False:
            continue
        resid = int(line.split()[5])



