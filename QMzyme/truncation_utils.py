import numpy as np
from QMzyme.QMzymeAtom import QMzymeAtom
from QMzyme.RegionBuilder import RegionBuilder


def cap_H(replace_atom, fixed_atom, bond_length=1.09):
    """
    :param replace_atom: The Atom to be converted to hydrogen.
    :param fixed_atom: The Atom bound to replace_atom that serves as reference point for bond vector calculation.
    """
    new_position = set_bond_length(replace_atom.position, fixed_atom.position, bond_length)
    new_atom_dict = {
            'element': 'H',
            'type': 'H',  
            'name': 'H1',
            'position': new_position, 
            'mass': 1.00794,
        }
    new_atom = create_new_atom(fixed_atom, new_atom_dict) # used fixed atom because sometimes replaced atom comes from a different residue
    return new_atom


def set_bond_length(mobile_coords, fixed_coords, new_length):
    M = new_length/np.linalg.norm(fixed_coords-mobile_coords)
    new_coords = fixed_coords-(M*(fixed_coords-mobile_coords))
    return new_coords


# def name_cap_H(region, resid, name='H1'):
#     residue_atoms = region.get_residue_atoms(resid)
#     i = 1
#     while name in [atom.name for atom in residue_atoms]:
#         i += 1
#         name = f'H{i}'
#     return name


def create_new_atom(atom, new_atom_dict):
    for key, val in atom.__dict__.items():
        if key.startswith('_QMzymeAtom__'):
            key = key.split('_QMzymeAtom__')[-1]
        if key not in new_atom_dict:
            new_atom_dict[key] = val
    new_atom = QMzymeAtom(**new_atom_dict)
    return new_atom


def get_preceding_Catom(region, resid):
    mda_atom = region.get_atom_group().universe.select_atoms(f'resid {resid-1} and name C')
    try:
        atom = RegionBuilder('temp', mda_atom).get_region().atoms[0]
    except:
        atom = None
    return atom


def get_following_Natom(region, resid):
    mda_atom = region.get_atom_group().universe.select_atoms(f'resid {resid+1} and name N')
    try:
        atom = RegionBuilder('temp', mda_atom).get_region().atoms[0]
    except:
        atom = None
    return atom


def has_Nterm_neighbor(atom):
    return atom.resid-1 in atom.region.resids


def has_Cterm_neighbor(atom):
    return atom.resid+1 in atom.region.resids
