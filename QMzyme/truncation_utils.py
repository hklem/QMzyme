###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
Module containing functions utilized by truncation_schemes.py.
"""

import numpy as np
from QMzyme.QMzymeAtom import QMzymeAtom
from QMzyme.RegionBuilder import RegionBuilder
from QMzyme.converters import *


def cap_H(replace_atom, fixed_atom, bond_length=1.09, base_atom=None):
    """
    :param replace_atom: The Atom to be converted to hydrogen.
    :param fixed_atom: The Atom bound to replace_atom that serves as reference point for bond vector calculation.
    """
    new_position = set_bond_length(replace_atom.position, fixed_atom.position, bond_length)
    new_atom_dict = {
            'element': 'H',
            'type': 'H',  
            'name': f'H{replace_atom.element}',
            'position': new_position, 
            'mass': 1.00794,
        }
    if base_atom is None:
        base_atom = replace_atom
        if fixed_atom.resname == 'PRO' and fixed_atom.name == 'N':
            new_atom_dict['charge'] = 0.0
            new_atom_dict['name'] = 'HN'
            base_atom = fixed_atom
    new_atom = create_new_atom(base_atom, new_atom_dict) # used fixed atom because sometimes replaced atom comes from a different residue
    return new_atom

def balance_charge(region, truncated_region):
    """
    To be used if original region has atom charge information, 
    so the QMzymeRegion guess_charge() and read_charges() methods do not get messed up.
    Within a residue, any atoms added should distribute the charges of any atoms removed. 
    """
    for res in region.residues:
        truncated_res = truncated_region.get_residue(res.resid)
        removed_atoms = []
        added_atoms = []
        chrg = 0
        for atom in res.atoms:
            if atom.name not in [atom.name for atom in truncated_res.atoms]:
                removed_atoms.append(atom)
        for atom in truncated_res.atoms:
            if atom.name not in [atom.name for atom in res.atoms]:
                added_atoms.append(atom)
        if len(removed_atoms) == 0:
            continue
        for atom in removed_atoms:
            chrg += atom.charge
        fractional_charge = chrg/len(added_atoms)
        for atom in added_atoms:
            atom.charge = fractional_charge


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


def create_new_atom(base_atom, new_atom_dict):
    for key, val in base_atom.__dict__.items():
        if key.startswith('_QMzymeAtom__'):
            key = key.split('_QMzymeAtom__')[-1]
        if key not in new_atom_dict:
            new_atom_dict[key] = val
    new_atom = QMzymeAtom(**new_atom_dict)
    return new_atom


def get_preceding_Catom(region, resid):
    if resid == 1:
        return None
    if region._universe != None:
        mda_atom = region._universe.select_atoms(f'resid {resid-1} and name C').atoms[0]
        atom = mda_atom_to_qmz_atom(mda_atom)
    return atom


def get_following_Natom(region, resid):
    try:
        if region._universe != None:
            mda_atom = region._universe.select_atoms(f'resid {resid+1} and name N').atoms[0]
        atom = mda_atom_to_qmz_atom(mda_atom)
    except:
        atom = None # covers if res is last protein res in universe
    return atom


def has_Nterm_neighbor(atom):
    return atom.resid-1 in atom.region.resids


def has_Cterm_neighbor(atom):
    return atom.resid+1 in atom.region.resids
