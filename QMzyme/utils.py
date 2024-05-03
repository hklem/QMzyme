import numpy as np
from functools import singledispatch
import QMzyme.MDAnalysisWrapper as MDAwrapper

protein_residues = ['ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'CYM', 'CYS', 'CYX',
                    'GLH', 'GLN', 'GLU', 'GLY', 'HIS', 'HID', 'HIE', 'HIP',
                    'HYP', 'ILE', 'LEU', 'LYN', 'LYS', 'MET', 'PHE', 'PRO',
                    'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HSE', 'HSD', 'HSP',
                    'SEC', 'PYL']

def set_bond_length(mobile_coords, fixed_coords, new_length):
    M = new_length/np.linalg.norm(fixed_coords-mobile_coords)
    new_coords = fixed_coords-(M*(fixed_coords-mobile_coords))
    return new_coords

@singledispatch
def translate_selection(selection, model):
    """
    Method to enable variable input comability: will return an MDA AtomGroup if 
    input was an MDA selection command str, or return the input if it was either 
    an MDA AtomGroup or QMzymeRegion.
    """
    return selection
    
@translate_selection.register(str)
def convert_MDA_selection(selection, model):
    return MDAwrapper.select_atoms(model.starting_structure, selection)
