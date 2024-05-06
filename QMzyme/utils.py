import numpy as np
from functools import singledispatch
import QMzyme.MDAnalysisWrapper as MDAwrapper
from QMzyme.QMzymeRegion import QMzymeRegion
from MDAnalysis.core.groups import AtomGroup


protein_residues = ['ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'CYM', 'CYS', 'CYX',
                    'GLH', 'GLN', 'GLU', 'GLY', 'HIS', 'HID', 'HIE', 'HIP',
                    'HYP', 'ILE', 'LEU', 'LYN', 'LYS', 'MET', 'PHE', 'PRO',
                    'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HSE', 'HSD', 'HSP',
                    'SEC', 'PYL']


@singledispatch
def translate_selection(selection, universe):
    """
    Method to enable variable input comability: will return an MDA AtomGroup if 
    input was an MDA selection command str, or return the input if it was either 
    an MDA AtomGroup or QMzymeRegion.
    """
    return selection

@translate_selection.register
def convert_MDA_selection(selection: str, universe):
    return MDAwrapper.select_atoms(universe, selection)

@translate_selection.register
def convert_MDA_selection(selection: QMzymeRegion, universe):
    return selection

@translate_selection.register
def convert_MDA_selection(selection: AtomGroup, universe):
    return selection
