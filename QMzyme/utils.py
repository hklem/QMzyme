###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@lsu.edu
###############################################################################

"""
Module containing broad functionality utilized throughout the package.
"""

from functools import singledispatch
from QMzyme.QMzymeModel import QMzymeModel
from QMzyme.RegionBuilder import RegionBuilder
import QMzyme.MDAnalysisWrapper as MDAwrapper
from MDAnalysis.core.groups import AtomGroup
from QMzyme.SelectionSchemes import SelectionScheme
from QMzyme.SelectionSchemes import DistanceCutoff
from abc import ABCMeta
from QMzyme.QMzymeRegion import QMzymeRegion
import numpy as np

def check_filename(filename, format):
    if filename.endswith(format):
        return filename
    if not format.startswith('.'):
            format = '.'+format
    return filename.split('.')[0]+format

@singledispatch
def make_selection(selection, model: QMzymeModel, name=None, **kwargs):
    """
    Method to enable variable input compatibility: will return an MDA AtomGroup if 
    input was an MDA selection command str, or return the input if it was either 
    an MDA AtomGroup or QMzymeRegion.
    """
    raise UserWarning(f"Invalid selection {selection}.")
    #print('make selection from: ', selection)
    #return selection

@make_selection.register
def MDA_str_selection(selection: str, model: QMzymeModel, name, **kwargs):
    region_builder = RegionBuilder(name)
    selection = MDAwrapper.select_atoms(model.universe, selection)
    region_builder.init_atom_group(selection)
    region = region_builder.get_region()
    return region

@make_selection.register
def MDA_AtomGroup_selection(selection: AtomGroup, model: QMzymeModel, name, **kwargs):
    region_builder = RegionBuilder(name)
    region_builder.init_atom_group(selection)
    region = region_builder.get_region()
    return region

@make_selection.register
def MDA_AtomGroup_selection(selection: QMzymeRegion, model: QMzymeModel, name, **kwargs):
    if name is not None:
         selection.name = name
    return selection

@make_selection.register
def scheme_selection(selection: ABCMeta, model: QMzymeModel, name, **kwargs):
    s = selection(model=model, name=name, **kwargs)
    region = s.return_region()
    return region

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

def vector_comparison(fixed, atom1, atom2, tolerance=1e-4):
    """
    Method to determine if two atoms share the same vector direction from a fixed atom.
    All atom selections are evaluated based on their Cartesian coordinate positions 
    to compute and compare normalized direction vectors.

    :param fixed: The central reference atom (Atom C) from which the vectors originate.
    :type fixed: `QMzyme.QMzymeAtom` or `MDAnalysis.Universe.Atom`
    :param atom1: The first target atom to check (Atom A).
    :type atom1: `QMzyme.QMzymeAtom` or `MDAnalysis.Universe.Atom`
    :param atom2: The second target atom to check (Atom B).
    :type atom2: `QMzyme.QMzymeAtom` or `MDAnalysis.Universe.Atom`
    :param tolerance: The absolute tolerance for comparing the dot product of the unit.
    :type tolerance: float, optional, default=1e-4
    :return: True if the vectors point in the same direction within the specified tolerance, False otherwise.
    :rtype: bool
    """
    # Extract coordinate positions
    pos_fixed = fixed.position
    pos_1 = atom1.position
    pos_2 = atom2.position
    
    # Calculate vectors originating from the fixed atom
    vec_1 = pos_1 - pos_fixed
    vec_2 = pos_2 - pos_fixed
    
    # Calculate magnitudes (norms)
    norm_1 = np.linalg.norm(vec_1)
    norm_2 = np.linalg.norm(vec_2)
    
    # Check for zero-length vectors to avoid division by zero
    if np.isclose(norm_1, 0.0) or np.isclose(norm_2, 0.0):
        raise ValueError("One of the target atoms is at the exact same position as the fixed atom.")
        
    # Normalize vectors to unit length
    unit_1 = vec_1 / norm_1
    unit_2 = vec_2 / norm_2
    
    # Compute the dot product (cos(theta))
    dot_product = np.dot(unit_1, unit_2)
    
    # If dot product is ~1.0, they share the exact same direction
    return np.isclose(dot_product, 1.0, atol=tolerance)