###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
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
    Method to enable variable input comability: will return an MDA AtomGroup if 
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
