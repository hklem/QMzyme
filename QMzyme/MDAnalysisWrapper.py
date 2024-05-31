###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
Code to integrate MDAnalysis utilities in QMzyme. 
"""
import numpy as np
import warnings
import MDAnalysis as mda
from MDAnalysis.lib.pkdtree import *


def init_universe(*args, frame=0, **kwargs):
    """
    Accepts all argument and key word arguments that :class:`~MDAnalysis.Universe`
    can accept to create a Universe instance. Note, you may need to pass the 
    format key word in some cases. 
    """
    u = mda.Universe(*args, **kwargs)
    if frame != 0:
            u.trajectory[frame]
    if not hasattr(u.atoms, "elements"):
        from MDAnalysis.topology.guessers import guess_types
        guessed_elements = guess_types(u.atoms.names)
        u.add_TopologyAttr("elements", guessed_elements)
        warnings.warn("Element information was missing from input. MDAnalysis.topology.guessers.guess_types was used to infer element types.", UserWarning)
    return u

def select_atoms(universe, selection):
    """
    :param universe: MDAnalysis Universe object.
    :param selection: Selection of atoms to be made- based on `MDAnalysis selection command language <https://docs.mdanalysis.org/stable/documentation_pages/selections.html>`_.
    :type selection: str, required
    """
    u = universe.select_atoms(selection)
    return u


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
    
    return sum(atoms)
