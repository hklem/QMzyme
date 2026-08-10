###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@lsu.edu
###############################################################################

"""
Code to integrate MDAnalysis utilities in QMzyme. 
"""
import numpy as np
import pickle
import warnings
import MDAnalysis as mda
from MDAnalysis.lib.pkdtree import *


def init_universe(*args, frame=0, **kwargs):
    """
    Accepts all argument and key word arguments that :class:`~MDAnalysis.Universe`
    can accept to create a Universe instance. Note, you may need to pass the 
    format key word in some cases. 
    """
    warnings.filterwarnings(
        action='ignore',
        category=DeprecationWarning
        )
    u = mda.Universe(*args, **kwargs)
    if frame != 0:
            u.trajectory[frame]
    if not hasattr(u.atoms, "elements"):
        from MDAnalysis.topology.guessers import guess_types
        guessed_elements = guess_types(u.atoms.names)
        u.add_TopologyAttr("elements", guessed_elements)
        warnings.warn("Element information was missing from input. MDAnalysis.topology.guessers.guess_types was used to infer element types.", UserWarning, stacklevel=2)
    if not hasattr(u.atoms, 'chainID') or u.atoms[0].chainID == '':
        u.add_TopologyAttr("chainID")
        u.atoms.chainIDs = 'X'
    return u

def select_atoms(universe, selection):
    """
    :param universe: MDAnalysis Universe object.
    :param selection: Selection of atoms to be made- based on `MDAnalysis selection command language <https://docs.mdanalysis.org/stable/documentation_pages/selections.html>`_.
    :type selection: str, required
    """
    atom_group = universe.select_atoms(selection)
    return atom_group

def universe_selection(universe, selection):
    sel = universe.select_atoms(selection)
    new_universe = mda.Merge(sel.atoms)

    if hasattr(universe, 'trajectory') and universe.trajectory.n_frames > 1:
        current_frame = universe.trajectory.ts.frame
        coordinates = np.array([
            sel.positions.copy() for ts in universe.trajectory
        ])
        universe.trajectory[current_frame]
        new_universe.load_new(coordinates, format=mda.coordinates.memory.MemoryReader)
        new_universe.trajectory[current_frame]

    return new_universe

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

def find(atom_group, name):
    """
    Look up a single atom by name within an AtomGroup.
 
    :param atom_group: MDAnalysis AtomGroup to search.
    :param name: Atom name to search for.
    :returns: MDAnalysisAtom object.
    """
    sel = atom_group.select_atoms(f"name {name}")
    return sel[0] if len(sel) else None