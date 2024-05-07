"""
Code to integrate MDAnalysis utilities in QMzyme. 
"""
import os
import numpy as np
import warnings
import MDAnalysis as mda
from MDAnalysis.lib.pkdtree import *
from QMzyme import utils
from MDAnalysis.core.universe import Universe


def init_universe(topology, **kwargs):
    u = mda.Universe(topology, **kwargs)
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
        
def get_atom(residue, atom):
    if type(atom) is str:
        return residue.atoms[list(residue.atoms.names).index(atom)]
    if type(atom) is int:
        return residue.atoms[list(residue.atoms.ids).index(atom)]


def alter_atom(atom, new_atom_dict):
    for key, val in new_atom_dict.items():
        if 'res' in key:
            setattr(atom.residue, key, val)
        elif hasattr(atom, key):
            setattr(atom, key, val)
    return atom


def name_cap_H(residue):
    if 'H1' not in residue.atoms.names:
        name = 'H1'
    elif 'H2' not in residue.atoms.names:
        name = 'H2'
    else:
        name = 'H3'
    return name


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
            'type': 'H',  
            'name': new_name, 
            'position': new_position, 
            'mass': 1.00794,
        }
    new_atom = alter_atom(atom, new_atom_dict)
    return new_atom


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
    N_atom = get_atom(get_next_residue(C_atom.residue), 'N')
    if cap == 'H':
        return cap_H(N_atom, C_atom, bond_length)
    

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
    

def build_universe_from_QMzymeRegion(region):
    
    n_atoms = region.n_atoms
    u = mda.Universe.empty(
        n_atoms=n_atoms,  
        n_residues=n_atoms, # Although this will make u.n_residues return a misleading number, 
                            #this won't matter after the group has been saved to a PDB and reloaded into a universe.
        atom_resindex=np.arange(n_atoms), # Doing it this way makes the attribute setting simpler
        trajectory=True) # Needs to be True so positions can be set

    # Store atom attributes
    atom_attributes = {}
    for atom in region.atoms:
        #for attr in dir(atom):
            #if attr.startswith('_') or attr.startswith('get'):
        for attr in atom.__dict__:
            if attr.startswith('_'):
                continue
            elif attr not in atom_attributes.keys():
                try:
                    atom_attributes[attr] = [getattr(atom, attr)]
                except:
                    pass
            else:
                atom_attributes[attr].append(getattr(atom, attr))
    
    exclude_attributes = [ #attributes that can't be set as topology level attributes
        'index', 'ix', 'ix_array', 'level', 'position', 'residue', 
        'resindex', 'segid', 'segindex', 'segment', 'universe', 'region'] 
    
    if 'chain' in atom_attributes:
        atom_attributes['chainID'] = atom_attributes['chain']
        del atom_attributes['chain']
    
    # Now load the attributes to the new Universe
    for attr, val in atom_attributes.items():
        if attr in exclude_attributes:
            continue
        u.add_TopologyAttr(attr, val)
    u.atoms.positions = atom_attributes['position']

    # Create AtomGroup and sort by resids
    atom_group = sum(list(u.atoms.sort(key='resids')))

    return atom_group


def build_universe(atom_list, save_pdb=False, filename=None):
    """
    Function to combine a list of atoms into one universe, regardless of if they are from different universes.

    Parameters
    -----------
    :param atom_list: List containing MDAnalysis Atom objects that can be from different universes. 
    :type atom_list: List[Atom]

    Returns
    --------
    MDAnalysis Universe
    """
    # suppress some MDAnalysis warnings when writing PDB files
    warnings.filterwarnings('ignore')

    # First, make sure no atoms repeat in atom_list
    atom_list = np.unique(atom_list)

    # Init empty Universe
    n_atoms = len(atom_list)
    u = mda.Universe.empty(
        n_atoms=n_atoms,  
        n_residues=n_atoms, # Although this will make u.n_residues return a misleading number, 
                            #this won't matter after the group has been saved to a PDB and reloaded into a universe.
        atom_resindex=np.arange(n_atoms), # Doing it this way makes the attribute setting simpler
        trajectory=True) # Needs to be True so positions can be set

    # Store atom attributes
    atom_attributes = {}
    for atom in atom_list:
        for attr in dir(atom):
            if attr.startswith('_') or attr.startswith('get'):
                continue
            elif attr not in atom_attributes.keys():
                try:
                    atom_attributes[attr] = [getattr(atom, attr)]
                except:
                    pass
            else:
                atom_attributes[attr].append(getattr(atom, attr))

    exclude_attributes = [ #attributes that can't be set as topology attributes
        'index', 'ix', 'ix_array', 'level', 'position', 'residue', 
        'resindex', 'segid', 'segindex', 'segment', 'universe'] 
    
    # Now load the attributes to the new Universe
    for attr, val in atom_attributes.items():
        if attr in exclude_attributes:
            continue
        u.add_TopologyAttr(attr, val)
    u.atoms.positions = atom_attributes['position']
    u.dimensions = atom_list[0].universe.dimensions # Avoids MDAnalysis raising a warning because PDB format requires this.

    # Create AtomGroup and sort by resids
    atom_group = sum(list(u.atoms.sort(key='resids')))

    # Save as PDB, then reload to fix n_res inconsistencies
    if filename is None:
        filename='temp.pdb'
    atom_group.write(filename)
    u = mda.Universe(filename)

    # Clean up
    if save_pdb is False:
        os.remove(filename)

    return u

def is_universe(object):
    return isinstance(object, Universe)

