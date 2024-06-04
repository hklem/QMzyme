import warnings
import numpy as np
from QMzyme.QMzymeAtom import QMzymeAtom
import MDAnalysis 
from MDAnalysis.core.topology import Topology
from MDAnalysis.topology.base import change_squash
from MDAnalysis.core.topologyattrs import (
    Atomids,
    Atomnames,
    Atomtypes,
    Elements,
    Masses,
    Charges,
    Bonds,
    Resids,
    Resnums,
    Resnames,
    Segids,
    ChainIDs,
    ICodes,
    Occupancies,
    Tempfactors,
)
def region_to_atom_group(region):
    
    n_atoms = region.n_atoms
    u = MDAnalysis.Universe.empty(
        n_atoms=n_atoms,  
        n_residues=n_atoms, # Although this will make u.n_residues return a misleading number, 
                            #this won't matter after the group has been saved to a PDB and reloaded into a universe.
        atom_resindex=np.arange(n_atoms), # Doing it this way makes the attribute setting simpler
        n_segments=n_atoms,
        trajectory=True) # Needs to be True so positions can be set

    # Store atom attributes
    atom_attributes = {}
    for atom in region.atoms:

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
    if 'chain' in atom_attributes:
        atom_attributes['chainID'] = atom_attributes['chain']
        del atom_attributes['chain']
    
    # Now load the attributes to the new Universe
    for attr, val in atom_attributes.items():
        # if attr in exclude_attributes:
        #     continue
        # u.add_TopologyAttr(attr, val)
        try:
            u.add_TopologyAttr(attr, val)
        except:
            pass
    u.atoms.positions = atom_attributes['position']

    # add segid info if provided
    if hasattr(region.atoms[0], "segid"):
        segments = list(set(region.segids))
        for segid in segments:
            segment = u.add_Segment(segid=segid)
            ag = []
            for atom in region.atoms:
                if atom.segid == segid:
                    ag.append(u.select_atoms(f'id {atom.id}')[0])
            sum(ag).residues.segments=segment
            #print(sum(ag).residues.segments)
            
        #segids = np.array([atom.segid for atom in region.atoms])
        #segids = region.segids
        #u.add_TopologyAttr("segids", segids)

    # Create AtomGroup and sort by resids
    atom_group = sum(list(u.atoms.sort(key='resids')))

    return atom_group


def mda_atom_to_qmz_atom(mda_atom):
    attrs = {}
    for attr in MDAATOMATTRIBUTES:
        if hasattr(mda_atom, attr):
            attrs[attr] = getattr(mda_atom, attr)
    if not hasattr(mda_atom, 'chainID'):
        attrs['chainID'] = 'X'
    atom = QMzymeAtom(**attrs)
    return atom
        
MDAATOMATTRIBUTES = [
    'chainID',
    'charge',
    'element',
    'force',
    'fragindex',
    'fragment',
    'icode',
    'id',
    'index',
    'ix',
    'ix_array',
    'mass',
    'name',
    'occupancy',
    'position',
    'radius',
    'record_type',
    'resid',
    'residue',
    'resindex',
    'resname',
    'resnum',
    'segid',
    'segindex',
    'segment',
    'tempfactor',
    'type',
    'universe',
    'velocity',
]
