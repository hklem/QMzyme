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
    """
    Converts a QMzymeRegion to an MDanalysis AtomGroup.
    """
    print("enetered")
    attrs = {}
    top_attrs = []
    n_atoms = region.n_atoms
    coordinates = np.empty((n_atoms, 3), dtype=np.float32)

    for i, atom in enumerate(region.atoms):
        coordinates[i] = atom.position
        for attr in MDAATOMATTRIBUTES:
            if not hasattr(atom, attr):
                continue
            if attr not in attrs.keys():
                attrs[attr] = []
            attrs[attr].append(getattr(atom, attr))

    # Atom attributes
    for vals, Attr, dtype in (
        (attrs['id'], Atomids, np.int32),
        (attrs['element'], Elements, object),
        (attrs['mass'], Masses, np.float32),
        (attrs['name'], Atomnames, object),
        (attrs['type'], Atomtypes, object),
    ):
        top_attrs.append(Attr(np.array(vals, dtype=dtype)))

    # Partial charges
    if 'charge' in attrs.keys():
        top_attrs.append(Charges(np.array(attrs['charge'], dtype=np.float32)))
    else:
        pass 

    # PDB only
    for vals, Attr, dtype in (
        (attrs['chainID'], ChainIDs, object),
        (attrs['occupancy'], Occupancies, np.float32),
        (attrs['tempfactor'], Tempfactors, np.float32),
    ):
        if vals:
            top_attrs.append(Attr(np.array(vals, dtype=dtype)))

    # Residue
    if any(attrs['resnum']) and not any(val is None for val in attrs['resnum']):
        resnums = np.array(attrs['resnum'], dtype=np.int32)
        resnames = np.array(attrs['resname'], dtype=object)
        segids = np.array(attrs['segid'], dtype=object)
        icodes = np.array(attrs['icode'], dtype=object)
        residx, (resnums, resnames, icodes, segids) = change_squash(
            (resnums, resnames, icodes, segids),
            (resnums, resnames, icodes, segids))
        n_residues = len(resnums)
        for vals, Attr, dtype in (
            (resnums, Resids, np.int32),
            (resnums.copy(), Resnums, np.int32),
            (resnames, Resnames, object),
            (icodes, ICodes, object),
        ):
            top_attrs.append(Attr(np.array(vals, dtype=dtype)))
    else:
        top_attrs.append(Resids(np.array([1])))
        top_attrs.append(Resnums(np.array([1])))
        residx = None
        n_residues = 1

    if any(segids) and not any(val is None for val in segids):
        segidx, (segids,) = change_squash((segids,), (segids,))
        n_segments = len(segids)
        top_attrs.append(Segids(segids))
    else:
        n_segments = 1
        top_attrs.append(Segids(np.array(['SYSTEM'], dtype=object)))
        segidx = None

    top = Topology(n_atoms, n_residues, n_segments,
                   attrs=top_attrs, atom_resindex=residx,
                   residue_segindex=segidx)
    
    u = MDAnalysis.Universe(top, coordinates)

    for traj_attr in ['forces', 'velocities', 'charges']:
        if traj_attr in attrs.keys():
            u.add_TopologyAttr(traj_attr, attrs[traj_attr])
    ag = u.select_atoms('all')

    return u.select_atoms('all')


def atomgroup_to_region(atomgroup):
    """
    Under development.
    """
    pass
#     """
#     Note that this will make the universe based on the region, all non-region atoms will be lost.
#     """
#     try:
#         elements = atomgroup.elements
#     except:
#         from MDAnalysis.topology.guessers import guess_types
#         u = atomgroup.universe
#         guessed_elements = guess_types(u.atoms.names)
#         u.add_TopologyAttr("elements", guessed_elements)
#         warnings.warn("Element information was missing from input. MDAnalysis.topology.guessers.guess_types was used to infer element types.", UserWarning)
    
#     attrs = mda_atom_to_qmz_atom(atomgroup)

#     for i, (atom, element) in enumerate(zip(atomgroup, guessed_elements)):
#         atom = QMzymeAtom()

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
