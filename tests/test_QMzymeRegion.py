"""
Tests for the QMzyme RegionBuilder.py ands QMzymeRegion.py code.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow


import pytest
import numpy as np
from QMzyme.QMzymeRegion import QMzymeRegion
from QMzyme.RegionBuilder import RegionBuilder
import MDAnalysis as mda
from importlib_resources import files

pdb_file = str(files('QMzyme.data').joinpath('1oh0.pdb'))
u = mda.Universe(pdb_file)
atom_group = u.select_atoms('resid 2-5')
id1 = atom_group.atoms[0].id

def test_RegionBuilder():
    rb1 = RegionBuilder(name='test')
    assert rb1.__repr__() == "<RegionBuilder: Current QMzymeRegion, test, "+\
                             "contains 0 atom(s) and 0 residue(s)>"
    
    rb1.init_atom_group(atom_group=atom_group)
    assert rb1.__repr__() == "<RegionBuilder: Current QMzymeRegion, test, "+\
                             "contains 61 atom(s) and 4 residue(s)>"

def test_QMzymeRegion():
    region_builder = RegionBuilder(name='test')
    assert region_builder.__repr__() == "<RegionBuilder: Current QMzymeRegion, test, "+\
                                        "contains 0 atom(s) and 0 residue(s)>"
    
    region_builder.init_atom_group(atom_group=atom_group)
    assert region_builder.__repr__() == "<RegionBuilder: Current QMzymeRegion, test, "+\
                                        "contains 61 atom(s) and 4 residue(s)>"
    
    # test region was populated as expected
    region = region_builder.get_region()
    assert region.__repr__() == "<QMzymeRegion test contains 61 atom(s) and 4 residue(s)>"
    assert region.has_atom(id=id1)
    assert region.has_residue(resid=3)
    assert region.atom_group == atom_group
    assert any(region.ids) == any(atom_group.ids)
    assert any(region.resids) == any(atom_group.resids)

    # add atom through region builder init_atom() method
    mda_atom = u.select_atoms('resid 1 and name CA').atoms[0]
    region_builder.init_atom(mda_atom)
    new_atom = region.get_atom(id=5)
    assert region.n_atoms == 62
    assert 1 in region.resids
    last_atom = region.atoms[-1]
    assert last_atom == new_atom

    # test sorting not in place
    sorted_region_atoms = region.sort(key='resid', in_place=False)
    assert sorted_region_atoms[0] == new_atom
    assert region.atoms[0] != new_atom

    # test sorting in place
    region.sort(key='resid', in_place=True)
    assert sorted_region_atoms[0] == new_atom
    assert region.atoms[0] == new_atom

    # now add the atom again, through QMzymeRegion add_atom() method- it will be changed because it was not unique.
    region.add_atom(new_atom)
    assert region.n_atoms == 63
    assert region.atoms[-1].id == new_atom.id+1
    assert region.atoms[-1].name == f"{new_atom.element}1"
    assert region.atoms[-1].id == max(region.get_residue(new_atom.resid).ids)

    # test setting fixed atoms
    ids = [20, 30, 60]
    region.set_fixed_atoms(ids=ids)
    for id in ids:
        assert region.get_atom(id).is_fixed

def test_QMzymeResidue():
    region_builder = RegionBuilder(name='test')
    region_builder.init_atom_group(atom_group=atom_group)
    region = region_builder.get_region()
    residue = region.residues[0]
    assert residue.__repr__() == "<QMzymeResidue resname: ASN, resid: 2, chain: A>"
    assert residue.get_atom('CA').__repr__() == "<QMzymeAtom 22: CA of resname ASN, resid 2>"
    assert residue.chain == 'A'

    
