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
    
    region = region_builder.get_region()
    assert region.__repr__() == "<QMzymeRegion test contains 61 atom(s) and 4 residue(s)>"
    assert region.atom_group == atom_group
    assert any(region.ids) == any(np.arange(20,81))
    assert any(region.resids) == any(np.arange(2,6))

    mda_atom = u.select_atoms('resid 1 and name CA').atoms[0]
    region_builder.init_atom(mda_atom)
    new_atom = region.get_atom(id=5)
    #last_atom = region.atoms[-1]
    assert 1 in region.resids
    last_atom = region.atoms[-1]
    assert region.atoms[-1] == new_atom

    sorted_region_atoms = region.sort(key='resid', in_place=False)
    assert sorted_region_atoms[0] == new_atom
    assert region.atoms[0] != new_atom

    region.sort(key='resid', in_place=True)
    assert sorted_region_atoms[0] == new_atom
    assert region.atoms[0] == new_atom





