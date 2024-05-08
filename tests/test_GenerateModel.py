"""
Tests for the QMzyme GenerateModel.py code.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow

import os
from QMzyme.GenerateModel import GenerateModel
import pytest
from QMzyme.QMzymeAtom import QMzymeAtom
from QMzyme.RegionBuilder import RegionBuilder
from MDAnalysis.core.universe import Universe
import MDAnalysis as MDA
from importlib_resources import files

pdb_file = str(files('QMzyme.data').joinpath('1oh0.pdb'))

def test_init(init_file=pdb_file):
    model = GenerateModel(file=init_file)
    assert model.name == '1oh0'
    assert model.starting_structure.__class__ == Universe
    assert model.filename == pdb_file
    assert model.regions == []

def test_set_catalytic_center(init_file=pdb_file, selection='resid 263'):
    model = GenerateModel(file=init_file)
    model.set_catalytic_center(selection)
    assert len(model.regions) == 1
    assert model.regions[0].name == 'catalytic_center'
    assert model.regions[0].n_atoms == 37


selection_str = 'resid 16 or resid 17'
@pytest.mark.parametrize(
        "Test, init_file, region_name, selection",
        [('Selection string as input', pdb_file, 'test', selection_str), 
         ('MDA AtomGroup as input', pdb_file, 'test', selection_str), 
         ('QMzymeRegion as input', pdb_file, 'test', selection_str),]
)
def test_set_region(Test, init_file, region_name, selection):
    model = GenerateModel(file=init_file)
    if Test == 'Selection string as input':
        model.set_region(region_name, selection)
    elif Test == 'MDA AtomGroup as input':
        mda_atomgroup = model.starting_structure.select_atoms(selection)
        model.set_region(region_name, mda_atomgroup)
    elif Test == 'QMzymeRegion as input':
        mda_atomgroup = model.starting_structure.select_atoms(selection)
        #qmz_region = RegionBuilder(region_name, mda_atomgroup).get_region()
        region_builder = RegionBuilder(region_name)
        region_builder.init_atom_group(mda_atomgroup)
        qmz_region = region_builder.get_region()
        model.set_region(region_name, qmz_region)
    assert len(model.regions) == 1
    assert model.regions[0].name == region_name
    assert model.regions[0].n_atoms == 40
    assert model.regions[0].n_residues == 2
