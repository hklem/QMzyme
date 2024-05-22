"""
Tests for the QMzyme GenerateModel.py code.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow

from QMzyme.GenerateModel import GenerateModel
import pytest
from QMzyme.RegionBuilder import RegionBuilder
from MDAnalysis.core.universe import Universe
#from importlib_resources import files
from QMzyme.data import PDB


def test_init():
    model = GenerateModel(PDB)
    assert model.__repr__() == "<QMzymeModel built from <Universe with 4258 atoms> contains 0 region(s)>"
    assert model.name == '1oh0'
    assert model.universe.__class__ == Universe
    assert model.filename == PDB
    assert model.regions == []

def test_set_catalytic_center(selection='resid 263'):
    model = GenerateModel(PDB)
    model.set_catalytic_center(selection)
    assert len(model.regions) == 1
    assert model.regions[0].name == 'catalytic_center'
    assert model.regions[0].n_atoms == 37

    model.remove_region("catalytic_center")
    assert len(model.regions) == 0


selection_str = 'resid 16 or resid 17'
@pytest.mark.parametrize(
        "Test, init_file, region_name, selection",
        [('Selection string as input', PDB, 'test', selection_str), 
         ('MDA AtomGroup as input', PDB, 'test', selection_str), 
         ('QMzymeRegion as input', PDB, 'test', selection_str),]
)
def test_set_region(Test, init_file, region_name, selection):
    model = GenerateModel(init_file)
    if Test == 'Selection string as input':
        model.set_region(region_name=region_name, selection=selection)
    elif Test == 'MDA AtomGroup as input':
        mda_atomgroup = model.universe.select_atoms(selection)
        model.set_region(region_name=region_name, selection=mda_atomgroup)
    elif Test == 'QMzymeRegion as input':
        mda_atomgroup = model.universe.select_atoms(selection)
        #qmz_region = RegionBuilder(region_name, mda_atomgroup).get_region()
        region_builder = RegionBuilder(region_name)
        region_builder.init_atom_group(mda_atomgroup)
        qmz_region = region_builder.get_region()
        model.set_region(region_name=region_name, selection=qmz_region)
    assert len(model.regions) == 1
    assert model.regions[0].name == region_name
    assert model.regions[0].n_atoms == 40
    assert model.regions[0].n_residues == 2
