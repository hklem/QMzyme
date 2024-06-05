"""
Tests for the QMzymeModel.py code.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow

import pytest
import os
import shutil
from QMzyme.QMzymeModel import QMzymeModel
from QMzyme.RegionBuilder import RegionBuilder
from QMzyme.data import PDB
from QMzyme.MDAnalysisWrapper import init_universe


u = init_universe(PDB)
original_contents = os.listdir()

def restore_directory():
    for name in os.listdir():
        if name not in original_contents:
            try:
                os.remove(name)
            except:
                shutil.rmtree(name)

def test_QMzymeModel():

    # instantiate model and ensure base attributes are set
    model = QMzymeModel(name='1oh0', universe=u)
    assert model.name == '1oh0'
    assert model.__repr__() == "<QMzymeModel 1oh0 built from <Universe with 4258 atoms> contains 0 region(s)>"
    assert model.filename == PDB
    assert model.n_regions == 0
    assert hasattr(model, 'universe')

    # build region and add to model
    rb1 = RegionBuilder(name='test_region')
    rb1.init_atom_group(atom_group=u.select_atoms('resid 263'))
    region = rb1.get_region()
    model.add_region(region)

    # check model has that region now
    assert model.has_region('test_region') 
    assert model.n_regions == 1

    # check region related methods
    assert model.get_region_names() == ['test_region']
    assert model.get_region(region_name='test_region') == region

    # negative test
    with pytest.raises(UserWarning):
        model.get_region('blah')

    # check writing of pymol visualization script
    model.pymol_visualize()
    assert f'QMzymeModel_{model.name}_visualize.py' in os.listdir()

    restore_directory()


