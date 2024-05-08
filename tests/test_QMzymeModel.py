"""
Tests for the QMzymeModel.py code.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow

import pytest
from QMzyme.QMzymeModel import QMzymeModel
from QMzyme.RegionBuilder import RegionBuilder
from importlib_resources import files
from QMzyme.MDAnalysisWrapper import init_universe

pdb_file = str(files('QMzyme.data').joinpath('1oh0.pdb'))
u = init_universe(pdb_file)

def test_QMzymeModel():
    model = QMzymeModel(name='test', universe=u)
    assert model.name == 'test'
    assert model.__repr__() == "<QMzymeModel test built from <Universe with 4258 atoms> contains 0 region(s)>"
    assert model.filename == pdb_file
    assert model.n_regions == 0

    rb1 = RegionBuilder(name='test_region')
    rb1.init_atom_group(atom_group=u.select_atoms('resid 263'))
    region = rb1.get_region()
    model.add_region(region)
    assert model.has_region('test_region') 
    assert model.n_regions == 1
    assert model.get_region_names() == ['test_region']
    assert model.get_region(region_name='test_region') == region
    with pytest.raises(UserWarning):
        model.get_region('blah')




