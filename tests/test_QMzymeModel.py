"""
Tests for the QMzymeModel.py code.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow

import pytest
import os
import shutil
import QMzyme
from QMzyme.QMzymeModel import QMzymeModel
from QMzyme.RegionBuilder import RegionBuilder
from QMzyme.data import PDB, PDB_xtal, LIG
from QMzyme.MDAnalysisWrapper import init_universe
from QMzyme.SelectionSchemes import DistanceCutoff


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

    # see the overview
    model.print_overview()

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


def test_import_region():
    model = QMzyme.GenerateModel(PDB_xtal)

    QMzyme.data.residue_charges.update({'EQU': -1}) # EQU ligand

    model.set_catalytic_center(selection='resname EQU and segid A')
    model.set_region(selection=DistanceCutoff, cutoff=3)

    assert len(model.regions) == 2

    region_19nt = model.import_region(LIG, name='region_19nt')
    QMzyme.data.residue_charges.update({'6VW': 0})

    assert len(model.regions) == 3
    assert model.region_19nt.n_residues == 1
    assert model.region_19nt.n_atoms == 46