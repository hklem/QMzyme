"""
Tests for the QMzyme GenerateModel.py code.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow

import QMzyme
from QMzyme.GenerateModel import GenerateModel
import pytest
import numpy as np
import numpy.testing as npt
from QMzyme.RegionBuilder import RegionBuilder
from MDAnalysis.core.universe import Universe
from QMzyme.TruncationSchemes import AlphaCarbon, TerminalAlphaCarbon, BetaCarbon
from QMzyme.data import PDB, PQR, DCD
import warnings
import os
import shutil

original_contents = os.listdir()

def restore_directory():
    for name in os.listdir():
        if name not in original_contents:
            try:
                os.remove(name)
            except:
                shutil.rmtree(name)

@pytest.mark.parametrize(
        "init_file, traj_file",
        [(PDB, None), 
         (PQR, DCD),]
)
def test_init(init_file, traj_file):
    if traj_file is not None:
        model = GenerateModel(init_file, traj_file)
        assert model.name == '1oh0_equ'
        model.set_region(name='test', selection='resid 263')
        assert hasattr(model.regions[-1].atoms[0], "charge")
        frame_0_positions = model.regions[-1].positions
        model = GenerateModel(init_file, traj_file, frame=2)
        model.set_region(name='test', selection='resid 263')
        #assert model.regions[-1].positions != frame_0_positions
        with pytest.raises(AssertionError):
            npt.assert_array_equal(model.regions[-1].positions, frame_0_positions)

    else:
        model = GenerateModel(init_file)
        assert model.name == '1oh0'
        model.set_region(name='test', selection='resid 263')
        assert not hasattr(model.regions[-1].atoms[0], "charge")
    assert model.__repr__() == "<QMzymeModel built from <Universe with 4258 atoms> contains 1 region(s)>"
    assert model.universe.__class__ == Universe
    assert model.filename == init_file

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
        model.set_region(name=region_name, selection=selection)
    elif Test == 'MDA AtomGroup as input':
        mda_atomgroup = model.universe.select_atoms(selection)
        model.set_region(name=region_name, selection=mda_atomgroup)
    elif Test == 'QMzymeRegion as input':
        mda_atomgroup = model.universe.select_atoms(selection)
        region_builder = RegionBuilder(region_name)
        region_builder.init_atom_group(mda_atomgroup)
        qmz_region = region_builder.get_region()
        model.set_region(name=region_name, selection=qmz_region)

        # check that warning is raised if you try to add a region with the same name of a region that already exists
        with pytest.raises(ValueError):
            model.set_region(name=region_name, selection=qmz_region)

    assert len(model.regions) == 1
    assert model.regions[0].name == region_name
    assert model.regions[0].n_atoms == 40
    assert model.regions[0].n_residues == 2

@pytest.mark.parametrize(
    "scheme_name, outcome",[
        (None, 'UserWarning'),
        (AlphaCarbon, None),
        (TerminalAlphaCarbon, None),
        (BetaCarbon, 'ValueError'),
    ]
)
def test_write_input(scheme_name, outcome):
    model = GenerateModel(PDB)
    model.set_region(name='test', selection='resid 3 or resid 4')
    region = model.test
    region.set_residue_method('QM')


    qm_method = QMzyme.QM_Method(
    basis_set='6-31G*',
    functional='wB97X-D3',
    qm_input='OPT FREQ',
    program='orca'
    )

    qm_method.assign_to_region(region=model.test)
    if scheme_name is not None:
        region.truncate(scheme=scheme_name, selection='resid 3 or resid 4')

    try:
        if outcome == 'UserWarning':
            with pytest.warns(UserWarning, match="The model has not been truncated"):
                model.write_input()
            assert not getattr(region, 'truncated', False)
        elif outcome == 'ValueError':
            with pytest.raises(ValueError, match="The model is only partially truncated"):
                model.write_input()
            assert not getattr(region, 'truncated', False)
        else:
            with warnings.catch_warnings():
                warnings.filterwarnings("error", message=r"The model (has not been|is only partially) truncated")
                model.write_input()
            assert region.truncated
    finally:
        restore_directory()
        assert 'QCALC' not in os.listdir()