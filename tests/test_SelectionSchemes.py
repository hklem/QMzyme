"""
Tests for the QMzyme GenerateModel.py code.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow

from QMzyme.GenerateModel import GenerateModel
import pytest
from QMzyme.RegionBuilder import RegionBuilder
from MDAnalysis.core.universe import Universe
from QMzyme.data import PDB, CSA_holo, CSA_apo, CSA_pkl, Cutoff_3
import QMzyme
import pickle
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

def test_DistanceCutoff():
    from QMzyme.SelectionSchemes import DistanceCutoff
    model = GenerateModel(PDB)
    with pytest.raises(UserWarning):
        DistanceCutoff(model=model, name=None, cutoff=3, include_whole_residues=True)
    with pytest.raises(UserWarning):
        model.set_region(selection=DistanceCutoff, name=None, cutoff=3, include_whole_residues=True)
    model.set_catalytic_center('resid 263')
    model.set_region(selection=DistanceCutoff, name=None, cutoff=3, include_whole_residues=True)
    assert len(model.regions) == 2
    assert model.cutoff_3.n_atoms == 275
    assert model.cutoff_3.n_residues == 19

    model.set_region(selection=DistanceCutoff, name=None, cutoff=5, include_whole_residues=True)
    assert len(model.regions) == 3
    assert model.cutoff_5.n_atoms == 427
    assert model.cutoff_5.n_residues == 33

def test_CSACutoff():
    from QMzyme.SelectionSchemes import CSACutoff
    model = GenerateModel(PDB)

    # Initial check to see the errors derived from lack of initial inputs
    with pytest.raises(UserWarning):
        CSACutoff(model=model, name=None, min_atoms=None, max_atoms=None, memory=None, nprocs=None, method=None)
    with pytest.raises(UserWarning):
        CSACutoff(model=model, name=None, method=None, holo_output_files=None, apo_output_files=None, pop="hirshfeld", charge_threshold=0.05)
    with pytest.raises(UserWarning):
        model.set_region(selection=CSACutoff, name=None, min_atoms=None, max_atoms=None, memory=None, nprocs=None, method=None)
    with pytest.raises(UserWarning):
        model.set_region(selection=CSACutoff, name=None, method=None, holo_output_files=None, apo_output_files=None, pop="hirshfeld", charge_threshold=0.05)

    # Error from lack of catalytic center
    qm_method = QMzyme.QM_Method(
        basis_set='6-31G*',
        functional='wB97XD',
        qm_input='pop=hirshfeld',
        program='gaussian'
    ) 
    with pytest.raises(UserWarning):
        model.set_region(selection=CSACutoff, name=None, min_atoms=950, max_atoms=1050, memory=None, nprocs=None, method=qm_method)

    # Error from lack of method
    model = GenerateModel(PDB)
    model.set_catalytic_center('resid 263')
    with pytest.raises(UserWarning):
        model.set_region(selection=CSACutoff, name=None, min_atoms=None, max_atoms=None, memory=None, nprocs=None, method=None)
    with pytest.raises(UserWarning):
        model.set_region(selection=CSACutoff, name=None, method=None, holo_output_files=CSA_holo, apo_output_files=CSA_apo, pop="hirshfeld", charge_threshold=0.05)

    qm_method = QMzyme.QM_Method(
        basis_set='6-31G*',
        functional='wB97XD',
        qm_input='pop=hirshfeld',
        program='gaussian'
    )    
    # Error from lack of method set, even with qm_method determined.
    with pytest.raises(UserWarning):
        model.set_region(selection=CSACutoff, name=None, min_atoms=None, max_atoms=None, memory=None, nprocs=None, method=None)

    # Error from catalytic center containing protein residues
    model = GenerateModel(PDB)
    model.set_catalytic_center('resid 10')
    qm_method = QMzyme.QM_Method(
        basis_set='6-31G*',
        functional='wB97XD',
        qm_input='pop=hirshfeld',
        program='gaussian'
    )  
    with pytest.raises(UserWarning):
        model.set_region(selection=CSACutoff, name=None, min_atoms=950, max_atoms=1050, memory=None, nprocs=None, method=qm_method)

    # Smaller system setup
    model = GenerateModel(Cutoff_3)
    model.set_catalytic_center('resid 263')
    qm_method = QMzyme.QM_Method(
        basis_set='6-31G*',
        functional='wB97XD',
        qm_input='pop=hirshfeld',
        program='gaussian'
    )  
    model.set_region(selection=CSACutoff, name=None, min_atoms=950, max_atoms=1050, memory=None, nprocs=None, method=qm_method)

    # Full system setup
    model = GenerateModel(PDB)
    model.set_catalytic_center('resid 263')
    qm_method = QMzyme.QM_Method(
        basis_set='6-31G*',
        functional='wB97XD',
        qm_input='pop=hirshfeld',
        program='gaussian'
    )  
    model.set_region(selection=CSACutoff, name=None, min_atoms=200, max_atoms=350, memory=None, nprocs=None, method=qm_method)

    # Error from using wrong step of CSACutoff
    with pytest.raises(UserWarning):
        model.set_region(selection=CSACutoff, name=None, method=None, holo_output_files=None, apo_output_files=None, pop="hirshfeld", charge_threshold=0.05)

    # Parameters for normal selection
    model = GenerateModel(PDB)
    model.set_catalytic_center('resid 263')
    model.set_region(selection=CSACutoff, name=None, min_atoms=950, max_atoms=1050, memory=None, nprocs=None, method=qm_method)
    assert len(model.regions) == 6
    assert model.catalytic_center.n_residues == 1
    assert model.catalytic_center.n_atoms == 37
    assert model.CSA_holo_truncated.n_residues == 88
    assert model.CSA_holo_truncated.n_atoms == 982
    assert model.CSA_apo_truncated.n_residues == 87
    assert model.CSA_apo_truncated.n_atoms == 945

    model = GenerateModel(PDB)
    model.set_catalytic_center('resid 263')

    # Parameters for D103A selection
    model.set_region(selection=CSACutoff, name=None, min_atoms=950, max_atoms=1050, memory='64GB', nprocs=200, method=qm_method, alanine_mutation="resid 40")
    assert len(model.regions) == 6
    assert model.catalytic_center.n_residues == 1
    assert model.catalytic_center.n_atoms == 37
    assert model.CSA_holo_truncated.n_residues == 88
    assert model.CSA_holo_truncated.n_atoms == 982
    assert model.CSA_apo_truncated.n_residues == 87
    assert model.CSA_apo_truncated.n_atoms == 943

    # Second part of the CSACutoff
    # With non-pickle file, error raised
    model = GenerateModel(PDB)
    model.set_catalytic_center('resid 263')

    with pytest.raises(UserWarning):
        model.set_region(selection=CSACutoff, name=None, method=None, holo_output_files=None, apo_output_files=None, pop="cm5", charge_threshold=0.05)

    with pytest.raises(UserWarning):
        model.set_region(selection=CSACutoff, name=None, method=qm_method, holo_output_files=CSA_holo, apo_output_files=CSA_apo, pop="cm5")

    with pytest.raises(UserWarning):
        model.set_region(selection=CSACutoff, method=qm_method, holo_output_files=CSA_holo, apo_output_files=CSA_apo, pop="cm5", charge_threshold=0.05)
    
    # With pickle file, charge threshold not stated
    with open(CSA_pkl, 'rb') as f:
        model = pickle.load(f)
    
    with pytest.raises(UserWarning):
        model.set_region(selection=CSACutoff, method=qm_method, holo_output_files=CSA_holo, apo_output_files=CSA_apo, pop="cm5")

    # Proper run of the second part
    model.set_region(selection=CSACutoff, method=qm_method, holo_output_files=CSA_holo, apo_output_files=CSA_apo, pop="cm5", charge_threshold=0.05)
    assert model.catalytic_center.n_residues == 1
    assert model.catalytic_center.n_atoms == 37
    assert model.CSA_cutoff_region.n_residues == 5
    assert model.CSA_cutoff_region.n_atoms == 77

    print(model.CSA_cutoff_region.creation_params)

    # Truncation check
    model.truncate()
    with pytest.raises(UserWarning):
        model.truncate()
    
    model = GenerateModel(PDB)
    model.set_catalytic_center('resid 263')
    with pytest.raises(UserWarning):
        model.truncate()

    restore_directory()
    assert 'QCALC' not in os.listdir()