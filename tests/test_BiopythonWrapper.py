"""
Tests for the QMzyme BiopythonWrapper.py module.
"""


import os
import shutil
import QMzyme
import pytest
from QMzyme.BiopythonWrapper import *
from importlib_resources import files

pdb_file = str(files('QMzyme.data').joinpath('1oh0_equ_from_amber_sim.pdb'))
original_files = os.listdir()

def restore_folder():
    for file in os.listdir():
        if file not in original_files:
            try:
                os.remove(file)
            except:
                shutil.rmtree(file)

@pytest.mark.parametrize(
        'test_type, pdb_file, restore',
        [
            ("Result after within_distance()", pdb_file, True),
            ("Result after truncate()", pdb_file, True),
        ]
)
def test_write_pdb(test_type, pdb_file, restore):
    model = QMzyme.GenerateModel(pdb_file)
    model.set_catalytic_center(resnumber=263)
    model.within_distance(5)

    if test_type == 'Result after within_distance()':
        write_pdb(model.models[-1], 'test.pdb')
        line0 = ['ATOM', '233', 'N', 'TYR', 'A', '16', '47.082', '33.457', '31.192', '1.00', '0.00', 'N']
        line100 = ['ATOM', '867', 'HD2', 'TYR', 'A', '57', '42.459', '42.329', '25.761', '1.00', '0.00', 'H']
    
    elif test_type == 'Result after truncate()':
        model.truncate()
        write_pdb(model.models[-1], 'test.pdb')
        line0 = ['ATOM', '233', 'Hcap', 'TYR', 'A', '16', '46.850', '33.818', '31.003', '1.00', '0.00', 'H']
        line100 = ['ATOM', '900', 'HB3', 'GLN', 'A', '59', '39.244', '47.428', '21.643', '1.00', '0.00', 'H']

    with open('test.pdb', 'r') as f:
        data = f.readlines()

    assert data[0].split() == line0
    assert data[100].split() == line100

    if restore is True:
        restore_folder()

def test_count_atoms():
    model = QMzyme.GenerateModel(pdb_file)
    assert count_atoms(model) == 4258

def test_count_residues():
    model = QMzyme.GenerateModel(pdb_file)
    assert count_residues(model) == 324

def test_count_nonprotein_residues():
    model = QMzyme.GenerateModel(pdb_file)
    assert count_nonprotein_residues(model) == 62
    
def test_count_chains():
    model = QMzyme.GenerateModel(pdb_file)

    assert count_chains(model) == 2

def test_remove_atom():
    model = QMzyme.GenerateModel(pdb_file)

    #remove all water atoms
    for atom in model.get_atoms():
        if atom.parent.resname == 'WAT':
            remove_atom(atom)
    
    assert count_atoms(model) == 4138

def test_get_element_idx():
    model = QMzyme.GenerateModel(pdb_file)
    exp = [12, 191, 475, 1036, 1201, 1232, 1338, 1444, 1563, 1750, 2231, 2410, 2694, 3255, 3420, 3451, 3557, 3663, 3782, 3969]

    assert get_element_idx(model, 'S') == exp

def test_get_atom_idx():
    model = QMzyme.GenerateModel(pdb_file)
    exp1 = [12, 191, 475, 1232, 1338, 1563, 1750, 2231, 2410, 2694, 3451, 3557, 3782, 3969]
    exp2 = [2003, 4222]

    assert get_atom_idx(model, 'SD') == exp1
    assert get_atom_idx(model, 'O1') == exp2

def test_make_model_dict():
    model = QMzyme.GenerateModel(pdb_file)
    model.set_catalytic_center(resnumber=263)
    model.within_distance(5)
    model.truncate()
    exp_keys = ['Residues', '_id', 'full_id', 'method']

    assert list(make_model_dict(model.models[-1]).keys()) == exp_keys



    





    





