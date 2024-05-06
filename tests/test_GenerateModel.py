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
        qmz_region = RegionBuilder(region_name, mda_atomgroup).get_region()
        model.set_region(region_name, qmz_region)
    assert len(model.regions) == 1
    assert model.regions[0].name == region_name
    assert model.regions[0].n_atoms == 40
    assert model.regions[0].n_residues == 2

@pytest.mark.parametrize(
    "Test, init_file, region_selection",[
        #('First and last residue in protein: MET1 GLN262', pdb_file, 'resid 1 or resid 262'),
        #('MET1 ASN2', pdb_file, 'resid 1 or resid 2'),
        #('MET1 LEU3', pdb_file, 'resid 1 or resid 3'),
        ('ASN2 THR5', pdb_file, 'resid 2 or resid 5'),
        ('ASN2 LEU3 THR5 ALA6', pdb_file, 'resid 2 or resid 3 or resid 5 or resid 6'),
        ('PRO4', pdb_file, 'resid 4'), 
        #('PRO4 THR5', pdb_file, 'resid 4 or resid 5'),
        #('LEU3 PRO4', pdb_file, 'resid 3 or resid 4'),
        #('Non protein residue: WAT265', pdb_file, 'resid 265'),
    ]
)
def test_truncate_region_CA_terminal(Test, init_file, region_selection, truncation_scheme="CA_terminal"):
    model = GenerateModel(file=init_file)
    original_region = model.set_region('test', region_selection)
    truncated_region = model.truncate_region(original_region, truncation_scheme)
    #First check that the original region didn't change:
    original_first_res = original_region.residues[0]
    truncated_first_res = truncated_region.residues[0]
    original_last_res = original_region.residues[-1]
    truncated_last_res = truncated_region.residues[-1]

    if original_first_res.resname != 'PRO':
        assert 'H' in [atom.name for atom in original_first_res.atoms]
        assert 'H' not in [atom.name for atom in truncated_first_res.atoms]
    
    if original_first_res.resname == 'PRO':
        assert 'N' in [atom.name for atom in original_first_res.atoms]
        assert 'H' not in [atom.name for atom in original_first_res.atoms]
        assert 'N' in [atom.name for atom in truncated_first_res.atoms]

    assert 'H1' in [atom.name for atom in truncated_first_res.atoms]
    assert 'C' in [atom.name for atom in original_last_res.atoms]
    assert 'C' not in [atom.name for atom in truncated_last_res.atoms]
    assert 'O' in [atom.name for atom in original_last_res.atoms]
    assert 'O' not in [atom.name for atom in truncated_last_res.atoms]


