"""
Unit and regression test for the QMzyme package.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow

import sys
import pytest
import QMzyme
from rdkit import Chem
import os

path = os.path.join('etc')

def grab_test_data(test_file):
    with try open(test_file, 'r') as f:
        data=f.readlines()
        raise Error("Test file {test_file} not found.")
    residue_sequence 
    number_of_residues
    atom_sequence
    number_of_atoms
    atom_coords
    resid_order

def test_QMzyme_imported():
    """
    Sample test, will always pass so long as import statement worked.
    """
    assert "QMzyme" in sys.modules

def test_catalytic_center(test_file):
    """
    tbd
    """
    test_file = os.path.join(path,test_file)
    with try open(test_file,'r'):
        print('It worked')
        raise Error('Test file {test_file} was not found.')

def test_active_site(test_file='test_catalytic_center_chainA_DNX_202_active_site_distance_cutoff_4.pdb'):
    """
    tbd
    """

def test_truncation(test_file='test_catalytic_center_chainA_DNX_202_truncated_active_site_distance_cutoff_4.pdb'):
    """
    tbd
    """

