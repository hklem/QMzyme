"""
Tests for the QMzyme utils.py code.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow

import sys
import os
import pytest
from QMzyme.utils import *
from importlib_resources import files, as_file

model = None
amber_file = str(files('QMzyme.data').joinpath('1oh0_equ_from_amber_sim.pdb'))
lig_file = str(files('QMzyme.data').joinpath('small_hetatm_ligand.pdb'))

#def test_download_pdb(pdb_code):
#        download_pdb('7ac8')

