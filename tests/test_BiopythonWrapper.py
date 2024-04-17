"""
Tests for the QMzyme BiopythonWrapper.py module.
"""


import sys
import QMzyme
from QMzyme.BiopythonWrapper import *
from importlib_resources import files

pdb_file = str(files('QMzyme.data').joinpath('1oh0_equ_xstal.pdb'))


def test_write_pdb(pdb_file):
    model = QMzyme.GenerateModel(pdb_file)




