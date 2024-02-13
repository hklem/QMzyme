"""
Tests for the QMzyme package.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow

import sys
import os
import pytest
import QMzyme
from importlib_resources import files, as_file

model = None
amber_file = str(files('QMzyme.data').joinpath('1oh0_equ_from_amber_sim.pdb'))

def test_generate_model(init_file=amber_file):
        assert "QMzyme" in sys.modules
        global model
        model = QMzyme.GenerateModel(protein_file=init_file, 
                                     save_json=False,
                                     verbose=False)

@pytest.mark.parametrize(
        'test_type, sel, res_name, res_number, chain, init_file',
        [
            ("mdanalysis_resnum", "resid 263", None, None, None, amber_file),
            ("mdanalysis_resname_and_segid", "resname EQU and segid A", None, None, None, amber_file),
            ("rdkit_selection", None, "EQU", None, "A", amber_file)
        ]
)
def test_catalytic_center(test_type, sel, res_name, res_number, chain, init_file):
        model.catalytic_center(sel=sel, res_name=res_name, res_number=res_number, chain=chain, save_file=False)
        line_first = 'ATOM      1  C1  EQU A 263      43.927  42.973  31.198  1.00  0.00           C'
        line_7 = 'ATOM      7  C6  EQU A 263      44.793  42.905  30.106  1.00  0.00           C'
        line_last = 'END'
        assert model.catalytic_center_pdb[0].strip() == line_first
        assert model.catalytic_center_pdb[6].strip() == line_7
        assert model.catalytic_center_pdb[-1].strip() == line_last


### Might use parameterized tests later on when I expand the scope of tests, but not there yet.
#@pytest.mark.parametrize(
#        'test_function, test_type, init_file, target',
#        [
#            ("init", None, amber_file, None),
#            ("catalytic_center", "selection_string1", amber_file, None),
#            ("catalytic_center", "selection_string2", amber_file, None),
#            ("catalytic_center", "selection_notstring", amber_file, None),
#            ("subsystem", None, amber_file, None),
#            ("truncate", "CA_terminal", amber_file, None),
#            (None, None, None, None)
#        ]
#    )

#def test_QMzyme_generate(test_function, test_type, init_file, target):
#    if test_function == "init":
#        assert "QMzyme" in sys.modules
#        global model
#        model = QMzyme.GenerateModel(protein_file=init_file, 
#                                     save_json=False,
#                                     verbose=False)

#    if test_function == "catalytic_center":
#        if test_type == "selection_string1":
#            model.catalytic_center(sel='resid 263',
#                                   save_file=False)
#        if test_type == "selection_string2":
#            model.catalytic_center(sel='resname EQU and segid A',
#                                   save_file=False)
#        if test_type == "selection_notstring":
#            model.catalytic_center(res_name='EQU',
#                                   chain='A',
#                                   save_file=False)
            
#        line_first = 'ATOM      1  C1  EQU A 263      43.927  42.973  31.198  1.00  0.00           C'
#        line_7 = 'ATOM      7  C6  EQU A 263      44.793  42.905  30.106  1.00  0.00           C'
#        line_last = 'END'
#        assert model.catalytic_center_pdb[0].strip() == line_first
#        assert model.catalytic_center_pdb[6].strip() == line_7
#        assert model.catalytic_center_pdb[-1].strip() == line_last

#    if test_function == "subsystem":
#        model.subsystem(distance_cutoff=4, save_file=False)








