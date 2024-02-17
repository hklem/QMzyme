"""
Tests for the QMzyme generate.py code.
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

def test_init(init_file=amber_file):
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
        model.catalytic_center(sel=sel, 
                               res_name=res_name, 
                               res_number=res_number, 
                               chain=chain, 
                               save_file=False)
        line_first = 'ATOM      1  C1  EQU A 263      43.927  42.973  31.198  1.00  0.00           C'
        line_7 = 'ATOM      7  C6  EQU A 263      44.793  42.905  30.106  1.00  0.00           C'
        line_last = 'END'
        assert model.catalytic_center_pdb[0].strip() == line_first
        assert model.catalytic_center_pdb[6].strip() == line_7
        assert model.catalytic_center_pdb[-1].strip() == line_last

@pytest.mark.parametrize(
        'test_type, distance_cutoff, init_file',
        [
            ('standard', 5, amber_file)
        ]
)
def test_subsystem(test_type, distance_cutoff, init_file):
        model.subsystem(distance_cutoff=distance_cutoff, 
                        save_file=False)
        line_first = 'ATOM      1  N   TYR A  16      47.082  33.457  31.192  1.00  0.00           N'
        line_7 = 'ATOM      7  HB3 TYR A  16      45.513  34.980  32.548  1.00  0.00           H'
        line_last = 'END'
        assert model.subsystem_pdb[0].strip() == line_first
        assert model.subsystem_pdb[0].strip() == line_7
        assert model.subsystem_pdb[0].strip() == line_last

@pytest.mark.parametrize(
        'test_type, truncation_scheme, constrain_atoms, init_file',
        [
            ('default', 'CA_terminal', ['CA'], amber_file),
            ('constrain_CA_and_CB', 'CA_terminal', ['CA', 'CB'], amber_file)
        ]
)
def test_truncate(test_type, truncation_scheme, constrain_atoms, init_file):
        model.truncate(scheme=truncation_scheme, 
                       constrain_atoms=constrain_atoms, 
                       save_file=False)
        
        if test_type == 'default':
                constrain_list = [2,23,40,54,74,78,97,114,121,138,152,160,175,193,207,222,239,253,264,279,287,309]
        if test_type == 'constrain_CA_and_CB':
                constrain_list = [2,4,23,25,40,42,54,56,71,74,78,80,97,99,114,121,123,138,140,152,154,160,162,175,177,193,195,207,209,222,224,239,241,253,255,264,266,279,281,287,289,309,311]

        line_first = 'ATOM      1  H*  TYR A  16      46.850  33.818  31.003  1.00  0.00           H'
        line_7 = 'ATOM      7  CG  TYR A  16      44.523  36.280  31.198  1.00  0.00           C'
        n_atoms = 391
        charge = -1

        assert model.truncated_subsystem_pdb[0].strip() == line_first
        assert model.truncated_subsystem_pdb[6].strip() == line_7
        assert model.model_atom_count == n_atoms
        assert model.subsystem_charge == charge
        assert model.constrain_atom_list == constrain_list
