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
        model = QMzyme.GenerateModel(structure=init_file)

@pytest.mark.parametrize(
        'test_type, resname, resnumber, chain, init_file, overwrite',
        [
            ("default", None, 103, None, amber_file, True),
            ("overwrite", None, 263, None, amber_file, False),
            ("default", None, 263, None, amber_file, True),
            
        ]
)
def test_catalytic_center(test_type, resname, resnumber, chain, init_file, overwrite):
        model.set_catalytic_center(resnumber=resnumber, overwrite=overwrite)
        if overwrite is True:
                if resnumber == 263:
                        assert len(model.catalytic_center) == 1
                        assert len(model.catalytic_center[0].list_atoms()) == 37
                elif resnumber == 103:
                        assert len(model.catalytic_center) == 1
                        assert len(model.catalytic_center[0].list_atoms()) == 13
        elif overwrite is False:
                assert len(model.catalytic_center) == 2
                l1 = len(model.catalytic_center[0].list_atoms())
                l2 = len(model.catalytic_center[1].list_atoms())
                assert l1 + l2 == 50
 
@pytest.mark.parametrize(
        'test_type, distance_cutoff, init_file',
        [
            ('default', 5, amber_file)
        ]
)
def test_within_distance(test_type, distance_cutoff, init_file):
        model.within_distance(distance_cutoff)
        
        assert len(model.models[-1].list_atoms()) == 427
        assert len(model.models[-1].list_residues()) == 33
        assert model.models[-1].method['type'] == 'within_distance'
        assert model.models[-1].method['cutoff'] == distance_cutoff
        assert model.models[-1].method['catalytic_center'] == model.catalytic_center




@pytest.mark.parametrize(
        'test_type, truncation_scheme, init_file',
        [
            ('default', 'CA_terminal', amber_file),
        ]
)
def test_truncate(test_type, truncation_scheme, init_file):
        model.truncate(scheme=truncation_scheme)
        
        assert len(model.models[-1].list_atoms()) == 391
        assert len(model.models[-1].list_residues()) == 33

        atom0 = model.models[-1].list_atoms()[0]
        assert atom0.element == 'H'
        assert atom0.name == 'Hcap'
        assert atom0.coord[0]== 46.84952163696289
        assert atom0.coord[1]== 33.818172454833984
        assert atom0.coord[2]== 31.003271102905273

        atom1 = model.models[-1].list_atoms()[1]
        assert atom1.name == 'CA'

def test_calculateQM(test_type, init_file)
