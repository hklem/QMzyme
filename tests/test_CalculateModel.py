from QMzyme.CalculateModel import QM_Method
from QMzyme.GenerateModel import GenerateModel
#from importlib_resources import files
from QMzyme.data import PDB
import os
import shutil
import pytest

original_contents = os.listdir()

def restore_directory():
    for name in os.listdir():
        if name not in original_contents:
            try:
                os.remove(name)
            except:
                shutil.rmtree(name)

@pytest.mark.parametrize(
    "Test, program",[
        ('Gaussian Test', 'gaussian'),
        ('Orca Test', 'orca'),
    ]
)
def test_QM_Method(Test, program):
    model = GenerateModel(PDB)
    region_name = 'cutoff_3'
    model.set_region(name=region_name, selection='byres around 3 resid 263')
    # check basis_set info is set for all atoms
    bs1 = '6-31g(d)'
    ids = model.cutoff_3.get_ids('name', 'CA')
    model.cutoff_3.set_fixed_atoms(ids)
    region, basis_set, functional = model.cutoff_3, '6-31G*', 'wB97X-D3'
    qm_method = QM_Method(
        basis_set=basis_set, 
        functional=functional, 
        qm_input='OPT FREQ', 
        program=program
    )
    region.set_charge(-100)
    qm_method.assign_to_region(region=region)
    
    # check region method:
    assert hasattr(region, "method")
    assert region.method["basis_set"] == basis_set
    assert region.method["functional"] == functional
    assert region.method["charge"] == -100
    assert region.method["mult"] == 1
    assert region.method["freeze_atoms"] == region.get_ix_array_from_ids(ids)

    model.write_input()
    restore_directory()
