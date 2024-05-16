from QMzyme.CalculateModel import QM, write_QM
from QMzyme.GenerateModel import GenerateModel
from importlib_resources import files
import os
import shutil
import pytest

pdb_file = str(files('QMzyme.data').joinpath('1oh0.pdb'))
model = GenerateModel(pdb_file)
region_name = 'cutoff_3'
model.set_region(region_name=region_name, selection='byres around 3 resid 263', layer=None)

original_contents = os.listdir()

def restore_directory():
    for name in os.listdir():
        if name not in original_contents:
            try:
                os.remove(name)
            except:
                shutil.rmtree(name)

@pytest.mark.parametrize(
    "Test, model, program",[
        ('Gaussian Test', model, 'gaussian'),
        ('Orca Test', model, 'orca'),
    ]
)
def test_QM(Test, model, program):
    # check basis_set info is set for all atoms
    bs1 = '6-31g(d)'
    ids = model.cutoff_3.get_ids('name', 'CA')
    model.cutoff_3.set_fixed_atoms(ids)
    region, basis_set, functional, charge, mult = model.cutoff_3, '6-31g(d)', 'wb97xd', -1, 1
    r1 = QM(region, basis_set, functional, charge, mult)
    
    # check region method:
    print(r1.__dict__)
    assert hasattr(region, "method")
    assert region.method["basis_set"] == basis_set
    assert region.method["functional"] == functional
    assert region.method["charge"] == charge
    assert region.method["mult"] == mult
    assert region.method["freeze_atoms"] == region.get_idxs_from_ids(ids)

    write_QM(region)
    assert 'QCALC' in os.listdir()
    restore_directory()
    assert 'QCALC' not in os.listdir()

