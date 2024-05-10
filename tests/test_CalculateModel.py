from QMzyme.CalculateModel import QM_Region
from QMzyme.GenerateModel import GenerateModel
from importlib_resources import files
import os
import shutil
import pytest

pdb_file = str(files('QMzyme.data').joinpath('1oh0.pdb'))
model = GenerateModel(pdb_file)
model.set_region('cutoff_3', 'byres around 3 resid 263')

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
def test_QM_Region(Test, model, program):
    # check basis_set info is set for all atoms
    bs1 = '6-31g(d)'
    ids = model.cutoff_3.get_ids('name', 'CA')
    model.cutoff_3.set_fixed_atoms(ids)
    r1 = QM_Region(model.cutoff_3)
    r1.set_basis_set(bs1)
    for bs in [atom.basis_set for atom in r1.atoms]:
        assert bs == bs1

    # check functional info is set for all atoms
    functional = 'wb97xd'
    r1.set_functional(functional)
    for fnc in [atom.functional for atom in r1.atoms]:
        assert fnc == functional

    r1.set_charge(-1)
    r1.set_multiplicity(1)

    calc_info = r1.write_qm_input(program=program)

    assert calc_info['freeze_atoms'] == model.cutoff_3.get_idxs_from_ids(ids)

    #restore_directory()

