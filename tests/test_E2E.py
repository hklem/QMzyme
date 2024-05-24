"""
End to end tests designed to match up with documentation examples.
"""
from QMzyme.GenerateModel import GenerateModel
from QMzyme.SelectionSchemes import *
from QMzyme.TruncationSchemes import *
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


def test_QM_only_calculation():
    from QMzyme.CalculateModel import QM_Method
    # Initialize model and define regions
    model = GenerateModel(PDB)
    model.set_catalytic_center(selection='resid 263')
    model.set_region(selection=DistanceCutoff, cutoff=5)
    assert model.n_regions == 2
    assert len(model.regions) == 2
    assert 'catalytic_center' in model.get_region_names()
    assert 'cutoff_5' in model.get_region_names()

    # Truncate region
    model.truncate_region(region=model.cutoff_5, scheme=CA_terminal)
    assert model.n_regions == 3
    assert 'cutoff_5_truncated' in model.get_region_names()

    # Designate coordinate constraints 
    c_alpha_atoms = model.cutoff_5_truncated.get_atoms(attribute='name', value='CA')
    model.cutoff_5_truncated.set_fixed_atoms(atoms=c_alpha_atoms)
    for atom in c_alpha_atoms:
        assert atom.is_fixed == True
    for atom in model.cutoff_5.atoms:
        assert atom.is_fixed == False

    # Build the QM Method and Assign to Region
    qm_method = QM_Method(
        basis_set='6-31G*', 
        functional='wB97X-D3', 
        qm_input='OPT FREQ', 
        program='orca'
    )

    qm_method.assign_to_region(region=model.cutoff_5_truncated)
    assert model.cutoff_5.method == None
    assert model.cutoff_5_truncated.method["qm_input"] == '6-31G* wB97X-D3 OPT FREQ'
    assert model.cutoff_5_truncated.charge == -1

    # Write QM Input File
    model.write_input()
    assert 'QCALC' in os.listdir()
    assert 'cutoff_5_truncated.inp' in os.listdir('QCALC')
    with open(os.path.join('QCALC', 'cutoff_5_truncated.inp'), 'r') as f:
        file = f.read()
    qm_input = file[file.find('!')+1:].split('\n')[0]
    assert qm_input.endswith(model.cutoff_5_truncated.method["qm_input"])
    assert f'%geom Constraints' in file
    assert f'%QMMM' not in file
    restore_directory()
    assert 'QCALC' not in os.listdir()

@pytest.mark.parametrize(
    "qm1_dict, qm2_dict",[
        ({'basis_set': '6-31+G**',
          "functional": 'wB97X-D3',
          "qm_input": 'OPT FREQ',
          "program": 'orca'},
          {'basis_set': '6-31G*',
          "functional": 'wB97X-D3'}),
    ]
)
def test_QMQM2_calculation(qm1_dict, qm2_dict):
    from QMzyme.CalculateModel import QM_Method
    model = GenerateModel(PDB)
    model.set_catalytic_center(selection='resid 263')
    model.catalytic_center.set_charge(-1)
    assert model.catalytic_center.charge == -1
    model.set_region(selection=DistanceCutoff, cutoff=5)
    model.truncate_region(region=model.cutoff_5)
    c_alpha_atoms = model.cutoff_5_truncated.get_atoms(attribute='name', value='CA')
    model.cutoff_5_truncated.set_fixed_atoms(atoms=c_alpha_atoms)

    qm1_method = QM_Method(**qm1_dict)
    qm2_method = QM_Method(**qm2_dict)
    qm1_method.assign_to_region(region=model.catalytic_center)
    qm2_method.assign_to_region(region=model.cutoff_5_truncated)
    for attr, val in qm1_dict.items():
        if attr == 'qm_input':
            continue
        assert model.catalytic_center.method[attr] == val
    for attr, val in qm2_dict.items():
        if attr == 'qm_input':
            continue
        assert model.cutoff_5_truncated.method[attr] == val
    
    model.write_input("temp")
    assert 'temp.inp' in os.listdir('QCALC')
    with open(os.path.join('QCALC', 'temp.inp'), 'r') as f:
        file = f.read()
    qm_input = file[file.find('!')+1:].split('\n')[0]
    assert qm_input.strip().startswith('QM/QM2')
    assert model.catalytic_center.method['qm_input'] in file
    assert f'%QMMM' in file
    qm2_input = file[file.find('%QMMM'):].split('\n')[0]
    val = model.cutoff_5_truncated.method["qm_input"]
    assert qm2_input.strip().endswith(f"QM2CUSTOMMETHOD '{val}'")
    assert "QMATOMS "+"{"+"0:37"+"}" in file
    assert len(file.split("%QMMM")) == 2

    # check that you cannot set a third QM method, or call gaussian with QMQM2
    model.set_region(selection=DistanceCutoff, cutoff=6)
    with pytest.raises(UserWarning):
        model.set_region(selection=DistanceCutoff, cutoff=6)
        qm2_method.assign_to_region(region=model.cutoff_6)
        model.catalytic_center.method["program"] = 'gaussian'

    model.write_input()

    restore_directory()

def test_QMXTB_calculation():
    from QMzyme.CalculateModel import CalculateModel, QM_Method, XTB_Method
    model = GenerateModel(PDB)
    model.set_catalytic_center(selection='resid 263')
    model.catalytic_center.set_charge(-1)
    model.set_region(selection=DistanceCutoff, cutoff=5)
    model.truncate_region(region=model.cutoff_5)
    c_alpha_atoms = model.cutoff_5_truncated.get_atoms(attribute='name', value='CA')
    model.cutoff_5_truncated.set_fixed_atoms(atoms=c_alpha_atoms)

    # check that warning is raised if a non-QM method is set before a QM method is set.
    with pytest.raises(UserWarning):
        assert XTB_Method().assign_to_region(region=model.cutoff_5_truncated)

    qm_method = QM_Method(basis_set='6-31G*', 
               functional='wB97X-D3', 
               qm_input='OPT FREQ', 
               program='orca')

    qm_method.assign_to_region(region=model.catalytic_center)
    XTB_Method().assign_to_region(region=model.cutoff_5_truncated)
    model.write_input("temp")
    assert CalculateModel.calculation['XTB'].n_atoms != model.cutoff_5_truncated.n_atoms
    with open(os.path.join('QCALC', 'temp.inp'), 'r') as f:
        file = f.read()
    assert "%QMMM\n" in file
    qm_input = file[file.find('!')+1:].split('\n')[0]
    assert qm_input.strip().startswith('QM/XTB')
    model.write_input("temp2")
    with open(os.path.join('QCALC', 'temp2.inp'), 'r') as f:
        file = f.read()
    assert len(file.split('%QMMM')) == 2
    assert len(file.split('QM/XTB')) == 2

    restore_directory()








