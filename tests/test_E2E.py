"""
End to end tests designed to match up with documentation examples.
"""
from QMzyme.GenerateModel import GenerateModel
from QMzyme.CalculateModel import CalculateModel
from QMzyme.SelectionSchemes import *
from QMzyme.TruncationSchemes import *
from QMzyme.data import PDB, PQR, protein_residues, residue_charges
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
    residue_charges.update({'EQU':-1})
    model = GenerateModel(PDB)
    model.set_catalytic_center(selection='resid 263')
    model.set_region(selection=DistanceCutoff, cutoff=5)
    assert model.n_regions == 2
    assert len(model.regions) == 2
    assert 'catalytic_center' in model.get_region_names()
    assert 'cutoff_5' in model.get_region_names()

    # Designate coordinate constraints 
    c_alpha_atoms = model.cutoff_5.get_atoms(attribute='name', value='CA')
    for atom in model.cutoff_5.atoms:
        assert atom.is_fixed == False
    model.cutoff_5.set_fixed_atoms(atoms=c_alpha_atoms)
    for atom in c_alpha_atoms:
        assert atom.is_fixed == True

    # Build the QM Method and Assign to Region
    qm_method = QM_Method(
        basis_set='6-31G*', 
        functional='wB97X-D3', 
        qm_input='OPT FREQ', 
        program='orca'
    )

    qm_method.assign_to_region(region=model.cutoff_5)
    assert model.cutoff_5.method["qm_input"] == '6-31G* wB97X-D3 OPT FREQ'
    assert model.cutoff_5.charge == -2
    for atom in c_alpha_atoms:
        assert CalculateModel.calculation['QM'].get_atom(id=atom.id).is_fixed == True

    # Truncate Model
    print(CalculateModel.calculation)
    print(CalculateModel.calc_type)
    model.truncate()

    # Write QM Input File
    model.write_input()
    assert 'QCALC' in os.listdir()
    assert 'cutoff_5_truncated.inp' in os.listdir('QCALC')
    with open(os.path.join('QCALC', 'cutoff_5_truncated.inp'), 'r') as f:
        file = f.read()
    qm_input = file[file.find('!')+1:].split('\n')[0]
    assert qm_input.endswith(model.cutoff_5.method["qm_input"])
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
    residue_charges.update({'EQU':-1})
    model = GenerateModel(PDB)
    model.set_catalytic_center(selection='resid 263')
    model.catalytic_center.set_charge(-1)
    assert model.catalytic_center.charge == -1
    model.set_region(selection=DistanceCutoff, cutoff=5)
    c_alpha_atoms = model.cutoff_5.get_atoms(attribute='name', value='CA')
    model.cutoff_5.set_fixed_atoms(atoms=c_alpha_atoms)

    qm1_method = QM_Method(**qm1_dict)
    qm2_method = QM_Method(**qm2_dict)
    print("Before any method assigning")
    qm1_method.assign_to_region(region=model.catalytic_center)
    print("After first method calc type: ",CalculateModel.calc_type)    
    qm2_method.assign_to_region(region=model.cutoff_5)
    print("After second method calc type: ",CalculateModel.calc_type)
    for attr, val in qm1_dict.items():
        if attr == 'qm_input':
            continue
        assert model.catalytic_center.method[attr] == val
    for attr, val in qm2_dict.items():
        if attr == 'qm_input':
            continue
        assert model.cutoff_5.method[attr] == val

    # Truncate Model
    #model.truncate()
    model.write_input("qmqm2")
    assert 'qmqm2.inp' in os.listdir('QCALC')
    with open(os.path.join('QCALC', 'qmqm2.inp'), 'r') as f:
        file = f.readlines()
    assert 'QM/QM2' in file[3]
    for component in model.catalytic_center.method['qm_input'].split('%QMMM')[0].split():
        assert component in file[3]

    assert '%QMMM' in file[4]
    assert file[4].strip() == "%QMMM QM2CUSTOMMETHOD '6-31G* wB97X-D3'"
    assert file[5].strip() == "QMATOMS {"+"360:396} END"
    assert file[6].strip() == "Charge_Total -2 END"
    assert file[7].strip() == f"%geom Constraints"

    restore_directory()


def test_QMXTB_calculation():
    from QMzyme.CalculateModel import CalculateModel, QM_Method, XTB_Method
    residue_charges.update({'EQU':-1})
    model = GenerateModel(PDB)
    model.set_catalytic_center(selection='resid 263')
    model.set_region(selection=DistanceCutoff, cutoff=5)
    c_alpha_atoms = model.cutoff_5.get_atoms(attribute='name', value='CA')
    model.cutoff_5.set_fixed_atoms(atoms=c_alpha_atoms)

    # check that warning is raised if a non-QM method is set before a QM method is set.
    # with pytest.raises(UserWarning):
    #     assert XTB_Method().assign_to_region(region=model.cutoff_5)

    qm_method = QM_Method(basis_set='6-31G*', 
               functional='wB97X-D3', 
               qm_input='OPT FREQ', 
               program='orca')

    qm_method.assign_to_region(region=model.catalytic_center)
    XTB_Method().assign_to_region(region=model.cutoff_5)
    model.truncate()
    model.write_input("test_qmxtb")
    assert CalculateModel.calculation['QMXTB'].n_atoms != model.cutoff_5.n_atoms
    with open(os.path.join('QCALC', 'test_qmxtb.inp'), 'r') as f:
        file = f.readlines()

    assert file[3].strip() == "! QM/XTB 6-31G* wB97X-D3 OPT FREQ"
    assert file[4].strip() == f"%QMMM"
    assert file[5].strip() == "QMATOMS {"+"324:360} END"
    assert file[6].strip() == "Charge_Total -2 END"
    assert file[32].strip() == "* xyz -1 1"

    # qm_input = file[file.find('!')+1:].split('\n')[0]
    # assert qm_input.strip().startswith('QM/XTB')

    # write again to make sure the method gets reset 
    model.write_input("test_qmxtb")
    with open(os.path.join('QCALC', 'test_qmxtb.inp'), 'r') as f:
        file = f.read()

    assert len(file.split('%QMMM')) == 2
    assert len(file.split('QM/XTB')) == 2

    restore_directory()

def test2_QM_XTB_calculation():
    from QMzyme.CalculateModel import QM_Method, XTB_Method
    qm_method = QM_Method(
        basis_set='6-31G*', 
        functional='wB97X-D3', 
        qm_input='OPT FREQ', 
        program='orca'
    )

    # Initialize model and define regions
    residue_charges.update({'EQU':-1})
    model = GenerateModel(PQR)
    #model = GenerateModel(PDB)

    model.set_catalytic_center(selection='resid 263 or resid 16')
    qm_method.assign_to_region(region=model.catalytic_center)
    model.set_region(name='xtb_region', selection=DistanceCutoff, cutoff=3)
    model.set_region(name='ile17', selection='resid 17')
    xtb_region = model.xtb_region.combine(model.ile17)
    c_alpha_atoms = xtb_region.get_atoms(attribute='name', value='CA')
    xtb_region.set_fixed_atoms(atoms=c_alpha_atoms)
    xtb_method = XTB_Method().assign_to_region(region=xtb_region)
    with pytest.raises(UserWarning):
        # you cannot add a region with the same name
        model.set_region(name='xtb_region', selection=xtb_region)
    model.remove_region('xtb_region') # replacing this region
    # need to add to model so combined regions can be truncated 
    model.set_region(name='xtb_region', selection=xtb_region)

    model.truncate()
    model.write_input('test2_qmxtb')

    with open(os.path.join('QCALC', 'test2_qmxtb.inp'), 'r') as f:
        file = f.readlines()

    assert file[3].strip() == "! QM/XTB 6-31G* wB97X-D3 OPT FREQ"
    assert file[4].strip() == f"%QMMM"
    assert file[5].strip() == "QMATOMS {"+"57:77} {413:449"+"} END"
    assert file[6].strip() == "Charge_Total -1 END"

    # check writing of pymol visualization script
    model.pymol_visualize()
    assert f'QMzymeModel_{model.name}_visualize.py' in os.listdir()
