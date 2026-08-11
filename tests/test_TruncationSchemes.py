"""
Tests for the QMzyme truncation_schemes.py and truncation_utils.py codes.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow

import QMzyme
from QMzyme.GenerateModel import GenerateModel
from QMzyme.TruncationSchemes import *
import pytest
from QMzyme.data import PDB
import os
import shutil

original_contents = os.listdir()

def restore_directory():
    for name in os.listdir():
        if name not in original_contents:
            try:
                os.remove(name)
            except:
                shutil.rmtree(name)

@pytest.mark.parametrize(
    "Test, region_selection, trunc_scheme, isolated_gly_ala",[
        ('MET1 ASN2', 'resid 1 or resid 2', TerminalAlphaCarbon, False),
        ('THR5 ALA6', 'resid 5 or resid 6', TerminalAlphaCarbon, False),
        ('PRO4 ALA6', 'resid 4 or resid 6', TerminalAlphaCarbon, True),
        ('GLN10 GLY11', 'resid 10 or resid 11', TerminalAlphaCarbon, False),
        ('VAL9 GLY11', 'resid 9 or resid 11', TerminalAlphaCarbon, True),
        ('MET1 ASN2', 'resid 1 or resid 2', AlphaCarbon, False),
        ('THR5 ALA6', 'resid 5 or resid 6', AlphaCarbon, True),
        ('PRO4 ALA6', 'resid 4 or resid 6', AlphaCarbon, True),
        ('GLN10 GLY11', 'resid 10 or resid 11', AlphaCarbon, True),
        ('VAL9 GLY11', 'resid 9 or resid 11', AlphaCarbon, True),
    ]
)
def test_check_gly_ala(Test, region_selection, trunc_scheme, isolated_gly_ala):
    model = GenerateModel(PDB)
    model.set_region(name='region', selection=region_selection)
    qm_method = QMzyme.QM_Method(
        basis_set='6-31G*', 
        functional='wB97X-D3', 
        qm_input='OPT FREQ', 
        program='orca'
    )

    qm_method.assign_to_region(region=model.region)
    if isolated_gly_ala == True:
        with pytest.raises(ValueError):
            model.truncate(scheme=trunc_scheme)
    if isolated_gly_ala == False:
        model.truncate(scheme=trunc_scheme)

@pytest.mark.parametrize(
    "Test, region_selection, trunc_scheme, remove_flag_name, remove_flag_value",[
        ('PRO4 ALA6, remove_ethane=True removes isolated ALA', 'resid 4 or resid 6', TerminalAlphaCarbon, 'remove_ethane', True),
        ('PRO4 ALA6, remove_ethane=False keeps isolated ALA + warns', 'resid 4 or resid 6', TerminalAlphaCarbon, 'remove_ethane', False),
        ('VAL9 GLY11, remove_methane=True removes isolated GLY', 'resid 9 or resid 11', TerminalAlphaCarbon, 'remove_methane', True),
        ('VAL9 GLY11, remove_methane=False keeps isolated GLY + warns', 'resid 9 or resid 11', TerminalAlphaCarbon, 'remove_methane', False),
        ('THR5 ALA6, remove_ethane=True removes isolated ALA (AlphaCarbon)', 'resid 5 or resid 6', AlphaCarbon, 'remove_ethane', True),
        ('THR5 ALA6, remove_ethane=False keeps isolated ALA + warns (AlphaCarbon)', 'resid 5 or resid 6', AlphaCarbon, 'remove_ethane', False),
    ]
)
def test_check_gly_ala_remove_flags(Test, region_selection, trunc_scheme, remove_flag_name, remove_flag_value, capsys):
    model = GenerateModel(PDB)
    model.set_region(name='region', selection=region_selection)
    qm_method = QMzyme.QM_Method(
        basis_set='6-31G*',
        functional='wB97X-D3',
        qm_input='OPT FREQ',
        program='orca'
    )
    qm_method.assign_to_region(region=model.region)

    isolated_resname = 'ALA' if remove_flag_name == 'remove_ethane' else 'GLY'

    assert any(res.resname == isolated_resname for res in model.region.residues)
    assert model.region.method is not None
    assert model.region.method["basis_set"] == '6-31G*'
    assert model.region.method["functional"] == 'wB97X-D3'

    kwargs = {remove_flag_name: remove_flag_value}
    model.truncate(scheme=trunc_scheme, **kwargs)

    isolated_present = any(res.resname == isolated_resname for res in model.QM_region.residues)

    if remove_flag_value is True:
        # remove_residue() actually pulled the isolated residue out of the region
        assert not isolated_present

    else:
        assert isolated_present
        captured = capsys.readouterr()
        organic_group = 'ethane' if remove_flag_name == 'remove_ethane' else 'methane'
        assert organic_group in captured.out
        assert "may not be an appropriate representation" in captured.out

@pytest.mark.parametrize(
    "Test, region_selection, trunc_scheme, isolated_gly_ala, extend_gly_ala_backbone",[
        ('THR5 ALA6', 'resid 5 or resid 6', TerminalAlphaCarbon, False, True),
        ('PRO4 ALA6', 'resid 4 or resid 6', TerminalAlphaCarbon, True, True),
        ('GLN10 GLY11', 'resid 10 or resid 11', TerminalAlphaCarbon, False, True),
        ('VAL9 GLY11', 'resid 9 or resid 11', TerminalAlphaCarbon, True, True),
        ('THR5 ALA6', 'resid 5 or resid 6', TerminalAlphaCarbon, False, False),
        ('PRO4 ALA6', 'resid 4 or resid 6', TerminalAlphaCarbon, True, False),
        ('GLN10 GLY11', 'resid 10 or resid 11', TerminalAlphaCarbon, False, False),
        ('VAL9 GLY11', 'resid 9 or resid 11', TerminalAlphaCarbon, True, False),
        ('THR5 ALA6', 'resid 5 or resid 6', AlphaCarbon, True, True),
        ('PRO4 ALA6', 'resid 4 or resid 6', AlphaCarbon, True, True),
        ('GLN10 GLY11', 'resid 10 or resid 11', AlphaCarbon, True, True),
        ('VAL9 GLY11', 'resid 9 or resid 11', AlphaCarbon, True, True),
        ('THR5 ALA6', 'resid 5 or resid 6', AlphaCarbon, True, False),
        ('PRO4 ALA6', 'resid 4 or resid 6', AlphaCarbon, True, False),
        ('GLN10 GLY11', 'resid 10 or resid 11', AlphaCarbon, True, False),
        ('VAL9 GLY11', 'resid 9 or resid 11', AlphaCarbon, True, False),
    ]
)
def test_extend_gly_ala_backbone(Test, region_selection, trunc_scheme, isolated_gly_ala, extend_gly_ala_backbone):
    model = GenerateModel(PDB)
    model.set_region(name='region', selection=region_selection)
    qm_method = QMzyme.QM_Method(
        basis_set='6-31G*',
        functional='wB97X-D3',
        qm_input='OPT FREQ',
        program='orca'
    )
    qm_method.assign_to_region(region=model.region)

    assert model.region.n_residues == 2
    assert model.region.method is not None
    assert model.region.method["basis_set"] == '6-31G*'

    # QM_region is only created with truncate(), so it shouldn't exist yet.
    assert not hasattr(model, 'QM_region')

    if isolated_gly_ala == True and trunc_scheme == TerminalAlphaCarbon and extend_gly_ala_backbone == True:
        model.truncate(scheme=trunc_scheme, extend_gly_ala_backbone=extend_gly_ala_backbone)
        for res in model.QM_region.residues:
            resname_by_resid = {res.resid: res.resname for res in model.QM_region.residues}
            atom_names = [atom.name for atom in res.atoms]

            if res.resname in ('GLY', 'ALA'):
                # the isolated residue flagged by _check_gly_ala: with
                # extend_gly_ala_backbone=True it gets bridged with ACE/NME
                # caps and is skipped by the main truncate() loop
                assert resname_by_resid.get(res.resid - 1) == 'ACE'
                assert resname_by_resid.get(res.resid + 1) == 'NME'

                # untouched by truncate(), so native backbone atoms remain
                assert 'N' in atom_names
                assert 'H' in atom_names
                assert 'C' in atom_names
                assert 'O' in atom_names
                assert 'HN' not in atom_names
                assert 'HC' not in atom_names

            elif res.resname == 'ACE':
                # cap_ACE: C/O from preceding residue's backbone, CH3 (was CA), and exactly 3 methyl H's named HH3n
                assert 'C' in atom_names
                assert 'O' in atom_names
                assert 'CH3' in atom_names
                methyl_hs = [n for n in atom_names if n.startswith('HH3')]
                assert len(methyl_hs) == 3
                assert all(atom.resname == 'ACE' for atom in res.atoms)
                assert all(atom.charge == 0 for atom in res.atoms)

            elif res.resname == 'NME':
                # cap_NME: N + CH3 (was CA) always present; amide H present unless the capped residue is Proline; exactly 3 methyl H's
                assert 'N' in atom_names
                assert 'CH3' in atom_names
                methyl_hs = [n for n in atom_names if n.startswith('HH3')]
                assert len(methyl_hs) == 3
                assert all(atom.resname == 'NME' for atom in res.atoms)
                assert all(atom.charge == 0 for atom in res.atoms)

            else:
                # every other protein residue truncated normally by the scheme
                if res.resid == 1:
                    assert 'N' in atom_names
                elif res.resname != 'PRO':
                    assert 'N' not in atom_names
                    assert 'HN' in atom_names
                else:
                    assert 'N' in atom_names
                    assert 'HN' in atom_names

                assert 'C' not in atom_names
                assert 'O' not in atom_names
                assert 'HC' in atom_names

        assert model.QM_region.truncated is True

    if isolated_gly_ala == True and extend_gly_ala_backbone == True and trunc_scheme == AlphaCarbon \
                and region_selection in ('resid 4 or resid 6', 'resid 9 or resid 11'):
            model.truncate(scheme=trunc_scheme, extend_gly_ala_backbone=extend_gly_ala_backbone)
            for res in model.QM_region.residues:
                resname_by_resid = {r.resid: r.resname for r in model.QM_region.residues}
                atom_names = [atom.name for atom in res.atoms]

                if res.resname in ('GLY', 'ALA'):
                    assert resname_by_resid.get(res.resid - 1) == 'ACE'
                    assert resname_by_resid.get(res.resid + 1) == 'NME'
                    assert 'N' in atom_names
                    assert 'H' in atom_names
                    assert 'C' in atom_names
                    assert 'O' in atom_names

                elif res.resname == 'ACE':
                    assert 'CH3' in atom_names
                    assert len([n for n in atom_names if n.startswith('HH3')]) == 3

                elif res.resname == 'NME':
                    assert 'CH3' in atom_names
                    assert len([n for n in atom_names if n.startswith('HH3')]) == 3

                else:
                    # normal residue, AlphaCarbon truncates both sides
                    if res.resid == 1:
                        assert 'N' in atom_names
                    elif res.resname != 'PRO':
                        assert 'N' not in atom_names
                        assert 'HN' in atom_names
                    else:
                        assert 'N' in atom_names
                        assert 'HN' in atom_names
                    assert 'C' not in atom_names
                    assert 'HC' in atom_names

            assert model.QM_region.truncated is True

    if isolated_gly_ala == False:
        model.truncate(scheme=trunc_scheme, extend_gly_ala_backbone=extend_gly_ala_backbone)

@pytest.mark.parametrize(
    "Test, region_selection, override_truncation",[
        ('Undecided (None): halts with ValueError', 'resid 1 or resid 2', None),
        ('override_truncation=False: skips already-truncated residues', 'resid 1 or resid 2', False),
        ('override_truncation=True: re-truncates already-truncated residues', 'resid 1 or resid 2', True),
    ]
)
def test_check_override_truncation(Test, region_selection, override_truncation):
    model = GenerateModel(PDB)
    model.set_region(name='region', selection=region_selection)
    qm_method = QMzyme.QM_Method(
        basis_set='6-31G*',
        functional='wB97X-D3',
        qm_input='OPT FREQ',
        program='orca'
    )
    qm_method.assign_to_region(region=model.region)

    assert not hasattr(model, 'QM_region')

    model.truncate(scheme=TerminalAlphaCarbon)
    truncation_params_before = {
        res.resid: res.truncation_params for res in model.QM_region.residues
    }
    assert all(v is not None for v in truncation_params_before.values())

    if override_truncation is None:
        with pytest.raises(ValueError):
            model.truncate(scheme=TerminalAlphaCarbon, override_truncation=override_truncation)

    if override_truncation is False:
        model.remove_region('QM_region')
        model.truncate(
            scheme=TerminalAlphaCarbon,
            override_truncation=override_truncation,
            override_capping=False,
        )
        for res in model.QM_region.residues:
            assert res.truncation_params == truncation_params_before[res.resid]

    if override_truncation is True:
        model.remove_region('QM_region')
        model.truncate(scheme=TerminalAlphaCarbon, override_truncation=override_truncation)
        for res in model.QM_region.residues:
            assert res.truncation_params is not None


@pytest.mark.parametrize(
    "Test, region_selection, override_capping",[
        ('Undecided (None): halts with ValueError', 'resid 4 or resid 6', None),
        ('override_capping=False: skips already-capped residue', 'resid 4 or resid 6', False),
        ('override_capping=True: removes existing caps and re-caps fresh', 'resid 4 or resid 6', True),
    ]
)
def test_check_override_capping(Test, region_selection, override_capping):
    model = GenerateModel(PDB)
    model.set_region(name='region', selection=region_selection)
    qm_method = QMzyme.QM_Method(
        basis_set='6-31G*',
        functional='wB97X-D3',
        qm_input='OPT FREQ',
        program='orca'
    )
    qm_method.assign_to_region(region=model.region)

    model.truncate(scheme=TerminalAlphaCarbon, extend_gly_ala_backbone=True)
    capped_res = next(res for res in model.QM_region.residues if res.resname == 'ALA')
    assert capped_res.capping_scheme is not None

    if override_capping is None:
        with pytest.raises(ValueError):
            model.truncate(
                scheme=TerminalAlphaCarbon,
                extend_gly_ala_backbone=True,
                override_truncation=True,
                override_capping=override_capping,
            )

    if override_capping is False:
        model.remove_region('QM_region')
        model.truncate(
            scheme=TerminalAlphaCarbon,
            extend_gly_ala_backbone=True,
            override_truncation=True,
            override_capping=override_capping,
        )
        still_capped = next(res for res in model.QM_region.residues if res.resname == 'ALA')
        assert still_capped.capping_scheme is not None

    if override_capping is True:
        model.remove_region('QM_region')
        model.truncate(
            scheme=TerminalAlphaCarbon,
            extend_gly_ala_backbone=True,
            override_truncation=True,
            override_capping=override_capping,
        )
        recapped = next(res for res in model.QM_region.residues if res.resname == 'ALA')
        assert recapped.capping_scheme is not None

def test_check_override_capping_partial(capsys):
    model = GenerateModel(PDB)
    model.set_region(name='region', selection='resid 2')
    qm_method = QMzyme.QM_Method(
        basis_set='6-31G*', functional='wB97X-D3', qm_input='OPT FREQ', program='orca'
    )
    qm_method.assign_to_region(region=model.region)

    model.region.add_N_terminus_ACE(2)
    model.truncate(
        scheme=TerminalAlphaCarbon,
        override_truncation=True,
        override_capping=False,
    )
    captured = capsys.readouterr()
    assert "only one terminus capped" in captured.out

    res_present = next(r for r in model.QM_region.residues if r.resid == 2)
    assert res_present.capping_scheme == 'cap_ACE, cap_H'

    model.QM_region.remove_region('QM_region') if False else model.remove_region('QM_region')
    model.truncate(
        scheme=TerminalAlphaCarbon,
        override_truncation=True,
        override_capping=True,
    )
    recapped = next(r for r in model.QM_region.residues if r.resid == 2)
    assert recapped.capping_scheme == 'cap_H'

@pytest.mark.parametrize(
    "Test, region_selection, setup, override_truncation, override_capping, expected_skip_resids",[
        ('Fresh region: nothing to skip', 'resid 1 or resid 2', 'none', None, None, set()),
        ('override_truncation=False: already-truncated resids are skipped', 'resid 1 or resid 2', 'pretruncate', False, None, {1, 2}),
        ('override_truncation=True: already-truncated resids are reset and reprocessed, not skipped', 'resid 1 or resid 2', 'pretruncate', True, None, set()),
        ('override_capping=False, only N-term capped (ACE): partial cap is not skipped', 'resid 2', 'partial_cap', True, False, set()),
        ('override_capping=False, both termini capped (ACE+NME): fully-capped resid is skipped', 'resid 2', 'full_cap', True, False, {2}),
    ]
)
def test_skip_resids(Test, region_selection, setup, override_truncation, override_capping, expected_skip_resids):

    model = GenerateModel(PDB)
    model.set_region(name='region', selection=region_selection)
    qm_method = QMzyme.QM_Method(
        basis_set='6-31G*', functional='wB97X-D3', qm_input='OPT FREQ', program='orca'
    )
    qm_method.assign_to_region(region=model.region)

    region = model.region

    if setup == 'none':
        # Nothing done yet; nothing should be truncated or capped.
        assert region.resids == [1, 2]
        for res in region.residues:
            assert res.truncation_params is None
            assert res.capping_scheme is None

    elif setup == 'pretruncate':
        region = TerminalAlphaCarbon(region, name='first_pass').return_region()

        assert region.resids == [1, 2]
        assert region.get_residue(1).truncation_params == 'TerminalAlphaCarbon'
        assert region.get_residue(2).truncation_params == 'TerminalAlphaCarbon'

    elif setup == 'partial_cap':
        region.add_N_terminus_ACE(2)

        assert region.resids == [1, 2]
        assert region.get_residue(1).resname == 'ACE'
        assert region.get_residue(2).capping_scheme == 'cap_ACE'
        assert region.get_residue(2).truncation_params is None

    elif setup == 'full_cap':
        region.add_N_terminus_ACE(2)
        region.add_C_terminus_NME(2)

        assert region.resids == [1, 2, 3]
        assert region.get_residue(1).resname == 'ACE'
        assert region.get_residue(3).resname == 'NME'
        assert region.get_residue(2).capping_scheme == 'cap_ACE, cap_NME'

    scheme = TerminalAlphaCarbon(
        region, name='second_pass',
        override_truncation=override_truncation,
        override_capping=override_capping,
    )

    assert scheme.skip_resids == expected_skip_resids

@pytest.mark.parametrize(
    "Test, init_file, region_selection",[
        ('First and last residue in protein: MET1 GLN262', PDB, 'resid 1 or resid 262'),
        ('MET1 ASN2', PDB, 'resid 1 or resid 2'),
        ('MET1 LEU3', PDB, 'resid 1 or resid 3'),
        ('ASN2 THR5', PDB, 'resid 2 or resid 5'),
        ('ASN2 LEU3 THR5 ALA6', PDB, 'resid 2 or resid 3 or resid 5 or resid 6'),
        ('PRO4 THR5', PDB, 'resid 4 or resid 5'),
        ('LEU3 PRO4', PDB, 'resid 3 or resid 4'),
        ('With Non protein residue: WAT265', PDB, 'resid 3 or resid 265'),
    ]
)

def test_AlphaCarbon(Test, init_file, region_selection, truncation_scheme=AlphaCarbon):
    model = GenerateModel(init_file)
    model.set_region(name='region', selection=region_selection)
    region_truncated = truncation_scheme(model.region, name="region_truncated", remove_ethane=False).return_region()
    assert region_truncated != model.region

    for i, res in enumerate(model.region.residues):
        truncated_res = region_truncated.residues[i]
        res_atoms = [atom.name for atom in res.atoms]
        truncated_res_atoms = [atom.name for atom in truncated_res.atoms]
        if res.resid == 1:
            assert 'N' in res_atoms
            assert 'N' in truncated_res_atoms
        elif res.resname != 'PRO':
            assert 'N' in res_atoms
            assert 'H' in res_atoms
            assert 'N' not in truncated_res_atoms
            assert 'HN' in truncated_res_atoms
        elif res.resname == 'PRO':
            assert 'N' in res_atoms
            assert 'H' not in res_atoms
            assert 'N' in truncated_res_atoms
            assert 'HN' in truncated_res_atoms
        assert 'C' in res_atoms
        assert 'O' in res_atoms
        if res.resid != 262:
            assert 'C' not in truncated_res_atoms
            assert 'O' not in truncated_res_atoms
            assert 'HC' in truncated_res_atoms
        elif res.resid == 262:
            assert 'C' in truncated_res_atoms
            assert 'O' in truncated_res_atoms


@pytest.mark.parametrize(
    "Test, init_file, region_selection",[
        ('First and last residue in protein: MET1 GLN262', PDB, 'resid 1 or resid 262'),
        ('MET1 ASN2', PDB, 'resid 1 or resid 2'),
        ('MET1 LEU3', PDB, 'resid 1 or resid 3'),
        ('ASN2 THR5', PDB, 'resid 2 or resid 5'),
        ('ASN2 LEU3 THR5 ALA6', PDB, 'resid 2 or resid 3 or resid 5 or resid 6'),
        ('PRO4 THR5', PDB, 'resid 4 or resid 5'),
        ('LEU3 PRO4', PDB, 'resid 3 or resid 4'),
        ('With Non protein residue: WAT265', PDB, 'resid 3 or resid 265'),
    ]
)
def test_TerminalAlphaCarbon(Test, init_file, region_selection, truncation_scheme=TerminalAlphaCarbon):
    model = GenerateModel(init_file)
    model.set_region(name='region', selection=region_selection)
    region_truncated = truncation_scheme(model.region, name="region_truncated").return_region()
    #model.truncate()
    #region_truncated = model.truncated
    assert region_truncated != model.region

    #First check that the original region didn't change:
    original_first_res = model.region.residues[0]
    truncated_first_res = region_truncated.residues[0]
    original_last_res = model.region.residues[-1]
    truncated_last_res = region_truncated.residues[-1]

    if original_first_res.resname != 'PRO':
        removed_atom_name = 'H'
        if original_first_res.resid == 1:
            assert 'H1' in [atom.name for atom in original_first_res.atoms]
            assert 'H2' in [atom.name for atom in original_first_res.atoms]
            assert 'H3' in [atom.name for atom in original_first_res.atoms]
            assert 'H2' in [atom.name for atom in truncated_first_res.atoms]
            assert 'H3' in [atom.name for atom in truncated_first_res.atoms]
        if original_first_res.resid != 1:
            assert 'H' in [atom.name for atom in original_first_res.atoms]
            assert 'N' in [atom.name for atom in original_first_res.atoms]
            assert 'H' not in [atom.name for atom in truncated_first_res.atoms]
            assert 'N' not in [atom.name for atom in truncated_first_res.atoms]
            assert 'HN' in [atom.name for atom in truncated_first_res.atoms]
    
    if original_first_res.resname == 'PRO':
        assert 'N' in [atom.name for atom in original_first_res.atoms]
        assert 'H' not in [atom.name for atom in original_first_res.atoms]
        assert 'N' in [atom.name for atom in truncated_first_res.atoms]
        assert 'HN' in [atom.name for atom in truncated_first_res.atoms]

    for i in range(model.region.n_residues-1):
        resid = model.region.resids[i]
        next_resid = model.region.resids[i+1] 
        original_res = model.region.get_residue(resid)
        original_next_res = model.region.get_residue(next_resid)
        truncated_res = region_truncated.get_residue(resid)
        truncated_next_res = region_truncated.get_residue(next_resid)
        assert 'C' in [atom.name for atom in original_res.atoms]
        assert 'O' in [atom.name for atom in original_res.atoms]
        assert 'N' in [atom.name for atom in original_next_res.atoms]
        if original_next_res.resname != 'PRO':
            assert 'H' in [atom.name for atom in original_next_res.atoms]
        if resid+1 == next_resid:
            #C term is not removed
            assert 'C' in [atom.name for atom in truncated_res.atoms]
            assert 'O' in [atom.name for atom in truncated_res.atoms]
            assert 'N' in [atom.name for atom in truncated_next_res.atoms]
            if original_next_res.resname != 'PRO':
                assert 'H' in [atom.name for atom in truncated_next_res.atoms]
        if resid+1 != next_resid:
            #C term is removed
            assert 'C' not in [atom.name for atom in truncated_res.atoms]
            assert 'O' not in [atom.name for atom in truncated_res.atoms]
            assert 'H' not in [atom.name for atom in truncated_next_res.atoms]
            if original_next_res.resname == 'PRO':
                assert 'N' in [atom.name for atom in truncated_next_res.atoms]
            elif original_next_res.resname != 'PRO':
                assert 'N' not in [atom.name for atom in truncated_next_res.atoms]
            assert 'HN' in [atom.name for atom in truncated_next_res.atoms]
            assert 'HC' in [atom.name for atom in truncated_res.atoms]
            names = [atom.name for atom in truncated_res.atoms]
            if 'N' not in names and 'C' not in names:
                print(original_res.atoms)
                assert 'HN' in names
                assert 'HC' in names

@pytest.mark.parametrize(
    "Test, init_file, region_selection, truncation_selection",[
        ('MET1', PDB, 'resid 1', 'resid 1'),
        ('MET1 LEU3', PDB, 'resid 1 or resid 3', 'resid 1 or resid 3'),
        ('PRO4 (Should Skip)', PDB, 'resid 4', 'resid 4'),
        ('GLY11 (Should Skip)', PDB, 'resid 11', 'resid 11'),
        ('Selective Truncation', PDB, 'resid 1:20', 'resid 5') 
    ]
)
def test_BetaCarbon(Test, init_file, region_selection, truncation_selection):
    model = GenerateModel(init_file)
    model.set_region(name='region', selection=region_selection)
    
    # Capture original atom names AND resnames BEFORE truncation
    original_state = {}
    original_resnames = {}
    for res in model.region.residues:
        original_state[res.resid] = [atom.name for atom in res.atoms]
        original_resnames[res.resid] = res.resname
    model.region.universe = model.universe

    model.set_region(name='trunc_sele', selection=truncation_selection)
    target_resids = [r.resid for r in model.trunc_sele.residues]

    # Run Truncation
    truncated = BetaCarbon(region=model.region, selection=truncation_selection, name='truncated_reg')
    region_truncated = truncated.return_region()

    assert region_truncated is not None
    
    for trunc_res in region_truncated.residues:
        resid = trunc_res.resid
        orig_resname = original_resnames[resid]
        orig_names = original_state[resid]
        trunc_names = [atom.name for atom in trunc_res.atoms]

        if resid in target_resids:
            assert trunc_res.truncation_params == 'BetaCarbon'
        else:
            assert trunc_res.truncation_params is None

        # Skip Logic (GLY/PRO/ALA all early-return in truncate(), or resid
        # simply wasn't part of the truncation selection)
        if orig_resname in ("GLY", "PRO", "ALA") or resid not in target_resids:
            assert len(trunc_names) == len(orig_names)
            assert trunc_res.resname == orig_resname
            continue

        assert len(trunc_names) != len(orig_names)

        # BetaCarbon relabels truncated residues (and all their atoms) as ALA
        assert trunc_res.resname == 'ALA'
        assert all(atom.resname == 'ALA' for atom in trunc_res.atoms)

        # Check Scaffold Preservation (N, CA, C, O, CB)
        for atom in ['N', 'CA', 'C', 'O', 'CB']:
            if atom in orig_names:
                assert atom in trunc_names
        
        # Verify that heavy atoms beyond CB are removed
        for side_chain_atom in ['CG', 'SD', 'OG', 'CD', 'CE']:
            if side_chain_atom in orig_names:
                assert side_chain_atom not in trunc_names

        # Check for Capping (Presence of new Hydrogens from cap_H)
        new_atoms = [n for n in trunc_names if n not in orig_names]
        assert len(new_atoms) > 0

def test_combine_truncation():
    model = GenerateModel(PDB)
    model.set_region(name='region', selection='resid 263 or resid 16 or resid 17 or resid 57')
    truncation = TerminalAlphaCarbon(region=model.region, name=None)
    truncate = truncation.return_region()
    assert {a.resid for a in truncate.atoms} == {263, 16, 17, 57}
    truncate.write('truncated_test')

    model.set_region(name='region1', selection='resid 263 or resid 16')
    model.set_region(name='region2', selection='resid 17 or resid 57')
    
    combined = model.region1.combine(model.region2)

    assert combined.n_atoms == model.region1.n_atoms + model.region2.n_atoms
 
    truncation = TerminalAlphaCarbon(region=combined, name=None)
    truncated_combined = truncation.return_region()
    assert {a.resid for a in truncated_combined.atoms} == {263, 16, 17, 57}
    truncated_combined.write('truncated_combined_test')

    assert truncate.n_atoms == truncated_combined.n_atoms
    assert {a.id for a in truncate.atoms} == {a.id for a in truncated_combined.atoms}
    for atom in truncate.atoms:
        assert atom.name == truncated_combined.get_atom(atom.id).name

    restore_directory()
    assert 'truncated_combined_test.pdb' not in os.listdir()
    assert 'truncated_test.pdb' not in os.listdir()