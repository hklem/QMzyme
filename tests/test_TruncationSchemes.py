"""
Tests for the QMzyme truncation_schemes.py and truncation_utils.py codes.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow

from QMzyme.GenerateModel import GenerateModel
from QMzyme.TruncationSchemes import *
import pytest
from QMzyme.data import PDB


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
    region_truncated = truncation_scheme(model.region, name=None).return_region()
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
    region_truncated = truncation_scheme(model.region, name=None).return_region()
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
    
    # Capture original atom names BEFORE truncation manually
    original_state = {}
    for res in model.region.residues:
        original_state[res.resid] = [atom.name for atom in res.atoms]
    
    ala_group = model.universe.select_atoms(truncation_selection)
    
    # Run Truncation
    truncator = BetaCarbon(region=model.region, ala_atom_group=ala_group, model=model, name='truncated_reg')
    region_truncated = truncator.return_region()

    assert region_truncated is not None
    
    for trunc_res in region_truncated.residues:
        resid = trunc_res.resid
        resname = trunc_res.resname
        
        orig_names = original_state[resid]
        trunc_names = [atom.name for atom in trunc_res.atoms]

        # Skip Logic (PRO/GLY or not in ala_group)
        if resname in ["GLY", "PRO"] or resid not in ala_group.resids:
            assert len(trunc_names) == len(orig_names)
            continue

        # Active Truncation
        # Detecting change in atom numbers
        assert len(trunc_names) != len(orig_names)

        # Check Scaffold Preservation (N, CA, C, O, CB)
        for atom in ['N', 'CA', 'C', 'O', 'CB']:
            if atom in orig_names:
                assert atom in trunc_names
        
        # Verify that heavy atoms beyond CB are removed
        for sc_atom in ['CG', 'SD', 'OG', 'CD', 'CE']:
            if sc_atom in orig_names:
                assert sc_atom not in trunc_names

        # Check for Capping (Presence of new Hydrogens from cap_H)
        new_atoms = [n for n in trunc_names if n not in orig_names]
        assert len(new_atoms) > 0