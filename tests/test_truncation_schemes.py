"""
Tests for the QMzyme truncation_schemes.py and truncation_utils.py codes.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow

from QMzyme.GenerateModel import GenerateModel
import pytest
from importlib_resources import files

pdb_file = str(files('QMzyme.data').joinpath('1oh0.pdb'))

@pytest.mark.parametrize(
    "Test, init_file, region_selection",[
        #('First and last residue in protein: MET1 GLN262', pdb_file, 'resid 1 or resid 262'),
        ('MET1 ASN2', pdb_file, 'resid 1 or resid 2'),
        #('MET1 LEU3', pdb_file, 'resid 1 or resid 3'),
        #('ASN2 THR5', pdb_file, 'resid 2 or resid 5'),
        #('ASN2 LEU3 THR5 ALA6', pdb_file, 'resid 2 or resid 3 or resid 5 or resid 6'),
        #('PRO4 THR5', pdb_file, 'resid 4 or resid 5'),
        #('LEU3 PRO4', pdb_file, 'resid 3 or resid 4'),
        #('Non protein residue: WAT265', pdb_file, 'resid 265'),
    ]
)
def test_truncate_region_CA_terminal(Test, init_file, region_selection, truncation_scheme="CA_terminal"):
    model = GenerateModel(file=init_file)
    original_region = model.set_region('test', region_selection)
    truncated_region = model.truncate_region(original_region, truncation_scheme)
    #First check that the original region didn't change:
    original_first_res = original_region.residues[0]
    truncated_first_res = truncated_region.residues[0]
    original_last_res = original_region.residues[-1]
    truncated_last_res = truncated_region.residues[-1]

    if original_first_res.resname != 'PRO':
        removed_atom_name = 'H'
        # if original_first_res.resid == 1:
        #     assert all(['H2', 'H3']) in [atom.name for atom in original_first_res.atoms]
        #     assert all(['H2', 'H3']) in [atom.name for atom in truncated_first_res.atoms]
        if original_first_res.resid != 1:
            assert removed_atom_name in [atom.name for atom in original_first_res.atoms]
            assert removed_atom_name not in [atom.name for atom in truncated_first_res.atoms]
    
    if original_first_res.resname == 'PRO':
        assert 'N' in [atom.name for atom in original_first_res.atoms]
        assert 'H' not in [atom.name for atom in original_first_res.atoms]
        assert 'N' in [atom.name for atom in truncated_first_res.atoms]
    
    assert 'H1' in [atom.name for atom in truncated_first_res.atoms]

    for i in enumerate(original_region.resids[:-1]):
        resid = original_region.resids[i]
        next_resid = original_region.resids[i+1] 
        original_res = original_region.get_residue(resid)
        original_next_res = original_region.get_residue(next_resid)
        truncated_res = truncated_region.get_residue(resid)
        truncated_next_res = truncated_region.get_residue(next_resid)
        assert 'C' in [atom.name for atom in original_res.atoms]
        assert 'O' in [atom.name for atom in original_res.atoms]
        assert 'N' in [atom.name for atom in original_next_res.atoms]
        if original_res.resname != 'PRO':
            assert 'H' in [atom.name for atom in original_next_res.atoms]
        if resid+1 == next_resid:
            #C term is not removed
            assert 'C' in [atom.name for atom in truncated_res.atoms]
            assert 'O' in [atom.name for atom in truncated_res.atoms]
            assert 'N' in [atom.name for atom in truncated_next_res.atoms]
            if original_res.resname != 'PRO':
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
            assert 'H1' in [atom.name for atom in truncated_next_res.atoms]
            assert 'H1' in [atom.name for atom in truncated_res.atoms]


