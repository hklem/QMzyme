import pytest
import QMzyme
from QMzyme.truncation_utils import cap_ACE, cap_NME, add_removed_atoms, remove_added_atoms
from QMzyme.GenerateModel import GenerateModel
from QMzyme.data import PDB
 
 
@pytest.mark.parametrize(
    "Test, region_selection, replace_atom, fixed_atom",[
        ('MET1, C_termina of protein', 'resid 1', 'resid 1 and name N', 'resid 1 and name CA'),
        ('GLN262, N-terminal of protein', 'resid 262', 'resid 262 and name N', 'resid 262 and name CA'),
        ('ASN2', 'resid 2', 'resid 2 and name N', 'resid 2 and name CA'),
        ('LEU3, with PRO4 as C-termini', 'resid 3', 'resid 3 and name N', 'resid 3 and name CA'),
        ('THR5, with PRO4 as N-termini', 'resid 5', 'resid 5 and name N', 'resid 5 and name CA'),
        ('GLN10, with GLY11 as C-termini', 'resid 10', 'resid 10 and name N', 'resid 10 and name CA'),
        ('GLN12, with GLY11 as N-termini', 'resid 12', 'resid 12 and name N', 'resid 12 and name CA'),
        ('LEU3, PRO4, THR5', 'resid 3 or resid 4 or resid 5', 'resid 5 and name N', 'resid 5 and name CA'),
    ]
)
def test_cap_ACE(Test, region_selection, replace_atom, fixed_atom):
    model = GenerateModel(PDB)
    model.set_region(name='region', selection=region_selection)

    resid = int(replace_atom.split()[1])
    universe = model.universe
    replace_atom = universe.select_atoms(replace_atom)[0]
    fixed_atom = universe.select_atoms(fixed_atom)[0]

    if resid == 1:
        with pytest.raises(ValueError, match="no preceding residue found"):
            cap_ACE(replace_atom, universe)
        return

    cap_atoms = cap_ACE(replace_atom, universe)
    cap_atom_names = [atom.name for atom in cap_atoms]

    assert 'C' in cap_atom_names
    assert 'O' in cap_atom_names
    assert 'CH3' in cap_atom_names

    methyl_h_names = [n for n in cap_atom_names if n.startswith('HH3')]
    assert len(methyl_h_names) == 3
    assert sorted(methyl_h_names) == ['HH31', 'HH32', 'HH33']

    assert all(atom.charge == 0 for atom in cap_atoms)
    assert all(atom.resname == "ACE" for atom in cap_atoms)


@pytest.mark.parametrize(
    "Test, region_selection, replace_atom, fixed_atom",[
        ('MET1, C_termina of protein', 'resid 1', 'resid 1 and name C', 'resid 1 and name CA'),
        ('GLN262, N-terminal of protein', 'resid 262', 'resid 262 and name C', 'resid 262 and name CA'),
        ('ASN2', 'resid 2', 'resid 2 and name C', 'resid 2 and name CA'),
        ('LEU3, with PRO4 as C-termini', 'resid 3', 'resid 3 and name C', 'resid 3 and name CA'),
        ('THR5, with PRO4 as N-termini', 'resid 5', 'resid 5 and name C', 'resid 5 and name CA'),
        ('GLN10, with GLY11 as C-termini', 'resid 10', 'resid 10 and name C', 'resid 10 and name CA'),
        ('GLN12, with GLY11 as N-termini', 'resid 12', 'resid 12 and name C', 'resid 12 and name CA'),
        ('LEU3, PRO4, THR5', 'resid 3 or resid 4 or resid 5', 'resid 3 and name C', 'resid 3 and name CA'),
    ]
)
def test_cap_NME(Test, region_selection, replace_atom, fixed_atom):
    model = GenerateModel(PDB)
    model.set_region(name='region', selection=region_selection)

    resid = int(replace_atom.split()[1])
    universe = model.universe
    replace_atom = universe.select_atoms(replace_atom)[0]
    fixed_atom = universe.select_atoms(fixed_atom)[0]

    if resid == 3:
        with pytest.raises(ValueError, match="Proline"):
            cap_NME(replace_atom, universe)
        return

    if resid == 262:
        with pytest.raises(ValueError, match="no following residue found"):
            cap_NME(replace_atom, universe)
        return

    cap_atoms = cap_NME(replace_atom, universe)
    cap_atom_names = [atom.name for atom in cap_atoms]

    assert 'N' in cap_atom_names
    assert 'H' in cap_atom_names
    assert 'CH3' in cap_atom_names

    methyl_h_names = [n for n in cap_atom_names if n.startswith('HH3')]
    assert len(methyl_h_names) == 3
    assert sorted(methyl_h_names) == ['HH31', 'HH32', 'HH33']

    assert all(atom.charge == 0 for atom in cap_atoms)
    assert all(atom.resname == "NME" for atom in cap_atoms)
 

def test_truncation_add_remove_atoms():
    model = GenerateModel(PDB)
    model.set_region(name='region', selection='resid 2')
    model.region.truncate(scheme=QMzyme.TruncationSchemes.TerminalAlphaCarbon, selection='resid 2')

    residue = model.region.get_residue(2)

    added = residue.added_atoms
    removed = residue.removed_atoms

    assert len(added) == 2
    assert len(removed) == 4

    atom_count_post_truncation = len(residue.atoms)

    # Test remove_added_atoms
    remove_added_atoms(residue)
    remaining_ids = {atom.id for atom in residue.atoms}
    assert all(atom.id not in remaining_ids for atom in added)
    assert len(residue.atoms) == atom_count_post_truncation - len(added)
    assert residue.added_atoms == []

    # Test add_removed_atoms
    atom_count_pre_readd = len(residue.atoms)
    add_removed_atoms(residue)
    assert len(residue.atoms) == atom_count_pre_readd + len(removed)
    assert residue.removed_atoms == []