"""
Tests for the QMzyme RegionBuilder.py ands QMzymeRegion.py code.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow

import pytest
import QMzyme
from QMzyme.RegionBuilder import RegionBuilder
import MDAnalysis as mda
from QMzyme.data import PDB
from QMzyme import GenerateModel

u = mda.Universe(PDB)
atom_group = u.select_atoms('resid 2-5')
id1 = atom_group.atoms[0].id

def test_RegionBuilder():
    rb1 = RegionBuilder(name='test')
    assert rb1.__repr__() == "<RegionBuilder: Current QMzymeRegion, test, "+\
                             "contains 0 atom(s) and 0 residue(s)>"
    
    rb1.init_atom_group(atom_group=atom_group)
    assert rb1.__repr__() == "<RegionBuilder: Current QMzymeRegion, test, "+\
                             "contains 61 atom(s) and 4 residue(s)>"

def test_QMzymeRegion():
    region_builder = RegionBuilder(name='test')
    assert region_builder.__repr__() == "<RegionBuilder: Current QMzymeRegion, test, "+\
                                        "contains 0 atom(s) and 0 residue(s)>"
    
    region_builder.init_atom_group(atom_group=atom_group)
    assert region_builder.__repr__() == "<RegionBuilder: Current QMzymeRegion, test, "+\
                                        "contains 61 atom(s) and 4 residue(s)>"
    
    # test region was populated as expected
    region = region_builder.get_region()
    assert region.__repr__() == "<QMzymeRegion test contains 61 atom(s) and 4 residue(s)>"
    assert region.has_atom(id=id1)
    assert region.has_residue(resid=3)
    #assert region.atom_group == atom_group
    assert any(region.ids) == any(atom_group.ids)
    assert any(region.resids) == any(atom_group.resids)

    # add atom through region builder init_atom() method
    mda_atom = u.select_atoms('resid 1 and name CA').atoms[0] # has id=5
    region_builder.init_atom(mda_atom)
    region = region_builder.get_region()
    qmz_atom = region.get_atom(id=5)
    assert region.n_atoms == 62
    assert 1 in region.resids
    assert 5 in region.ids

    region_builder.init_atom(mda_atom)
    with pytest.raises(UserWarning): 
        region = region_builder.get_region() # Because this atom already exists in region.
    # remove that problem atom from region_builder atoms to continue with testing
    region_builder.atoms = region_builder.atoms[:-1]
    
    res = region.get_residue(resid=qmz_atom.resid)
    assert f"{qmz_atom.element}1" not in [atom.name for atom in res.atoms] # check it doesn't exist first
    # now add the atom again- it will be changed because it was not unique and not an mda_atom with immutable id.
    region_builder.init_atom(qmz_atom)
    assert qmz_atom != region_builder.atoms[-1] # atom has changed because it was already there
    assert qmz_atom == region_builder.atoms[-2]
    region = region_builder.get_region()
    assert region.n_atoms == 63 
    res = region.get_residue(resid=qmz_atom.resid)
    assert f"{qmz_atom.element}1" in [atom.name for atom in res.atoms]
    assert region_builder.atoms[-1].id == max(region.get_residue(qmz_atom.resid).ids)

    # test getting atom ids for all CA atoms
    ids = region.get_ids(attribute='name', value='CA')
    assert sorted(ids) == [5, 22, 36, 63, 69]

    # test setting fixed atoms
    region.set_fixed_atoms(ids=ids)
    for id in ids:
        assert region.get_atom(id).is_fixed == True

def test_add_regions():
    model = QMzyme.GenerateModel(PDB)
    model.set_region(name='r1', selection='resid 263')
    model.set_region(name='r2', selection='resid 103')
    model.set_region(name='r1_r2', selection='resid 103 or resid 263')

    assert model.r1.n_atoms + model.r2.n_atoms == model.r1_r2.n_atoms
    assert model.r1 + model.r2 == model.r1_r2
    assert model.r1_r2 + model.r1 == model.r1_r2
    model.r1.set_fixed_atoms(ids=model.r1.ids)

    r3 = model.r1_r2 + model.r1
    for atom in r3.atoms:
        assert atom.is_fixed == False

    r3 = model.r1 + model.r1_r2
    for atom in r3.atoms:
        if atom.id in model.r1.ids:
            assert atom.is_fixed == True
        else:
            assert atom.is_fixed == False

def test_subtract_regions():
    model = QMzyme.GenerateModel(PDB)
    model.set_region(name='r1', selection='resid 263')
    model.set_region(name='r2', selection='resid 103')
    model.set_region(name='r1_r2', selection='resid 103 or resid 263')

    assert model.r1_r2 - model.r1 == model.r2
    assert (model.r1 - model.r1_r2).n_atoms == 0


def test_equal_regions():
    model = QMzyme.GenerateModel(PDB)
    model.set_region(name='r1', selection='resid 263')
    model.set_region(name='r2', selection='resid 263')
    assert model.r1 == model.r2

    setattr(model.r1.atoms[0], "name", "X")
    assert model.r1 != model.r2


def test_QMzymeResidue():
    region_builder = RegionBuilder(name='test')
    region_builder.init_atom_group(atom_group=atom_group)
    region = region_builder.get_region()
    residue = region.residues[0]
    assert residue.__repr__() == "<QMzymeResidue resname: ASN, resid: 2, chain: A>"
    assert residue.get_atom('CA').__repr__() == "<QMzymeAtom 22: CA of resname ASN, resid 2>"
    assert residue.chain == 'A'
    assert not residue.has_atom(100000000)

def test_get_overlapping_atoms():
    model = GenerateModel(PDB)
    model.set_region(name='reg_A', selection='resid 1-3')
    model.set_region(name='reg_B', selection='resid 3-5')
    
    reg_A = model.reg_A
    reg_B = model.reg_B

    # Test overlapping atoms (Residue 3 intersection)
    overlap = reg_A.get_overlapping_atoms(reg_B)
    assert len(overlap) > 0
    assert any(a.resid == 3 for a in overlap)

    # Test summary dictionary (Fixed keys based on FAILED output)
    summary = reg_A.summarize()
    assert 'Resname' in summary
    assert 'Resid' in summary
    assert 'Charge' in summary
    assert 'MET' in summary['Resname']

    # Test residue-level access and backbone retrieval
    res = reg_A.residues[0]
    bb = res.get_backbone_atoms()
    # Standard protein backbone consists of N, CA, C, O
    assert len(bb) >= 4 
    
    # Verify atom naming
    assert res.get_atom("CA").name == "CA"

def test_guess_charge():

    model = GenerateModel(PDB)
    model.set_region(name='test', selection='resid 1')
    res = model.test.residues[0] 

    assert "MET" in repr(res)
    res.set_chain("A")
    assert res.get_atom("CA").name == "CA"
    
    # Touch backbone and charge guessing
    assert len(res.get_backbone_atoms()) >= 4 
    res.guess_charge(verbose=False)
    assert res.charge is not None

def test_find_nearby_residues(capsys):
    model = GenerateModel(PDB)
    model.set_region(selection="all", name="full_region")
    model.set_region(name='test', selection='resid 1')
    model.test.find_nearby_residues(model.full_region, 1.5)

    captured = capsys.readouterr()

    assert "is within 1.5" in captured.out
    assert "resid: 2" in captured.out

def test_removed_atoms():
    model = GenerateModel(PDB)
    model.set_region(name='test', selection='resid 3')

    residue = model.test.residues[0]
    removed_names = [a.name for a in residue.removed_atoms]
    assert removed_names == []

    for res in model.test.residues:
        assert res.removed_atoms == []

    model.test.truncate(scheme=QMzyme.TruncationSchemes.TerminalAlphaCarbon, selection="all")

    residue = model.test.residues[0]
    removed_names = [a.name for a in residue.removed_atoms]
    assert removed_names == ['N', 'H', 'C', 'O']

    for res in model.test.residues:
        removed = res.removed_atoms
        assert len(removed) == 4
        assert [a.name for a in removed] == ['N', 'H', 'C', 'O']
        assert all(a.resname == 'LEU' and a.resid == 3 for a in removed)
    
 
def test_added_atoms():
    model = GenerateModel(PDB)
    model.set_region(name='test', selection='resid 3')

    residue = model.test.residues[0]
    added_names = [a.name for a in residue.added_atoms]
    assert added_names == []

    for res in model.test.residues:
        assert res.added_atoms == []

    model.test.truncate(scheme=QMzyme.TruncationSchemes.TerminalAlphaCarbon, selection="all")

    residue = model.test.residues[0]
    added_names = [a.name for a in residue.added_atoms]
    assert added_names == ['HN', 'HC']

    for res in model.test.residues:
        added = res.added_atoms
        assert len(added) == 2
        assert [a.name for a in added] == ['HN', 'HC']
        assert all(a.resid == 3 for a in added)
 
def test_set_residue_method():
    # This test was made to chcek if the change in the attribute assignment for the method
    # from segid to method_type.
    model = GenerateModel(PDB)
    model.set_region(name='test', selection='resid 1-3')
    region = model.test
    original_segids = {atom.segid for atom in region.atoms}
    for residue in region.residues:
        assert residue.method_type != 'XTB'
 
    region.set_residue_method('XTB')
 
    for residue in region.residues:
        assert residue.method_type == 'XTB'
    assert {atom.segid for atom in region.atoms} == original_segids

def test_truncate():
    from QMzyme.TruncationSchemes import TerminalAlphaCarbon, AlphaCarbon, BetaCarbon

    model = GenerateModel(PDB)
    model.set_region(name='test', selection='resid 1-5')
    region = model.test
    region.set_residue_method('QM')
    n_atoms_before = region.n_atoms
 
    region.truncate(scheme=TerminalAlphaCarbon, selection='resid 1-5')
 
    # Truncation should not add atoms beyond what the original selection had
    assert model.test.n_atoms < n_atoms_before
    assert model.test.n_atoms == 79
    assert all(val == 'TerminalAlphaCarbon' for val in model.test._residue_truncation_attr.values())
    assert model.test._residue_capping_attr[5] == 'cap_H'
    assert all(val == 'QM' for val in model.test._residue_method_attr.values())

    model = GenerateModel(PDB)
    model.set_region(name='test', selection='resid 1-20')

    model.test.truncate(scheme=TerminalAlphaCarbon, selection='resid 1-5')
    model.test.truncate(scheme=AlphaCarbon, selection='resid 6-10', remove_ethane=False)
    model.test.truncate(scheme=BetaCarbon, selection='resid 11-15')

    region = model.test

    # TerminalAlphaCarbon: resid 1-5
    for resid in range(1, 6):
        res = region.get_residue(resid)
        assert res.removed_atoms == []
        assert res.truncation_params == 'TerminalAlphaCarbon'

    # AlphaCarbon: resid 6-10
    for resid in range(6, 11):
        res = region.get_residue(resid)
        assert res.truncation_params == 'AlphaCarbon'
        assert res.capping_scheme == 'cap_H'
        assert res.removed_atoms != []

    # BetaCarbon: resid 11-15
    res_11 = region.get_residue(11)
    assert res_11.truncation_params == 'BetaCarbon'
    assert res_11.removed_atoms == []
    assert res_11.added_atoms == []
    assert res_11.capping_scheme is None

    for resid in [12, 13, 14, 15]:
        res = region.get_residue(resid)
        assert res.truncation_params == 'BetaCarbon'

    for resid in [12, 13, 15]:
        res = region.get_residue(resid)
        assert res.removed_atoms != []
        assert res.capping_scheme == 'cap_H'

    res_14 = region.get_residue(14)
    assert res_14.removed_atoms == []
    assert res_14.capping_scheme is None

    # Untouched: resid 16-20
    for resid in range(16, 21):
        res = region.get_residue(resid)
        assert res.truncation_params is None
        assert res.removed_atoms == []
        assert res.added_atoms == []
        assert res.capping_scheme is None

def test_add_C_terminus_NME():
    model = GenerateModel(PDB)
    model.set_region(name='test', selection='resid 16')
    region = model.test
    n_before = region.n_atoms

    region.add_C_terminus_NME(16)

    # NME cap keeps neighbor's N, H (2 atoms) + adds CH3, HH31, HH32, HH33 (4 atoms)
    assert region.n_atoms == n_before + 6
    neighbor_atoms = [a for a in region.atoms if a.resid == 17]
    assert len(neighbor_atoms) == 6
    assert {a.name for a in neighbor_atoms} == {'N', 'H', 'CH3', 'HH31', 'HH32', 'HH33'}

    # Terminal C case. neighbor_resid doesn't exist in the universe at all
    model_cterm = GenerateModel(PDB)
    model_cterm.set_region(name='test', selection='resid 262')
    region_cterm = model_cterm.test
    n_before_cterm = region_cterm.n_atoms

    with pytest.warns(UserWarning, match="true chain terminus"):
        region_cterm.add_C_terminus_NME(262)

    # It should be skipped entirely
    assert region_cterm.n_atoms == n_before_cterm
    assert 263 not in {a.resid for a in region_cterm.atoms}
    assert 262 not in getattr(region_cterm, '_cap_flags', {})

    # Proline case. cap_NME raises a warning and nothing gets added
    model_pro = GenerateModel(PDB)
    model_pro.set_region(name='test', selection='resid 3')
    region_pro = model_pro.test
    n_before_pro = region_pro.n_atoms

    with pytest.warns(UserWarning) as record:
        region_pro.add_C_terminus_NME(3)
    messages = [str(w.message) for w in record]
    assert any("NME cap skipped for resid 3" in m for m in messages)
    assert any("Proline" in m for m in messages)
    assert region_pro.n_atoms == n_before_pro
    assert region_pro._cap_flags[3]['NME'] is False

    # Neighbor already present in the region but as a normal (uncapped) residue
    model_existing = GenerateModel(PDB)
    model_existing.set_region(name='test', selection='resid 10 or resid 11')
    region_existing = model_existing.test
    n_before_existing = region_existing.n_atoms

    region_existing.add_C_terminus_NME(10)

    assert region_existing.n_atoms == n_before_existing
    assert 10 not in getattr(region_existing, '_cap_flags', {})


def test_add_N_terminus_ACE():
    model = GenerateModel(PDB)
    model.set_region(name='test', selection='resid 40')
    region = model.test
    n_before = region.n_atoms

    region.add_N_terminus_ACE(40)

    # ACE cap keeps neighbor's C, O (2 atoms) + adds CH3, HH31, HH32, HH33 (4 atoms)
    assert region.n_atoms == n_before + 6
    neighbor_atoms = [a for a in region.atoms if a.resid == 39]
    assert len(neighbor_atoms) == 6
    assert {a.name for a in neighbor_atoms} == {'C', 'O', 'CH3', 'HH31', 'HH32', 'HH33'}

    # Terminal N case. neighbor_resid (resid - 1) doesn't exist in the universe
    model_nterm = GenerateModel(PDB)
    model_nterm.set_region(name='test', selection='resid 1')
    region_nterm = model_nterm.test
    n_before_nterm = region_nterm.n_atoms

    with pytest.warns(UserWarning, match="true chain terminus"):
        region_nterm.add_N_terminus_ACE(1)

    assert region_nterm.n_atoms == n_before_nterm
    assert 0 not in {a.resid for a in region_nterm.atoms}
    assert 1 not in getattr(region_nterm, '_cap_flags', {})

    # Neighbor already present in the region but as a normal (uncapped) residue
    model_existing = GenerateModel(PDB)
    model_existing.set_region(name='test', selection='resid 11 or resid 12')
    region_existing = model_existing.test
    n_before_existing = region_existing.n_atoms

    region_existing.add_N_terminus_ACE(12)

    assert region_existing.n_atoms == n_before_existing
    assert 12 not in getattr(region_existing, '_cap_flags', {})

    # Both caps attempted on the same resid but neighbor doesn't exist for one side
    model_lone = GenerateModel(PDB)
    model_lone.set_region(name='test', selection='resid 40')
    region_lone = model_lone.test
    region_lone.add_N_terminus_ACE(40)
    assert region_lone._cap_flags[40]['ACE'] is True
    assert 'NME' not in region_lone._cap_flags[40]


def test_insert_bridge_residue(capsys):
    model = GenerateModel(PDB)
    model.set_region(name='test', selection='resid 7 or resid 9')
    region = model.test

    # Cap resid 7's C-terminus; neighbor resid 8 becomes an NME-only cap fragment (6 atoms)
    region.add_C_terminus_NME(7)
    cap_only_atoms = [a for a in region.atoms if a.resid == 8]
    assert len(cap_only_atoms) == 6
    n_before_bridge = region.n_atoms

    # Cap resid 9's N-terminus; its neighbor (8) already holds NME atoms, which
    # triggers _insert_bridge_residue, replacing the 6 cap-only atoms with the real residue
    region.add_N_terminus_ACE(9)

    bridged_atoms = [a for a in region.atoms if a.resid == 8]
    assert 'CH3' not in [a.name for a in bridged_atoms]
    assert bridged_atoms[0].resname != 'NME'
    assert len(bridged_atoms) != 6
    assert region.n_atoms != n_before_bridge

    # Cap ACE first, then NME, exercising the mirror branch in add_C_terminus_NME (neighbor_resname == 'ACE') instead of add_N_terminus_ACE's
    model_rev = GenerateModel(PDB)
    model_rev.set_region(name='test', selection='resid 19 or resid 21')
    region_rev = model_rev.test

    region_rev.add_N_terminus_ACE(21)
    cap_only_rev = [a for a in region_rev.atoms if a.resid == 20]
    assert len(cap_only_rev) == 6
    assert cap_only_rev[0].resname == 'ACE'
    n_before_bridge_rev = region_rev.n_atoms

    region_rev.add_C_terminus_NME(19)
    bridged_rev = [a for a in region_rev.atoms if a.resid == 20]
    assert 'CH3' not in [a.name for a in bridged_rev]
    assert bridged_rev[0].resname != 'ACE'
    assert len(bridged_rev) != 6
    assert region_rev.n_atoms != n_before_bridge_rev

    # Bridge residue is Proline: should be kept whole, bit run through
    # BetaCarbon truncation, so ring atoms (CB/CG/CD) survive intact
    model_pro = GenerateModel(PDB)
    model_pro.set_region(name='test', selection='resid 3')
    region_pro = model_pro.test

    region_pro._insert_bridge_residue(4)  # PRO4, inserted directly as bridge
    pro_atoms = [a for a in region_pro.atoms if a.resid == 4]
    pro_atom_names = {a.name for a in pro_atoms}
    assert pro_atoms[0].resname == 'PRO'
    # full proline ring retained -> not reduced to backbone-only via BetaCarbon
    assert {'CB', 'CG', 'CD'}.issubset(pro_atom_names)

    pro_residue = next(r for r in region_pro.residues if r.resid == 4)
    assert pro_residue.truncation_params is None

    # Bridge residue triggered via a failed capping attempt (not a pre-existing
    # cap fragment): resid 3's NME cap fails because neighbor (PRO4) can't be capped,
    # then resid 3's ACE cap succeeds, so both flags are set and a bridge is attempted
    # on resid 3 itself rather than on the neighbor
    model_fail = GenerateModel(PDB)
    model_fail.set_region(name='test', selection='resid 3')
    region_fail = model_fail.test

    region_fail.add_C_terminus_NME(3)  # fails (Proline neighbor), flags NME=False
    capsys.readouterr()  # clear captured output from the NME failure
    region_fail.add_N_terminus_ACE(3)  # succeeds, flags ACE=True, both flags now set

    assert region_fail._cap_flags[3]['NME'] is False
    assert region_fail._cap_flags[3]['ACE'] is True
    # bridge attempt fired on resid 3 (no-op here since resid 3 isn't cap-only,
    # but confirms no exception was raised and the real residue atoms remain)
    resid3_atoms = [a for a in region_fail.atoms if a.resid == 3]
    assert resid3_atoms[0].resname == 'LEU'