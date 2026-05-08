"""
Tests for the QMzyme RegionBuilder.py ands QMzymeRegion.py code.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow


import numpy as np
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
    model.set_region(name='test_reg', selection='resid 1')
    res = model.test_reg.residues[0] 

    assert "MET" in repr(res)
    res.set_chain("A")
    assert res.get_atom("CA").name == "CA"
    
    # Touch backbone and charge guessing
    assert len(res.get_backbone_atoms()) >= 4 
    res.guess_charge(verbose=False)
    assert res.charge is not None