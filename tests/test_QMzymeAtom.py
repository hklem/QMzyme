"""
Tests for the QMzyme QMzymeAtom.py code.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow

import sys
import shutil
import os
import pytest
from QMzyme.QMzymeAtom import QMzymeAtom
from MDAnalysis.core.groups import Atom
import MDAnalysis as MDA
from importlib_resources import files

pdb_file = str(files('QMzyme.data').joinpath('1oh0.pdb'))

atom_dict = {
    "name": "CA",
    "element": "C",
    "position": [1.000, 2.000, 3.000],
    "resid": "1",
    "resname": "VAL",
    "id": 100,
    "region": "placeholder",
    "random_property": "placeholder",
}


@pytest.mark.parametrize(
        "test_type, atom_input",
        [("Testing required args", atom_dict), ("Testing optional args", atom_dict)]
)
def test_init(test_type, atom_input):
    if test_type == "Testing required args":
        name = atom_input["name"]
        element = atom_input["element"]
        position = atom_input["position"]
        resid = atom_input["resid"]
        resname = atom_input["resname"]
        qmz_atom = QMzymeAtom(name, element, position, resid, resname)

        assert qmz_atom.name == name
        assert qmz_atom.element == element
        assert qmz_atom.position == position
        assert qmz_atom.resid == resid
        assert qmz_atom.resname == resname
        assert str(qmz_atom.__repr__()) == f'<QMzymeAtom {qmz_atom.id}: {name} of resname {resname}, resid {resid}>'

    if test_type == "Testing optional args":
        qmz_atom = QMzymeAtom(**atom_dict)

        assert qmz_atom.name == atom_input["name"]
        assert qmz_atom.element == atom_input["element"]
        assert qmz_atom.position == atom_input["position"]
        assert qmz_atom.resid == atom_input["resid"]
        assert qmz_atom.resname == atom_input["resname"]
        assert qmz_atom.id == atom_input["id"]
        assert qmz_atom.region == atom_input["region"]
        assert qmz_atom.random_property == atom_input["random_property"]
        assert str(qmz_atom.__repr__()) == f'<QMzymeAtom {atom_input["id"]}: {atom_input["name"]} of resname {atom_input["resname"]}, resid {atom_input["resid"]}>'
                
def test_set_neighbor():
    qmz_atom = QMzymeAtom(**atom_dict)
    assert hasattr(qmz_atom, "is_neighbor") == False
    qmz_atom.set_neighbor(True)
    assert hasattr(qmz_atom, "is_neighbor") == True
    assert qmz_atom.is_neighbor == True
