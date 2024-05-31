"""
Tests for the QMzymeAtom.py code.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow

import pytest
from QMzyme.QMzymeAtom import QMzymeAtom


atom_dict = {
    "name": "CA",
    "element": "C",
    "position": [1.000, 2.000, 3.000],
    "resid": "1",
    "resname": "VAL",
    "id": 100,
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
        assert qmz_atom.random_property == atom_input["random_property"]
        assert str(qmz_atom.__repr__()) == f'<QMzymeAtom {atom_input["id"]}: {atom_input["name"]} of resname {atom_input["resname"]}, resid {atom_input["resid"]}>'
                
def test_setters():
    qmz_atom = QMzymeAtom(**atom_dict)

    # name
    # with pytest.raises(AttributeError):
    #     qmz_atom.name = 'F'
    assert qmz_atom.name != 'F'
    #qmz_atom.set_name('F')
    qmz_atom.name = 'F'
    assert qmz_atom.name == 'F'

    # element
    # with pytest.raises(AttributeError):
    #     qmz_atom.element = 'P'
    assert qmz_atom.element != 'P'
    #qmz_atom.set_element('P')
    qmz_atom.element = 'P'
    assert qmz_atom.element == 'P'

    # region
    assert qmz_atom.region is None
    with pytest.raises(AttributeError):
        qmz_atom.region = 'placeholder'
    assert qmz_atom.region != 'placeholder'
    qmz_atom._set_region('placeholder')
    assert qmz_atom.region == 'placeholder'

    # is_fixed
    assert not qmz_atom.is_fixed
    qmz_atom.set_fixed()
    assert qmz_atom.is_fixed

    # is_point_charge
    assert qmz_atom.is_point_charge == False
    with pytest.raises(UserWarning):
        qmz_atom.set_point_charge()
    #assert qmz_atom.is_point_charge

    # is_neighbor
    with pytest.raises(AttributeError):
        f = qmz_atom.is_neighbor
    qmz_atom.set_neighbor()
    assert qmz_atom.is_neighbor


