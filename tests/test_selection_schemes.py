"""
Tests for the QMzyme GenerateModel.py code.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow

from QMzyme.GenerateModel import GenerateModel
import pytest
from QMzyme.RegionBuilder import RegionBuilder
from MDAnalysis.core.universe import Universe
from QMzyme.SelectionSchemes import DistanceCutoff
from QMzyme.data import PDB


def test_DistanceCutoff():
    from QMzyme.SelectionSchemes import DistanceCutoff
    model = GenerateModel(PDB)
    with pytest.raises(UserWarning):
        DistanceCutoff(model=model, name=None, cutoff=3, include_whole_residues=True)
    with pytest.raises(UserWarning):
        model.set_region(selection=DistanceCutoff, name=None, cutoff=3, include_whole_residues=True)
    model.set_catalytic_center('resid 263')
    model.set_region(selection=DistanceCutoff, name=None, cutoff=3, include_whole_residues=True)
    assert len(model.regions) == 2
    assert model.cutoff_3.n_atoms == 275
    assert model.cutoff_3.n_residues == 19

    model.set_region(selection=DistanceCutoff, name=None, cutoff=5, include_whole_residues=True)
    assert len(model.regions) == 3
    assert model.cutoff_5.n_atoms == 427
    assert model.cutoff_5.n_residues == 33
