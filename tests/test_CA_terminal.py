
from QMzyme.GenerateModel import GenerateModel
from QMzyme.TruncationSchemes import *
from QMzyme.truncation_utils import *
import pytest
from QMzyme.data import PDB, protein_residues

def test_CA_terminal():
    model = GenerateModel(PDB)
    model.set_region(name='region', selection='resid 263 or resid 16 or resid 17 or resid 57')
    truncation = TerminalAlphaCarbon(region=model.region, name=None)
    truncate = truncation.return_region()
    truncate.write('truncated_test')


    model.set_region(name='region1', selection='resid 263 or resid 16')
    model.set_region(name='region2', selection='resid 17 or resid 57')
    combined = model.region1.combine(model.region2)
    truncation = TerminalAlphaCarbon(region=combined, name=None)
    truncated_combined = truncation.return_region()
    truncated_combined.write('truncated_combined_test')

    assert truncate.n_atoms == truncated_combined.n_atoms
    for atom in truncate.atoms:
        assert atom.name == truncated_combined.get_atom(atom.id).name

