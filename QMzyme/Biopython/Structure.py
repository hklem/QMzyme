# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""The structure class, representing a macromolecular structure."""

from QMzyme.Biopython.Entity import Entity

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from QMzyme.Biopython.Model import Model


class Structure(Entity[None, "Model"]):
    """The Structure class contains a collection of Model instances."""

    def __init__(self, id):
        """Initialize the class."""
        self.level = "S"
        Entity.__init__(self, id)

    def __repr__(self):
        """Return the structure identifier."""
        return f"<Structure id={self.get_id()}>"

    def get_models(self):
        """Return models."""
        yield from self

    def get_chains(self):
        """Return chains from models."""
        for m in self.get_models():
            yield from m

    def get_residues(self):
        """Return residues from chains."""
        for c in self.get_chains():
            yield from c

    def get_atoms(self):
        """Return atoms from residue."""
        for r in self.get_residues():
            yield from r
