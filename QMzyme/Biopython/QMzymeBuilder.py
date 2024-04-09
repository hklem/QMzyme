# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Consumer class that builds a Structure object.

This is used by the PDBParser and MMCIFparser classes.
"""

from typing import Optional
import numpy as np
import warnings

# SMCRA hierarchy
from QMzyme.Biopython.Structure import Structure
from QMzyme.Biopython.Model import Model
from QMzyme.Biopython.Chain import Chain
from QMzyme.Biopython.Residue import Residue, DisorderedResidue
from QMzyme.Biopython.PDBExceptions import PDBConstructionException
from QMzyme.Biopython.PDBExceptions import PDBConstructionWarning
from QMzyme.Biopython.StructureBuilder import StructureBuilder

class QMzymeBuilder(StructureBuilder):
    """Modified class to constructe QMzyme Structure object.

    The QMzymeBuilder class is used by the PDBParser classes to
    translate a file to a Structure object to be used in QMzyme.
    """

    def __init__(self):
        """Initialize this instance."""
        StructureBuilder.__init__(self)

    def init_model(self, model_id: int, serial_num: Optional[int] = None):
        """Create a new Model object with given id.

        Arguments:
         - id - int
         - serial_num - int
        """
        self.model = Model(model_id, serial_num)
        setattr(self.structure, 'parent', self.model)
        #self.structure.add(self.model)

    def init_chain(self, chain_id: str):
        """Create a new Chain object with given id.

        Arguments:
         - chain_id - string
        """
        if self.parent.has_id(chain_id):
            self.chain = self.parent[chain_id]
        else:
            self.chain = Chain(chain_id)
            self.parent.add(self.chain)

    def init_residue(self, resname, field, resseq, icode):
        """Create a new Residue object.

        Arguments:
         - resname - string, e.g. "ASN"
         - field - hetero flag, "W" for waters, "H" for
           hetero residues, otherwise blank.
         - resseq - int, sequence identifier
         - icode - string, insertion code

        """
        if field != " ":
            if field == "H":
                # The hetero field consists of H_ + the residue name (e.g. H_FUC)
                #field = "H_" + resname
                field = "HETATM"
        #res_id = (field, resseq, icode)
        res_id = (resname, resseq, self.chain.id)

        if field == " ":
            if self.chain.has_id(res_id):
                # There already is a residue with the id (field, resseq, icode).
                # This only makes sense in the case of a point mutation.
                warnings.warn(
                    "WARNING: Residue ('%s', %i, '%s') redefined at line %i."
                    % (resname, resseq, self.chain.id, self.line_counter),
                    PDBConstructionWarning,
                )
                duplicate_residue = self.chain[res_id]
                if duplicate_residue.is_disordered() == 2:
                    # The residue in the chain is a DisorderedResidue object.
                    # So just add the last Residue object.
                    if duplicate_residue.disordered_has_id(resname):
                        # The residue was already made
                        self.residue = duplicate_residue
                        duplicate_residue.disordered_select(resname)
                    else:
                        # Make a new residue and add it to the already
                        # present DisorderedResidue
                        new_residue = Residue(res_id, resname, self.segid)
                        duplicate_residue.disordered_add(new_residue)
                        self.residue = duplicate_residue
                        return
                else:
                    if resname == duplicate_residue.resname:
                        warnings.warn(
                            "WARNING: Residue ('%s', %i, '%s') already defined "
                            "with the same name at line  %i."
                            % (resname, resseq, self.chain.id, self.line_counter),
                            PDBConstructionWarning,
                        )
                        self.residue = duplicate_residue
                        return
                    # Make a new DisorderedResidue object and put all
                    # the Residue objects with the id (field, resseq, icode) in it.
                    # These residues each should have non-blank altlocs for all their atoms.
                    # If not, the PDB file probably contains an error.
                    if not _is_completely_disordered(duplicate_residue):
                        # if this exception is ignored, a residue will be missing
                        self.residue = None
                        raise PDBConstructionException(
                            "Blank altlocs in duplicate residue %s ('%s', %i, '%s')"
                            % (resname, field, resseq, icode)
                        )
                    self.chain.detach_child(res_id)
                    new_residue = Residue(res_id, resname, self.segid)
                    disordered_residue = DisorderedResidue(res_id)
                    self.chain.add(disordered_residue)
                    disordered_residue.disordered_add(duplicate_residue)
                    disordered_residue.disordered_add(new_residue)
                    self.residue = disordered_residue
                    return
        self.residue = Residue(res_id, resname, self.segid)
        setattr(self.residue, 'resid', resseq)
        delattr(self.residue, xtra)
        self.chain.add(self.residue)