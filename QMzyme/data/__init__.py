from importlib import resources
import copy
from QMzyme.configuration import protein_residues, residue_charges, backbone_atoms

_data_ref = resources.files('QMzyme.data')

PDB = (_data_ref / '1oh0.pdb').as_posix()
TOP = (_data_ref / '1oh0_equ.prmtop').as_posix()
RST = (_data_ref / '1oh0_equ.rst7').as_posix()
DCD = (_data_ref / '1oh0_equ.prod_1.stripped.dcd').as_posix()
PQR = (_data_ref / '1oh0_equ.prod_1.stripped.pqr').as_posix()

# protein_residues = QMzyme.configuration.protein_residues
# residue_charges = QMzyme.configuration.residue_charges
# backbone_atoms = QMzyme.configuration.backbone_atoms
