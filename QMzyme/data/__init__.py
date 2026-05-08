from importlib import resources
import copy
from QMzyme.configuration import protein_residues, residue_charges, backbone_atoms, element_name_to_atomic_number

_data_ref = resources.files('QMzyme.data')

PDB = (_data_ref / '1oh0.pdb').as_posix()
TOP = (_data_ref / '1oh0_equ.prmtop').as_posix()
RST = (_data_ref / '1oh0_equ.rst7').as_posix()
DCD = (_data_ref / '1oh0_equ.prod_1.stripped.dcd').as_posix()
PQR = (_data_ref / '1oh0_equ.prod_1.stripped.pqr').as_posix()
CSA_holo = (_data_ref / 'CSA_holo.log').as_posix()
CSA_apo = (_data_ref / 'CSA_apo.log').as_posix()
CSA_pkl = (_data_ref / 'CSA_holo_apo.pkl').as_posix()
Cutoff_3 = (_data_ref / '1OH0_cutoff3.pdb').as_posix()

# protein_residues = QMzyme.configuration.protein_residues
# residue_charges = QMzyme.configuration.residue_charges
# backbone_atoms = QMzyme.configuration.backbone_atoms
