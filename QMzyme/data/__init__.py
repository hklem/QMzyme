from importlib import resources
import copy

_data_ref = resources.files('QMzyme.data')

PDB = (_data_ref / '1oh0.pdb').as_posix()
TOP = (_data_ref / '1oh0_equ.prmtop').as_posix()
RST = (_data_ref / '1oh0_equ.rst7').as_posix()
DCD = (_data_ref / '1oh0_equ.prod_1.stripped.dcd').as_posix()
PQR = (_data_ref / '1oh0_equ.prod_1.stripped.pqr').as_posix()


protein_residues = {
    'ALA': 0,
    'ARG': 1,
    'ASH': 0,
    'ASN': 0,
    'ASP': -1,
    'CYM': -1,
    'CYS': 0,
    'CYX': 0,
    'GLH': 0,
    'GLN': 0,
    'GLU': -1,
    'GLY': 0,
    'HID': 0,
    'HIE': 0,
    'HIP': 1,
    'HYP': 0,
    'ILE': 0,
    'LEU': 0,
    'LYN': 0,
    'LYS': 1,
    'MET': 0,
    'PHE': 0,
    'PRO': 0,
    'SER': 0,
    'THR': 0,
    'TRP': 0,
    'TYR': 0,
    'VAL': 0,
 }

residue_charges = copy.copy(protein_residues)
residue_charges['WAT'] = 0
