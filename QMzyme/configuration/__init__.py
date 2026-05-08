###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
Module containing the following information used by QMzyme that users can 
modify to their specific needs. 
    - Protein residue three letter names and associated formal charges 
    - Any residue three letter names and associated formal charges 
    - Single letter atom names of amino-acid backbone atoms
"""

import copy


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
    'HIS': 0,
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
    'VAL': 0
    }
"""
The entries in protein_residues are used in TruncationSchemes.py to ignore non
amino acid residues, since those will not have backbone atoms and would raise
errors during truncation. If you want to add a residue and charge but the residue
is not an amino-acid, or should not be considered in truncation procedures, please 
add that residue to the residue_charges dictionary, rather than the protein_residues
dictionary. All residues in the protein_residues dictionary are automatically copied
into the residue_charges dictionary, but not vice-versa.
"""


residue_charges = copy.copy(protein_residues)
"""
If there is a non-native amino acid you will encounter often you can add its 
name and charge to the protein_residues dictionary so you do not have 
to manually add that information each time you run QMzyme. 
"""
residue_charges['WAT'] = 0
residue_charges['HOH'] = 0
residue_charges['Na+'] = 1
residue_charges['Cl-'] = -1


backbone_atoms = {
    'N':'N', 
    'H':'H', 
    'CA':'CA', 
    'HA':'HA', 
    'C':'C', 
    'O':'O'
    }
"""
The entries in backbone_atoms are used in TruncationSchemes.py to decide what atoms to
remove and replace with hydrogen. The keys should not be changed in this dictionary, only
the values assigned to each key. The default backbone_atoms dictionary is designed to match 
the conventions used by the `chemical component dictionary <http://www.wwpdb.org/data/ccd>`_.


If your amino acid backbone atom names do not match the chemical component dictionary 
convention ('N', 'H', 'CA', 'HA', 'C', 'O') you can edit that here globally, or you can
redefine them at the start of your QMzyme run by re-defining the QMzyme.configuration.backbone_atoms
dictionary accordingly. This might be necessary if you are using a structure generated from a 
force-field topology that, for example, names the H atom bound to the backbone N atom 'HN' instead
of 'H'. In this case, you could simply set QMzyme.configuration.backbone_atoms['H'] = 'HN'. 
"""

element_name_to_atomic_number = {
    'H': 1,
    'He': 2,
    'Li': 3,
    'Be': 4,
    'B': 5,
    'C': 6,
    'N': 7,
    'O': 8,
    'F': 9,
    'Ne': 10,
    'Na': 11,
    'Mg': 12,
    'Al': 13,
    'Si': 14,
    'P': 15,
    'S': 16,
    'Cl': 17,
    'Ar': 18,
    'K': 19,
    'Ca': 20,
    'Sc': 21,
    'Ti': 22,
    'V': 23,
    'Cr': 24,
    'Mn': 25,
    'Fe': 26,
    'Co': 27,
    'Ni': 28,
    'Cu': 29,
    'Zn': 30,
    'Ga': 31,
    'Ge': 32,
    'As': 33,
    'Se': 34,
    'Br': 35,
    'Kr': 36,
    'Rb': 37,
    'Sr': 38,
    'Y': 39,
    'Zr': 40,
    'Nb': 41,
    'Mo': 42,
    'Tc': 43,
    'Ru': 44,
    'Rh': 45,
    'Pd': 46,
    'Ag': 47,
    'Cd': 48,
    'In': 49,
    'Sn': 50,
    'Sb': 51,
    'Te': 52,
    'I': 53,
    'Xe': 54
}

"""
This module is used to convert the element name to the atomic number. This is especially useful
when trying to use in conjunction with cclib, due to cclib only providing atomic number.
"""