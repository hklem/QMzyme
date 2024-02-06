###############################################################################
# Code written by Heidi Klem while at
# Colorado State University as a graduate student
# in the Paton and McCullagh groups and at the
# National Institute of Standards and Technology
# as an NRC Postdoc (Fed).
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

import os
import numpy as np
from QMzyme.utils import get_outlines
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
solvent_list=['HOH','WAT','T3P','SOL']
protein_residues =  ['ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'CYM', 'CYS', 'CYX',
                     'GLH', 'GLN', 'GLU', 'GLY', 'HIS', 'HID', 'HIE', 'HIP',
                     'HYP', 'ILE', 'LEU', 'LYN', 'LYS', 'MET', 'PHE', 'PRO',
                     'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HSE', 'HSD', 'HSP' ]

def store_mol_pdb(mol):
    Chem.MolToPDBFile(mol,'temp2.pdb')
    data = get_outlines('temp2.pdb')
    os.remove('temp2.pdb')
    return data

def atom_coords(mol, atom):
    '''
    Parameters
    ----------
    mol : rdkit mol object.
    atom : atom from rdkit mol object.

    Returns
    -------
    numpy array of the atomic cartesian coorinates.
    '''
    return np.asarray(mol.GetConformer().GetAtomPosition(atom.GetIdx()))

def mol_from_pdb(file):
    return Chem.MolFromPDBFile(file, removeHs=False, sanitize=False)

def centroid_coords(mol):
    return np.asarray(Chem.rdMolTransforms.ComputeCentroid((mol.GetConformer())))

def get_atoms(mol):
    return mol.GetAtoms()

def check_pdb_rdkit(file):
    '''
    Function to assess PDB format, identify any formatting incompatibilities
    with the QMzyme code, and gather useful basic information.

    Parameters
    ----------
    file : str
        PDB file.

    Returns
    -------
    protein_res_count : TYPE
        DESCRIPTION.
    h_present : TYPE
        DESCRIPTION.
    non_protein_residues : TYPE
        DESCRIPTION.
    non_protein_res_count : TYPE
        DESCRIPTION.
    non_protein_chemical_name : TYPE
        DESCRIPTION.

    '''
    wat_count, protein_res_count, non_protein_res_count = 0,0,0
    data, protein_seq, non_protein_seq = [], [], []
    h_present, no_chain_info = False, False
    previous_res = None
    edits_made = False
    residue_count = 0
    res_num_previous = None
    non_protein_residues = {}
    non_protein_residues['Chain'] = []
    non_protein_residues['Name'] = []
    non_protein_residues['Number'] = []
    non_protein_chemical_name = {}

    mol = Chem.MolFromPDBFile(file,
                              removeHs=False,
                              sanitize=False)

    for atom in mol.GetAtoms():
        if rdkit_info(atom,'Atom name').startswith('H'):
            h_present = True
        if rdkit_info(atom,'Chain') == ' ':
            no_chain_info = True
            atom.GetPDBResidueInfo().SetChainId('X')
        current_res = define_residue(atom)
        if previous_res == current_res:
            continue
        previous_res = current_res
        if current_res['Residue name'] in solvent_list:
            wat_count += 1
            continue
        if atom.GetPDBResidueInfo().GetIsHeteroAtom() is True:
            non_protein_seq.append(current_res)
            continue
        if current_res['Residue name'] in protein_residues:
            protein_seq.append(current_res)
            protein_res_count += 1
        else:
            non_protein_seq.append(current_res)

    print("Information in PDB file: {}".format(file))
    if h_present == True:
        print("Hydrogens are present")
    else:
        print("Hydrogens are not present.")
    if no_chain_info == True:
        print("Chain IDs not defined were set to 'X'.")
    print("Total number of atoms: {}".format(mol.GetNumAtoms()))
    print("Water molecules: {}".format(wat_count))
    print("Standard amino acid residues: {}".format(protein_res_count))
    print("The following {} non-protein residues were detected:"
          .format(non_protein_res_count))
    if non_protein_res_count > 0:
        chain = non_protein_residues['Chain']
        name = non_protein_residues['Name']
        number = non_protein_residues['Number']
        for i in range(non_protein_res_count):
            chemical_name = non_protein_chemical_name[name[i]]
            print("Chain: {}".format(chain[i]),
                  " Residue Name: {}".format(name[i]),
                  " Chemical Name: {}".format(chemical_name),
                  " Residue Number: {}".format(number[i]))

    return protein_res_count, h_present, non_protein_residues, non_protein_res_count, non_protein_chemical_name

def fix_h_bond(mol, h_atom, length=1.00):
    rdMolTransforms.SetBondLength(mol.GetConformer(),
                                  h_atom.GetNeighbors()[0].GetIdx(),
                                  h_atom.GetIdx(), length)
    return mol

def h_cap(new_mol, atom_id, cap_name=' H* '):
    cap_atom = new_mol.GetAtomWithIdx(atom_id)
    cap_atom.SetAtomicNum(1)
    cap_atom.GetPDBResidueInfo().SetName(cap_name)
    return new_mol

def remove_atoms(mol, indices):
    indices = list(set(indices))
    indices.sort(reverse=True)
    for i in indices:
        mol.RemoveAtom(i)
    return mol

def rdkit_info(atom,info='Atom name'):
    if type(atom) is list:
        return [rdkit_dispatch[info](at) for at in atom]
    else:
        return rdkit_dispatch[info](atom)


def define_residue(atom):
    res = {}
    for key in rdkit_dispatch:
        if key != 'Atom name':
            res[key] = rdkit_dispatch[key](atom)
    return res

rdkit_dispatch = {
    'Residue name': lambda atom: atom.GetPDBResidueInfo().GetResidueName().split()[0],
    'Atom name': lambda atom: atom.GetPDBResidueInfo().GetName().split()[0],
    'Residue number': lambda atom: atom.GetPDBResidueInfo().GetResidueNumber(),
    'Chain': lambda atom: check_chain(atom)
    }

def check_chain(atom):
    try:
        chain = atom.GetPDBResidueInfo().GetChainId().split()[0]
    except:
        chain = ' '
    return chain

def show_mol(mol):
    mol = Chem.RWMol(mol)
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Chem import rdDepictor
    rdDepictor.SetPreferCoordGen(True)
    from IPython.display import SVG
    from rdkit.Chem import rdCoordGen
    rdCoordGen.AddCoords(mol)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'C':
            atom.SetProp("atomLabel",atom.GetSymbol())
    drawer = rdMolDraw2D.MolDraw2DSVG(400,400)
    drawer.drawOptions().addStereoAnnotation = False
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    return SVG(drawer.GetDrawingText())
