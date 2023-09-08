################################################
# Code written by Heidi Klem while at
# Colorado State University as a graduate student
# in the Paton and McCullagh groups and at the
# National Institute of Standards and Technology
# as an NRC Postdoc (Fed).
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
#################################################

"""Generate QM-based enzyme model."""


def canvas(with_attribution=True):
	"""
	Placeholder function to show example docstring (NumPy format).

	Replace this function and doc string for your own project.

	Parameters
	----------
	with_attribution : bool, Optional, default: True
		Set whether or not to display who the quote is from.

	Returns
	-------
	quote : str
		Compiled string including quote and optional attribution.
	"""

	quote = "The code is but a canvas to our imagination."
	if with_attribution:
		quote += "\n\t- Adapted from Henry David Thoreau"
	return quote

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
#from openbabel import openbabel as ob
from urllib.request import urlopen
from rdkit.Chem import rdMolTransforms

protein_residues =	["ALA", "ARG", "ASH", "ASN", "ASP", "CYM",
 "CYS", "CYX", "GLH", "GLN", "GLU", "GLY", "HIS", "HID", "HIE",
 "HIP", "HYP", "ILE", "LEU", "LYN", "LYS", "MET", "PHE", "PRO",
 "SER", "THR", "TRP", "TYR", "VAL"]

def res_info(atom,info):
	''' Options for info are 'chain', 'res_name', 'res_number', 'atom_name'.
	'''
	if info == 'res_name':
		return atom.GetPDBResidueInfo().GetResidueName()
	if info == 'atom_name':
		return atom.GetPDBResidueInfo().GetName()
	if info == 'res_number':
		return atom.GetPDBResidueInfo().GetResidueNumber()
	if info == 'chain':
		return atom.GetPDBResidueInfo().GetChainId()

def define_residue(atom):
	''' Returns a tuple that completely defines the residue that atom belongs to: ('chain','residue name', 'residue number')
	'''
	return (res_info(atom,'chain'),res_info(atom,'res_name'),res_info(atom,'res_number'))


def atom_coords(mol,atom):
	return np.asarray(mol.GetConformer().GetAtomPosition(atom.GetIdx()))

def download(pdb_list):
	baseUrl = 'http://www.pdb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId='
	for structure in pdb_list.split():
		pdbUrl = baseUrl + structure[:4]
		outFileName = structure[:4] + '.pdb'
		print(pdbUrl)

		with urlopen(pdbUrl) as response, open(outFileName, 'wb') as outFile:

			data = response.read()
			outFile.write(data)
			print('Downloading {} as {}.'.format(structure[:4], outFileName))

	return(data)

def PDB_file_info(pdb_file):
	'''Input: PDB file name. Prints helpful information about what is present in the structure.
	'''
	base_mol = Chem.MolFromPDBFile(pdb_file,removeHs=False,sanitize=False)
	wat_count, protein_res_count, non_protein_res_count, anion_count, cation_count= 0, 0, 0, 0, 0
	protein_seq, non_protein_seq = [], []
	h_present, no_chain_info = False, False
	previous_res = None

	for atom in base_mol.GetAtoms():
		if res_info(atom,'atom_name') == ' H  ':
			h_present = True
		if res_info(atom,'chain') == ' ':
			no_chain_info = True
			atom.GetPDBResidueInfo().SetChainId('X')
		current_res = (res_info(atom,'chain'),res_info(atom,'res_name'),res_info(atom,'res_number'))
		if previous_res == current_res:
			continue
		previous_res = current_res
		if res_info(atom,'res_name') in ["WAT","HOH"]:
			wat_count += 1
			continue
		if "+" in res_info(atom,'res_name'):
			cation_count += 1
			cation = res_info(atom,'res_name')
			continue
		if "-" in res_info(atom,'res_name'):
			anion_count += 1
			anion = res_info(atom,'res_name')
			continue
		if atom.GetPDBResidueInfo().GetIsHeteroAtom() is True:
			non_protein_seq.append(current_res)
			non_protein_res_count += 1
			continue
		if res_info(atom,'res_name') in protein_residues:
			protein_seq.append(current_res)
			protein_res_count += 1
		else:
			non_protein_seq.append(current_res)
			non_protein_res_count += 1

	print("Information on PDB file: {}".format(pdb_file))
	if h_present == True:
		print("Hydrogens are present")
	else:
		print("Hydrogens are not present. Structure will be automatically protonated during truncation step.")
	if no_chain_info == True:
		print("Chain IDs not defined were set to 'X'.")
	print("----------------------------------------------")
	print("Total number of atoms: {}".format(base_mol.GetNumAtoms()))
	print("Water molecules: {}".format(wat_count))
	if cation_count > 0:
		print("{} ions: {}".format(cation,cation_count))
	if cation_count == 0:
		if anion_count == 0:
			print("Ions: 0")
	if anion_count > 0:
		print("{} ions: {}".format(anion,anion_count))
	print("Standard amino acid residues: {}".format(protein_res_count))
	print("----------------------------------------------")
	print("The following {} ligands were detected:".format(non_protein_res_count))
	if non_protein_res_count > 0:
		for i in range(non_protein_res_count):
			print("{}:".format(i+1))
			print("Chain: {}".format(non_protein_seq[i][0])," Name: {}".format(non_protein_seq[i][1])," ID: {}".format(non_protein_seq[i][2]))
			print(".........")
	print("----------------------------------------------")

def separate_protein(pdb_file,add_protein_residue=[]):
	""" Useful to identify non-amino acid molecules in the pdb file. This can be used to identify bound ligands, substrates, cofactors, etc. The code automatically removes water molecules for simplicity. You can add any three letter residue name to be considered as protein. For example, you might want to include water molecules in the protein mol object, so you would include the argument add_protein_residue=['WAT'] or add_protein_residue=['HOH'] depending on how it is defined in the PDB. This function automatically generates structure files in pdb format of the separated protein and non_protein residues."""

	if len(add_protein_residue)>0:
		for i in range(len(add_protein_residue)):
			protein_residues.append(add_protein_residue[i])

	base_mol = Chem.MolFromPDBFile(pdb_file,removeHs=False,sanitize=False)
	protein_mol = Chem.RWMol(base_mol)
	non_protein_mol = Chem.RWMol(base_mol)
	non_protein_id = []
	protein_id = []
	non_protein_residues = []

	for atom in base_mol.GetAtoms():
		if res_info(atom,'res_name') in protein_residues:
			protein_id.append(atom.GetIdx())
			continue
		elif res_info(atom,'res_name') == "WAT":
			protein_id.append(atom.GetIdx())
			non_protein_id.append(atom.GetIdx())
		elif res_info(atom,'res_name') == "HOH":
			protein_id.append(atom.GetIdx())
			non_protein_id.append(atom.GetIdx())
		else:
			x = '{}{}'.format(res_info(atom,'res_name'),res_info(atom,'res_number'))
			non_protein_residues.append(x)
			non_protein_id.append(atom.GetIdx())

	non_protein_residues = list(set(non_protein_residues))
	protein_id.sort(reverse=True)
	non_protein_id.sort(reverse=True)

	for atom in protein_id:
		non_protein_mol.RemoveAtom(atom)
	for atom in non_protein_id:
		protein_mol.RemoveAtom(atom)

	if len(non_protein_residues) == 0:
		print('No non protein residues found in PDB.')

	if len(non_protein_residues) > 0:
		Chem.MolToPDBFile(non_protein_mol,pdb_file.split('.')[0]+'_non_protein_residues'+'.pdb')
		print('Non protein residue(s) found in PDB:')
		for res in non_protein_residues:
			print(res)

	Chem.MolToPDBFile(protein_mol,pdb_file.split('.')[0]+'_protein_residues'+'.pdb')
	return protein_mol, non_protein_mol

def catalytic_center(pdb_file, catalytic_center=[],definition=['res_number','res_name','chain']):
	if len(catalytic_center) != len(definition):
		print("Catalytic center is defined by {} components: {} but {} were given: {}".format(len(definition),definition,len(catalytic_center),catalytic_center))
		print("Either change definition or change catalytic_center.")
	else:
		protein_mol = Chem.MolFromPDBFile(pdb_file,removeHs=False,sanitize=False)
		catalytic_center_mol = Chem.RWMol(protein_mol)
		count = 0
		remove_ids = []
		for atom in reversed(protein_mol.GetAtoms()):
			for component in definition:
				if res_info(atom,component) not in catalytic_center:
					if atom.GetIdx() not in remove_ids:
						remove_ids.append(atom.GetIdx())
				else:
					count +=1
		if count == 0:
			print("WARNING: No atoms found matching catalytic center definition.")
		else:
			remove_ids.sort(reverse=True)
			for atom in remove_ids:
				catalytic_center_mol.RemoveAtom(atom)
			print('Catalytic center contains {} atoms.'.format(catalytic_center_mol.GetNumAtoms()))
			print("Structure saved as catalytic_center.pdb")
			Chem.MolToPDBFile(catalytic_center_mol,'catalytic_center.pdb')
			return catalytic_center_mol, protein_mol

def residue_shell(center_mol,radius,pdb_file=None,base_mol=None,centroid=False,include_residues=[]):
	if base_mol == None:
		if pdb_file == None:
			print("Error: Must define on of the following: pdb_file or base_mol")
		base_mol = Chem.MolFromPDBFile(pdb_file,removeHs=False,sanitize=False)

	res_name, res_number, res_chain, res_atom, atom_type = [], [], [], [], []
	N_termini_interactions = []
	side_chain_interactions = []
	C_termini_interactions = []
	current_res = 'nothing'
	backbone_atoms=[' O  ',' C	',' N  ',' CA ']

	if centroid == True:
		centroid_coords = np.asarray(Chem.rdMolTransforms.\
			ComputeCentroid((center_mol.GetConformer())))
		distances = [np.linalg.norm(np.asarray(center_mol.GetConformer().\
			GetAtomPosition(atom.GetIdx()))-centroid_coords) for atom in \
			center_mol.GetAtoms()]
		radius_buffer = np.max(distances)

		for atom1 in base_mol.GetAtoms():
			if res_info(atom1,'atom_name') == ' H  ':
				continue
			if str(res_info(atom1,'res_name'))+str(res_info(atom1,'res_number')) in include_residues:
				atomic_distance = radius
			else:
				coords1 = np.asarray(base_mol.GetConformer().GetAtomPosition(atom1.GetIdx()))
				atomic_distance = np.linalg.norm(coords1-centroid_coords)
			if atomic_distance < radius+radius_buffer:
				res_name.append(res_info(atom1,'res_name'))
				res_number.append(res_info(atom1,'res_number'))
				res_chain.append(res_info(atom1,'chain'))
				res_atom.append(res_info(atom1,'atom_name'))
				if res_info(atom1,'atom_name') in backbone_atoms:
					atom_type.append('Backbone')
				else:
					atom_type.append('Sidechain')
	already_added = []
	if centroid == False:
		for i,atom1 in enumerate(base_mol.GetAtoms()):
			if str(res_info(atom1,'res_name'))+str(res_info(atom1,'res_number')) in include_residues:
				keep_atom = True
			else:
				keep_atom = False
				coords1 = np.asarray(base_mol.GetConformer().GetAtomPosition(atom1.GetIdx()))
			for j,atom2 in enumerate(center_mol.GetAtoms()):
				if keep_atom == True:
					atomic_distance = radius-1
				if keep_atom == False:
					coords2 = np.asarray(center_mol.GetConformer().GetAtomPosition(atom2.GetIdx()))
					atomic_distance = np.linalg.norm(coords1-coords2)
				if atomic_distance < radius:
					if str(res_info(atom1,'atom_name'))+str(res_info(atom1,'res_number')) in already_added:
						continue
					res_name.append(res_info(atom1,'res_name'))
					res_number.append(res_info(atom1,'res_number'))
					res_chain.append(res_info(atom1,'chain'))
					res_atom.append(res_info(atom1,'atom_name'))
					already_added.append(str(res_info(atom1,'atom_name'))+str(res_info(atom1,'res_number')))
					if res_info(atom1,'atom_name') in backbone_atoms:
						atom_type.append('Backbone')
					else:
						atom_type.append('Sidechain')
	new_mol = Chem.RWMol(base_mol)

	res_dict = {'Residue Name':res_name,'Residue Number':res_number,'Residue Atom':res_atom,'Atom Type':atom_type,'Residue Chain':res_chain}

	for atom in reversed(base_mol.GetAtoms()):
		if res_info(atom,'res_number') in res_number:
			if res_info(atom,'chain') == res_chain[res_number.index(res_info(atom,'res_number'))]:
				continue
			else:
				new_mol.RemoveAtom(atom.GetIdx())
		else:
			new_mol.RemoveAtom(atom.GetIdx())
	if centroid is True:
		print('Initial pass results in {} atoms.'.format(new_mol.GetNumAtoms()))

	if centroid is False:
		print('Final pass results in {} atoms.'.format(new_mol.GetNumAtoms()))
		print("Structure saved as active_site_radius_{}.pdb".format(radius))
		Chem.MolToPDBFile(new_mol,'active_site_radius_{}.pdb'.format(radius))

	return new_mol, res_dict


def truncate_new(base_mol, scheme='CA_terminal', skip_residues=['HOH','WAT'], skip_resnumbers=[], remove_resnumbers=[], remove_atoms=[], remove_atom_ids=[], remove_sidechains=[], constrain_atoms=[' CA '],radius=None):

	new_mol = Chem.RWMol(base_mol)
	proline_count = 0
	constrain_list = []
	res_num = []
	res_name = []
	N_terminus = []
	C_terminus = []
	bb_atom_count = 0
	previous_res = None
	remove_ids = []
	fix_bonds = []

##########################################################
# Record information about main chain atom bonds
# Record information about main chain atom bonds
	for atom in reversed(base_mol.GetAtoms()):
		#current_res = (res_info(atom,'chain'),res_info(atom,'res_name'),res_info(atom,'res_number'))
		current_res = define_residue(atom)
		if res_info(atom,'atom_name') in remove_atoms:
			remove_ids.append(atom.GetIdx())
			continue
		if atom.GetIdx() in remove_atom_ids:
			remove_ids.append(atom.GetIdx())
		if current_res[1] in skip_residues:
			continue
		if current_res[2] in skip_resnumbers:
			continue
		if current_res[2] in remove_resnumbers:
			remove_ids.append(atom.GetIdx())
			continue
		if current_res != previous_res:
			if res_info(atom,'atom_name') == ' C  ':
				bb_atom_count += 1
				C_id = atom.GetIdx()
			if res_info(atom,'atom_name') == ' O  ':
				bb_atom_count += 1
				O_id = atom.GetIdx()
			if res_info(atom,'atom_name') == ' N  ':
				bb_atom_count += 1
				N_id = atom.GetIdx()
			if bb_atom_count == 3:
				C_atom = base_mol.GetAtomWithIdx(C_id)
				N_atom = base_mol.GetAtomWithIdx(N_id)
				C_bonds = [res_info(x,'atom_name') for x in C_atom.GetNeighbors()]
				C_bonds_atoms = [x for x in C_atom.GetNeighbors()]
				N_bonds = [res_info(x,'atom_name') for x in N_atom.GetNeighbors()]
				N_bonds_atoms = [x for x in N_atom.GetNeighbors()]
			if bb_atom_count !=3:
				continue
##########################################################
# need to make CA_all capping scheme code
##########################################################
		# Get atom IDs that will be removed according to the
		# CA_terminal capping scheme
		if scheme == 'CA_terminal':
			if res_info(atom,'res_name') == 'PRO':
				N_terminus.append('Keep')
				proline_count += 1
			if res_info(atom,'res_name') != 'PRO':
				if ' C  ' in N_bonds:
					N_terminus.append('Keep')
				else:
					N_terminus.append('Remove')
					cap_atom = new_mol.GetAtomWithIdx(N_id)
					cap_atom.SetAtomicNum(1)
					cap_atom.GetPDBResidueInfo().SetName(' H  ')
					fix_bonds.append(cap_atom.GetIdx())
					for a in range(len(N_bonds)):
						if res_info(N_bonds_atoms[a],'atom_name') == ' H  ':
							remove_ids.append(N_bonds_atoms[a].GetIdx())
			if ' N  ' in C_bonds:
				C_terminus.append('Keep')
			else:
				C_terminus.append('Remove')
				remove_ids.append(O_id)
				cap_atom = new_mol.GetAtomWithIdx(C_id)
				cap_atom.SetAtomicNum(1)
				cap_atom.GetPDBResidueInfo().SetName(' H  ')
				fix_bonds.append(cap_atom.GetIdx())
				for a in range(len(C_bonds)):
					if res_info(C_bonds_atoms[a],'atom_name') == ' H  ':
						remove_ids.append(C_bonds_atoms[a].GetIdx())
			bb_atom_count = 0
			previous_res = current_res

##########################################################
# Need to consolidate this into previous loop for efficiency
##########################################################
	# Remove sidechains specified when the function was called
	for atom in reversed(base_mol.GetAtoms()):
		current_res = (res_info(atom,'chain'),res_info(atom,'res_name'),res_info(atom,'res_number'))
		if current_res[-1] not in remove_sidechains:
			continue
		if res_info(atom,'res_name') in ["PRO"]:
			continue
		if res_info(atom,'atom_name') in [' N  ',' C  ',' O  ',' CA ']:
			continue
		if res_info(atom,'atom_name') != ' H  ':
			for x in atom.GetNeighbors():
				if res_info(x,'atom_name') == ' H  ':
					remove_ids.append(x.GetIdx())
			if res_info(atom,'atom_name') == ' CB ':
				cap_atom = new_mol.GetAtomWithIdx(atom.GetIdx())
				cap_atom.SetAtomicNum(1)
				cap_atom.GetPDBResidueInfo().SetName(' H  ')
				continue
			remove_ids.append(atom.GetIdx())
##########################################################

	# Now actually create the truncated and capped mol object
	for atom in reversed(new_mol.GetAtoms()):
		if res_info(atom,'atom_name') == ' CA ':
			lone_methyl = True
			bound_atoms = []
			for x in atom.GetNeighbors():
				bound_atoms.append(x.GetIdx())
				if res_info(x,'atom_name') != ' H  ':
					lone_methyl = False
			if lone_methyl is True:
				remove_ids.append(atom.GetIdx())
				for a in bound_atoms:
					remove_ids.append(a)
	remove_ids.sort(reverse=True)
	for a in remove_ids:
		new_mol.RemoveAtom(a)
	
	Chem.SanitizeMol(new_mol)
	for atom in reversed(new_mol.GetAtoms()):
		if res_info(atom,'atom_name') == ' H  ':
			if res_info(atom.GetNeighbors()[0],'atom_name') == ' CA ':
				if rdMolTransforms.GetBondLength(new_mol.GetConformer(), atom.GetNeighbors()[0].GetIdx(), atom.GetIdx()) > 1.01:
					rdMolTransforms.SetBondLength(new_mol.GetConformer(),atom.GetNeighbors()[0].GetIdx(),atom.GetIdx(), 0.970)
##########################################################
# Record CA ids to add to constrain list
	for atom in new_mol.GetAtoms():
		if res_info(atom,'atom_name') in constrain_atoms:
			constrain_list.append(atom.GetIdx()+1)

	print('Final active site model contains {} atoms.'.format(new_mol.GetNumAtoms()))
	print("Structure saved as truncated_active_site_radius_{}.pdb".format(radius))
	Chem.MolToPDBFile(new_mol,'truncated_active_site_radius_{}.pdb'.format(radius))

	if proline_count > 0:
		print('WARNING: active site model contains {} proline(s). N atom was kept.'.format(proline_count))
##########################################################
	return new_mol, constrain_list
##########################################################
##########################################################
def obabel_protonate(active_site_pdb_file):

	ob_conv = ob.OBConversion()
	ob_conv.SetInFormat('pdb')
	ob_mol = ob.OBMol()
	ob_conv.ReadFile(ob_mol, str(active_site_pdb_file))

	ob_mol.AddHydrogens(False, True, 7.0) # polaronly, correctForPH, pH
	conv = ob.OBConversion()
	conv.SetInAndOutFormats('pdb', 'pdb')

	his_ids = []

	active_site_mol = Chem.MolFromPDBFile(str(active_site_pdb_file),removeHs=False,sanitize=False)
	for atom in active_site_mol.GetAtoms():
		if res_info(atom,'res_name') == 'HIS':
			his_ids.append(res_info(atom,'res_number'))

	if len(his_ids) > 0:
		print('WARNING: There are {} HIS residues, user should check protonation state and adjust pdb output file accordingly. Default is epsilon_N.'.format(len(np.unique(his_ids))))

	file_name_base = str(active_site_pdb_file).split('.pdb')[0]
	conv.WriteFile(ob_mol, str(file_name_base)+'_protonated'+'.pdb')
	print("Protonated structure saved as "+ str(file_name_base)+'_protonated'+'.pdb')

def show_mol(base_mol):
	mol = Chem.RWMol(base_mol)
	from rdkit.Chem import Draw
	from rdkit.Chem.Draw import rdMolDraw2D
	from rdkit.Chem import rdDepictor
	rdDepictor.SetPreferCoordGen(True)
	from rdkit.Chem.Draw import IPythonConsole
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

if __name__ == "__main__":
	# Do something if this file is invoked on its own
	print(canvas())
