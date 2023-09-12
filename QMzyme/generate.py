###############################################################################
# Code written by Heidi Klem while at
# Colorado State University as a graduate student
# in the Paton and McCullagh groups and at the
# National Institute of Standards and Technology
# as an NRC Postdoc (Fed).
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""Generate QM-based enzyme model."""

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

delimeter="---------------------------------------------------------"

###############################################################################
def res_info(atom,info='atom_name'):
	''' 
	Options for info are 'chain', 'res_name', 'res_number', 'atom_name'.
	'''
	if info == 'res_name':
		return atom.GetPDBResidueInfo().GetResidueName()
	if info == 'atom_name':
		return atom.GetPDBResidueInfo().GetName()
	if info == 'res_number':
		return atom.GetPDBResidueInfo().GetResidueNumber()
	if info == 'chain':
		return atom.GetPDBResidueInfo().GetChainId()
	else:
		print("ERROR: Improper residue info selection. Options are 'chain',"
			  " 'res_name', 'res_number' and 'atom_name'.")

###############################################################################
def define_residue(atom):
	''' 
	Returns a tuple that completely defines the residue that 
	atom belongs to: ('chain','residue name', 'residue number')
	'''
	return (res_info(atom,'chain'),res_info(atom,'res_name'),
			res_info(atom,'res_number'))

###############################################################################
def atom_coords(mol,atom):
	return np.asarray(mol.GetConformer().GetAtomPosition(atom.GetIdx()))

###############################################################################
def download(pdb_list):
	baseUrl = 'http://www.pdb.org/pdb/download/downloadFile.do?fileFormat'+\
	 '=pdb&compression=NO&structureId='
	for structure in pdb_list.split():
		pdbUrl = baseUrl + structure[:4]
		outFileName = structure[:4] + '.pdb'
		print(pdbUrl)

		with urlopen(pdbUrl) as response, open(outFileName, 'wb') as outFile:
			data = response.read()
			outFile.write(data)
			print('Downloading {} as {}.'.format(structure[:4], outFileName))
	print(delimeter)
	return(data)

###############################################################################
def PDB_file_info(pdb_file):
	'''
	Input: PDB file name. Prints helpful information about what is present in 
	the structure.
	'''
	
	base_mol = Chem.MolFromPDBFile(pdb_file,removeHs=False,sanitize=False)
	wat_count, protein_res_count, non_protein_res_count = 0,0,0 
	anion_count, cation_count = 0, 0
	protein_seq, non_protein_seq = [], []
	h_present, no_chain_info = False, False
	previous_res = None

	for atom in base_mol.GetAtoms():
		if ' H' in res_info(atom,'atom_name'):
			h_present = True
		if res_info(atom,'chain') == ' ':
			no_chain_info = True
			atom.GetPDBResidueInfo().SetChainId('X')
		current_res = define_residue(atom)
		if previous_res == current_res:
			continue
		previous_res = current_res
		if current_res[1] in ["WAT","HOH"]:
			wat_count += 1
			continue
		if "+" in current_res[1]:
			cation_count += 1
			cation = current_res[1]
			continue
		if "-" in current_res[1]:
			anion_count += 1
			anion = current_res[1]
			continue
		if atom.GetPDBResidueInfo().GetIsHeteroAtom() is True:
			non_protein_seq.append(current_res)
			non_protein_res_count += 1
			continue
		if current_res[1] in protein_residues:
			protein_seq.append(current_res)
			protein_res_count += 1
		else:
			non_protein_seq.append(current_res)
			non_protein_res_count += 1

	print("Information on PDB file: {}".format(pdb_file))
	if h_present == True:
		print("Hydrogens are present")
	else:
		print("Hydrogens are not present.")
	if no_chain_info == True:
		print("Chain IDs not defined were set to 'X'.")
	print(delimeter)
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
	print(delimeter)
	print("The following {} ligands were detected:"
		  .format(non_protein_res_count))
	if non_protein_res_count > 0:
		for i in range(non_protein_res_count):
			print("{}:".format(i+1))
			print("Chain: {}".format(non_protein_seq[i][0]),
				  " Residue Name: {}".format(non_protein_seq[i][1]),
				  " Residue Number: {}".format(non_protein_seq[i][2]))
			print(".................................................")
	print(delimeter)

###############################################################################
def catalytic_center(pdb_file, res_name=None, res_number=None, chain=None):
	'''
	Currently only supports definition of 1 residue (i.e., lengths 
	of res_name, res_number, and chain inputs must not exceed 1. 
	Future versions will adapt for this functionality. If you wish 
	to use multiple residues right now you can bypass this by 
	manipulating the PDB file: for example, change the residue 
	names to LIG for all residues you want in the catalytic center, 
	and set res_name='LIG'.)
	'''

	count = 0
	definition, catalytic_center=[],[]
	file_suffix = '_'
	if chain is not None:
		definition.append('chain')
		catalytic_center.append(chain)
		file_suffix=file_suffix+'chain{}_'.format(chain)
	if res_name is not None:
		definition.append('res_name')
		catalytic_center.append(res_name)
		file_suffix=file_suffix+'{}_'.format(res_name)
	if res_number is not None:
		definition.append('res_number')
		catalytic_center.append(res_number)
		file_suffix=file_suffix+'{}_'.format(res_number)
	file_suffix = file_suffix+pdb_file.split('.pdb')[0]
	output_file = 'catalytic_center{}.pdb'.format(file_suffix)
	
	protein_mol = Chem.MolFromPDBFile(pdb_file,removeHs=False,sanitize=False)
	catalytic_center_mol = Chem.RWMol(protein_mol)

	for atom in reversed(protein_mol.GetAtoms()):
		current_res = []
		for component in definition:
			current_res.append(res_info(atom,component))
		if current_res != catalytic_center:
			catalytic_center_mol.RemoveAtom(atom.GetIdx())
			#remove_ids.append(atom.GetIdx())
		else:
			count += 1
			continue

	if count==0:
		print("WARNING: No atoms found matching catalytic center definition.")
		print(delimeter)
	if count!=0:
		# sanity check
		raise_issue = False
		previous_res = None
		catalytic_center_residues=[]
		for atom in catalytic_center_mol.GetAtoms():
			current_res = define_residue(atom)
			if previous_res is None:
				previous_res=current_res
				catalytic_center_residues.append(current_res)
			if current_res!=previous_res:
				catalytic_center_residues.append(current_res)
				previous_res=current_res
				raise_issue = True
		if raise_issue is True:
			print("WARNING: Your catalytic center definition is not unique"+
				  " and multiple residues were therefore included: {}."\
				  .format(catalytic_center_residues)+" Please ensure this"+
				  " is what you intended!")
			
		print('Catalytic center contains {} atoms.'
			  .format(catalytic_center_mol.GetNumAtoms()))
		print("Structure saved as {}".format(output_file))
		Chem.MolToPDBFile(catalytic_center_mol,output_file)

		print(delimeter)
		return catalytic_center_mol, protein_mol
###############################################################################
def residue_shell(center_mol,radius,pdb_file=None,base_mol=None,
				  centroid=False,include_residues=[]):
	'''
	Selects all residues that have at least one atom within the 
	cutoff radius from the predefined center_mol object. The radius 
	can be calculated from the center_mol centroid (by setting 
	centroid=True), or from all center_mol atoms (by setting 
	centroid=False). Calculating from all center_mol atoms requires 
	X more loops in the code where X is the number of center_mol 
	atoms. It is recommended to first run a centroid=True pass to 
	initially reduce the total number of atoms, then run a 
	centroid=False pass on the reduced mol object for efficiency. 
	If centroid=False is selected, the code will add a buffer 
	distance equal to the largest distance between the center_mol 
	centroid and a center_mol atom. This will ensure that all atoms 
	within the cutoff radius will be included in this initial pass. 
	You will get the same shell of residues either way, but doing it 
	in this two-stage manner saves time, especially depending on the 
	size of your center_mol.  
	'''

	if base_mol == None:
		if pdb_file == None:
			print("Error: Must define on of the following:" +
				  " pdb_file or base_mol")
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
			if 'H' in res_info(atom1,'atom_name'):
				continue
			if str(res_info(atom1,'res_name'))+\
			   str(res_info(atom1,'res_number')) in include_residues:
				atomic_distance = radius
			else:
				coords1 = np.asarray(base_mol.GetConformer()
									 .GetAtomPosition(atom1.GetIdx()))
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
			if str(res_info(atom1,'res_name'))+\
			   str(res_info(atom1,'res_number')) in include_residues:
				keep_atom = True
			else:
				keep_atom = False
				coords1 = np.asarray(base_mol.GetConformer().
									 GetAtomPosition(atom1.GetIdx()))
			for j,atom2 in enumerate(center_mol.GetAtoms()):
				if keep_atom == True:
					atomic_distance = radius-1
				if keep_atom == False:
					coords2 = np.asarray(center_mol.GetConformer().
										 GetAtomPosition(atom2.GetIdx()))
					atomic_distance = np.linalg.norm(coords1-coords2)
				if atomic_distance < radius:
					if str(res_info(atom1,'atom_name'))+\
					   str(res_info(atom1,'res_number')) in already_added:
						continue
					res_name.append(res_info(atom1,'res_name'))
					res_number.append(res_info(atom1,'res_number'))
					res_chain.append(res_info(atom1,'chain'))
					res_atom.append(res_info(atom1,'atom_name'))
					already_added.append(str(res_info(atom1,'atom_name'))+\
										 str(res_info(atom1,'res_number')))
					if res_info(atom1,'atom_name') in backbone_atoms:
						atom_type.append('Backbone')
					else:
						atom_type.append('Sidechain')
	new_mol = Chem.RWMol(base_mol)

	res_dict = {'Residue Name':res_name,'Residue Number':res_number,
				'Residue Atom':res_atom,'Atom Type':atom_type,
				'Residue Chain':res_chain}

	for atom in reversed(base_mol.GetAtoms()):
		if res_info(atom,'res_number') in res_number:
			if res_info(atom,'chain') == \
			   res_chain[res_number.index(res_info(atom,'res_number'))]:
				continue
			else:
				new_mol.RemoveAtom(atom.GetIdx())
		else:
			new_mol.RemoveAtom(atom.GetIdx())
	if centroid is True:
		print('Initial pass results in {} atoms.'
			  .format(new_mol.GetNumAtoms()))

	if centroid is False:
		print('Final pass results in {} atoms.'.format(new_mol.GetNumAtoms()))
		print("Structure saved as active_site_radius_{}.pdb".format(radius))
		Chem.MolToPDBFile(new_mol,'active_site_radius_{}.pdb'.format(radius))
	
	print(delimeter)
	return new_mol, res_dict

###############################################################################
def truncate_new(base_mol, scheme='CA_terminal', skip_residues=['HOH','WAT'], 
				 skip_resnumbers=[], remove_resnumbers=[], 
				 remove_atom_ids=[], remove_sidechains=[], 
				 keep_backbones=[], constrain_atoms=[' CA '], radius=None):
	
	'''
	This function is called to prepare truncated enzyme models. 
	It will remove atoms based on the defined scheme and any other 
	specifications when the function is called.
	'''
	
	### Heidi's to do: add CA capping scheme, create capping summary as a 
	# returned item, allow capping summary to be a function input so users
	# can specify the cap they want for each residue. 

	new_mol = Chem.RWMol(base_mol)
	proline_count,bb_atom_count = 0,0
	constrain_list,res_num,res_name,N_terminus,C_terminus=[],[],[],[],[]
	previous_res = None
	remove_ids = []

	for atom in reversed(base_mol.GetAtoms()):
		current_res = define_residue(atom)
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
			atom_name = res_info(atom,'atom_name')
			if atom_name == ' C  ':
				bb_atom_count += 1
				C_id = atom.GetIdx()
			if atom_name == ' O  ':
				bb_atom_count += 1
				O_id = atom.GetIdx()
			if atom_name == ' N  ':
				bb_atom_count += 1
				N_id = atom.GetIdx()
			if bb_atom_count == 3:
				C_atom = base_mol.GetAtomWithIdx(C_id)
				N_atom = base_mol.GetAtomWithIdx(N_id)
				C_bonds = [res_info(x,'atom_name') \
						   for x in C_atom.GetNeighbors()]
				C_bonds_atoms = [x for x in C_atom.GetNeighbors()]
				N_bonds = [res_info(x,'atom_name') \
						   for x in N_atom.GetNeighbors()]
				N_bonds_atoms = [x for x in N_atom.GetNeighbors()]
			if bb_atom_count !=3:
				continue

		# CA_all capping scheme
		if scheme == 'CA_all':
			print('The C-alpha only capping scheme is currently'+\
				  '  under development.')		
		
		# CA_terminal capping scheme
		if scheme == 'CA_terminal':
			if current_res[1]  == 'PRO':
				N_terminus.append('Keep')
				proline_count += 1
			if current_res[1] != 'PRO':
				if ' C  ' in N_bonds:
					N_terminus.append('Keep')
				else:
					N_terminus.append('Remove')
					cap_atom = new_mol.GetAtomWithIdx(N_id)
					cap_atom.SetAtomicNum(1)
					cap_atom.GetPDBResidueInfo().SetName(' H* ')
					for a in range(len(N_bonds)):
						if 'H' in res_info(N_bonds_atoms[a],'atom_name'):
							remove_ids.append(N_bonds_atoms[a].GetIdx())
			if ' N  ' in C_bonds:
				C_terminus.append('Keep')
			else:
				C_terminus.append('Remove')
				remove_ids.append(O_id)
				cap_atom = new_mol.GetAtomWithIdx(C_id)
				cap_atom.SetAtomicNum(1)
				cap_atom.GetPDBResidueInfo().SetName(' H* ')
				for a in range(len(C_bonds)):
					if 'H' in res_info(C_bonds_atoms[a],'atom_name'):
						remove_ids.append(C_bonds_atoms[a].GetIdx())
			bb_atom_count = 0
			previous_res = current_res

	if len(remove_sidechains) > 0:
		for atom in reversed(base_mol.GetAtoms()):
			current_res = define_residue(atom)
			if current_res[2] in remove_sidechains:
				atom_name = res_info(atom, 'atom_name')
				if atom_name in (' CA ',' C  ',' O  ',' N  '):
					continue
				if atom_name == ' CB ':
					cap_atom = new_mol.GetAtomWithIdx(atom.GetIdx())
					cap_atom.SetAtomicNum(1)
					cap_atom.GetPDBResidueInfo().SetName(' H* ')
					continue
				if atom_name == ' H  ':
					continue
				if ' CA ' not in [res_info(x,'atom_name') for x in \
				   atom.GetNeighbors()]:
					remove_ids.append(atom.GetIdx())
			else:
				continue
	
	for atom in reversed(new_mol.GetAtoms()):
		if res_info(atom,'atom_name') == ' CA ':
			bound_atoms = []
			lone_methyl = True
			for x in atom.GetNeighbors():
				bound_atoms.append(x.GetIdx())
				if 'H' not in res_info(x,'atom_name'):
					lone_methyl = False
			if lone_methyl is True:
				remove_ids.append(atom.GetIdx())
				for a in bound_atoms:
					remove_ids.append(a)
		
	# Now actually create the truncated and capped mol object
	for a in reversed(np.unique(remove_ids)):
		new_mol.RemoveAtom(int(a))

	Chem.SanitizeMol(new_mol)
	for atom in reversed(new_mol.GetAtoms()):
		if ' H* ' in res_info(atom,'atom_name'):
			if rdMolTransforms.GetBondLength(new_mol.GetConformer(), 
			   atom.GetNeighbors()[0].GetIdx(), atom.GetIdx()) > 1.01:
				# This changes cap H bond length to something more physical
				rdMolTransforms.SetBondLength(new_mol.GetConformer(),
				atom.GetNeighbors()[0].GetIdx(),atom.GetIdx(), 0.970)
		# Record atom ids to add to constrain list
		if res_info(atom,'atom_name') in constrain_atoms:
			constrain_list.append(atom.GetIdx()+1)

	print('Final active site model contains {} atoms.'
		  .format(new_mol.GetNumAtoms()))
	print("Structure saved as truncated_active_site_radius_{}.pdb"
		  .format(radius))
	
	#testing add Hs
	#new_mol = Chem.rdmolops.AddHs(new_mol,addCoords=True)
	Chem.MolToPDBFile(new_mol,'truncated_active_site_radius_{}.pdb'
					  .format(radius))

	if proline_count > 0:
		print("WARNING: Active site model contains {}".format(proline_count) +\
			  " proline(s). The N backbone atom was automatically kept.")
	

	print(delimeter)
	return new_mol, constrain_list

###############################################################################
def obabel_protonate(active_site_pdb_file):

	ob_conv = ob.OBConversion()
	ob_conv.SetInFormat('pdb')
	ob_mol = ob.OBMol()
	ob_conv.ReadFile(ob_mol, str(active_site_pdb_file))

	ob_mol.AddHydrogens(False, True, 7.0) # polaronly, correctForPH, pH
	conv = ob.OBConversion()
	conv.SetInAndOutFormats('pdb', 'pdb')

	his_ids = []

	active_site_mol = Chem.MolFromPDBFile(str(active_site_pdb_file),
										  removeHs=False,sanitize=False)
	for atom in active_site_mol.GetAtoms():
		if res_info(atom,'res_name') == 'HIS':
			his_ids.append(res_info(atom,'res_number'))

	if len(his_ids) > 0:
		print("WARNING: There are {} HIS residues, user should check"+
			  "  protonation state and adjust pdb output file accordingly."+
			  "  Default is epsilon_N.".format(len(np.unique(his_ids))))

	file_name_base = str(active_site_pdb_file).split('.pdb')[0]
	conv.WriteFile(ob_mol, str(file_name_base)+'_protonated'+'.pdb')
	print("Protonated structure saved as "+ str(file_name_base)+\
		  "_protonated"+".pdb")
	print(delimeter)

###############################################################################
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
	print('You know nothing Jon Snow')
