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
import os

protein_residues =	["ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", 
					 "GLH", "GLN", "GLU", "GLY", "HIS", "HID", "HIE", "HIP", 
					 "HYP", "ILE", "LEU", "LYN", "LYS", "MET", "PHE", "PRO",
					 "SER", "THR", "TRP", "TYR", "VAL", "HSE", "HSD", "HSP" ]

elements = ['H','He','Li','Be','B','C','N','O','F','Ne',
           'Na','Mg','Al','Si','P','S','Cl','Ar','K', 'Ca',
           'Sc', 'Ti', 'V','Cr', 'Mn', 'Fe', 'Co', 'Ni',
           'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
           'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
           'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
           'I', 'Xe','Cs', 'Ba','La', 'Ce', 'Pr', 'Nd', 'Pm',
           'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm',
           'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir',
           'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
           'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am',
           'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
           'Rf', 'Db', 'Sg', 'Bh','Hs', 'Mt', 'Ds', 'Rg', 'Cn',
           'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

delimeter="---------------------------------------------------------"

###############################################################################
def res_info(atom,info='atom_name'):
	''' 
	The atom input is an atom from an rdkit mol object. 
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
	The atom input is an atom from an rdkit mol object. 
	'''
	return (res_info(atom,'chain'),res_info(atom,'res_name'),
			res_info(atom,'res_number'))

###############################################################################
def atom_coords(mol,atom):
	'''
	The mol input is an rdkit mol object. 
	The atom input is an atom from an rdkit mol object. 
	'''
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
	#return(data)

###############################################################################
def get_PDB_info(pdb_file):
	print("The function 'get_PDB_info()' is deprecated. Please use "
		  "'check_pdb()' instead.")

def parse_pdb(line,data=''):
    return_data = ''
    if data=='record_type':
        return_data=str(line[0:6]).split()[0]
    if data in ['atom_num','atom_number','atom_serial_number','atom_id']:
        return_data=str(line[6:11])
    if data in ['atom','atom_name']:
        return_data=str(line[12:16])
    if data in ['altloc','alt_loc','alternate_location_indicator']:
        return_data=str(line[16])
    if data in ['residue_name','res','residue','res_name','resname']:
        return_data=str(line[17:20])
    if data in ['chain_identifier','chain_id','chain']:
        return_data=str(line[21])
    if data in ['resid','res_id','residue_id','res_num',
				'residue_number','resnum','residue_sequence_number']:
        return_data=str(line[22:26])
    if data=='occupancy':
        return_data=str(line[54:60])
    if data=='temperature_factor':
        return_data=str(line[60:66])
    if data in ['segid','seg_id','segment_identifier']:
        return_data=str(line[72:76])
    if data in ['element_symbol','element','element_name']:
        return_data=str(line[76:78])
    if data=='charge':
        return_data=str(line[78:80])
    if return_data in ['','\n']:
        raise ValueError('When parsing the pdb for {} nothing was found.'\
						 .format(str(data)))
        pass
    return return_data


###############################################################################
def check_pdb(pdb_file, clean=True):
	'''
	Input: PDB file name. Prints helpful information about what is present in 
	the structure. Will identify if there are any atom type formatting issues,
	that would cause the code to break in preceding steps. If issues are found
	a new pdb file with the suffix "_cleaned.pdb" will be created.
	'''

	wat_count, protein_res_count, non_protein_res_count = 0,0,0 
	anion_count, cation_count = 0, 0
	data, protein_seq, non_protein_seq = [], [], []
	h_present, no_chain_info = False, False
	previous_res = None
	edits_made = False
	residue_count = 0
	res_sequence = []
	res_num_previous = None
	

	#TO DO: run initial check over PDB to identify what sections are there and what aren't
	# for example, some PDB files might not have chain info.  
	# once that is done, we can make sure all the included sections are formatted correctly.
	# for example, that the section is right justified if it's supposed to be. 
	# Ran into this issue with 4 digit residue IDs (from solvent molecules) that 
	# were left justified instead of right justified. So RDKit thought residue
	# ' 133' was the same as residue ' 1330'. The residue id of the 4 digit res
	# needs to be shifted to the left (i.e., the space needs to be removed. 
	# The correct spacing is as follows '   1', '  10', ' 100', '1000' according
	# to the columnspace of the residue_id in PDB format (columns 23-26, 1-indexing)) 

	pdb_format_sections=['record_type', 'atom_number','atom_name','alt_loc',
	    				  'residue_name','chain_id','residue_number',
						  'occupancy','temperature_factor','segment_id',
						  'element_symbol','charge']
	residues_reordered=False
	if clean==True:
		with open(pdb_file,'r') as f:
			data = f.readlines()
		for i,line in enumerate(data):
			if parse_pdb(line,data='record_type') in ['ATOM','HETATOM']:
				atom_type = parse_pdb(line,data='atom_name')
				res_num = parse_pdb(line,data='res_num')
				res_name = parse_pdb(line,data='res_name')
				if res_num != res_num_previous:
					if res_name not in ['WAT','HOH','TIP']:
						residue_count += 1
						if int(res_num.split()[0]) != residue_count:
							residues_reordered = True
							edits_made=True
							if residue_count>999:
								res_num = str(residue_count)
							elif residue_count>99:
								res_num = ' '+str(residue_count)
							elif residue_count>9:
								res_num = '  '+str(residue_count)
							else:
								res_num = '   '+str(residue_count)
							data[i]=line[0:22]+res_num+line[26:]
							line=line[0:22]+res_num+line[26:]
							if parse_pdb(data[i+1],data='res_name') == res_name:
								residue_count -= 1
				res_num_previous = res_num
				try:
					element = parse_pdb(line,data='element')
				except ValueError:
					element=None
				if atom_type[0]!=' ':
					if atom_type[0]=='H':
						if len(atom_type)==4:
							if element is None:
								edits_made=True
								data[i]=line[0:76]+' H'+line[78:]
							continue
				if atom_type[:2] in elements:
					print('yes')
					if element is None:
						edits_made=True
						data[i]=line[0:76]+atom_type[:2]+line[78:]
					elif atom_type[:2] == element:
						continue
				if atom_type[1] in elements:
					if element is None:
						edits_made=True
						data[i]=line[0:76]+atom_type[:2]+line[78:]
					elif atom_type[:2] == element:
						continue
					#else:
				else:
					if atom_type[1]+atom_type[2].lower() in elements:
						edits_made=True
						data[i]=line[0:12]+atom_type[1:]+' '+line[16:76]+\
								  atom_type[1:3]+line[78:]
						print('SUSPECTED ISSUE FOUND ON ATOM{}: Element symbol'\
							  .format(parse_pdb(line,data='atom_number')) + 
							  ' with two letters ({}) must begin in the 13th'\
							  .format(atom_type[1:3]) + ' columnspace according' +
							  ' to PDB format.')                    
					elif atom_type[2] in elements:
						edits_made=True
						unk_sym = atom_type[:2].split()[0]
						if len(unk_sym) == 2:
							data[i] = line[0:12]+' '+atom_type[2:]+' '+\
									  line[16:76]+' '+atom_type[2]+line[78:]
						if len(unk_sym) == 1:
							data[i] = line[0:12]+' '+atom_type[2:]+' '+\
									  line[16:76]+' '+atom_type[2]+line[78:]
						print('SUSPECTED ISSUE FOUND ON ATOM{}'\
							  .format(parse_pdb(line,data='atom_number')) +
							  ': Atom type ({}) '.format(atom_type) + 
							  'is preceeded by unknown symbols ({}). These'\
							  .format(atom_type[:2]) +' will be removed.')
												   
		if residues_reordered==True:
			print('WARNING: Some residues were renumbered.')
		if edits_made is True:
			new_file = pdb_file.split('.pdb')[0]+'_original'+'.pdb'
			cmd = "cp {} {}".format(pdb_file,new_file)
			os.system(cmd)
			print('Saved original pdb file as {}.'.format(new_file)+
				  ' Suspected issue(s) have been fixed in {}.'.format(pdb_file))
			with open(pdb_file, 'w+') as g:
				for i,line in enumerate(data):
					g.writelines(data[i])
		print(delimeter)	
	mol = Chem.MolFromPDBFile(pdb_file,removeHs=False,sanitize=False)
	for atom in mol.GetAtoms():
		if ' H' in res_info(atom,'atom_name'):
			h_present = True
		if res_info(atom,'chain') == ' ':
			no_chain_info = True
			atom.GetPDBResidueInfo().SetChainId('X')
		current_res = define_residue(atom)
		if previous_res == current_res:
			continue
		previous_res = current_res
		if current_res[1] in ["WAT","HOH","TIP"]:
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
	print("Total number of atoms: {}".format(mol.GetNumAtoms()))
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
	print("The following {} non-protein residues were detected:"
		  .format(non_protein_res_count))
	if non_protein_res_count > 0:
		for i in range(non_protein_res_count):
			#print("{}:".format(i+1))
			print("Chain: {}".format(non_protein_seq[i][0]),
				  " Residue Name: {}".format(non_protein_seq[i][1]),
				  " Residue Number: {}".format(non_protein_seq[i][2]))
			#print(".................................................")
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
def residue_shell(catalytic_center_mol=None, catalytic_center_pdb=None,
				  distance_cutoff=0, extended_pdb=None, extended_mol=None,
				  centroid=False,include_residues=[]):
	'''
	Selects all residues that have at least one atom within the 
	cutoff distance from the predefined center_mol object. The distance 
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
	within the cutoff distance will be included in this initial pass. 
	You will get the same shell of residues either way, but doing it 
	in this two-stage manner saves time, especially depending on the 
	size of your catalytic_center_mol.  
	'''
	if distance_cutoff==0:
		raise ValueError('Please specify a distance cutoff by distance_cutoff=int')	

	if catalytic_center_mol == None:
		if catalytic_center_pdb == None:
			raise ValueError(' Must define on of the following:' +
							 'catalytic_center_pdb or catalytic_center_mol.')
		else:
			catalytic_center_mol = Chem.MolFromPDBFile(extended_pdb,
													   removeHs=False,
													   sanitize=False)

	if extended_mol == None:
		if extended_pdb == None:
			raise ValueError("Must define on of the following:" +
							 " pdb_file or extended_mol")
		extended_mol = Chem.MolFromPDBFile(extended_pdb,
										   removeHs=False,
										   sanitize=False)

	res_name, res_number, res_chain, res_atom, atom_type = [], [], [], [], []
	N_termini_interactions = []
	side_chain_interactions = []
	C_termini_interactions = []
	current_res = None
	backbone_atoms=[' O  ',' C	',' N  ',' CA ']

	if centroid == True:
		centroid_coords = np.asarray(Chem.rdMolTransforms.\
			ComputeCentroid((catalytic_center_mol.GetConformer())))
		distances = [np.linalg.norm(np.asarray(catalytic_center_mol.GetConformer().\
			GetAtomPosition(atom.GetIdx()))-centroid_coords) for atom in \
			catalytic_center_mol.GetAtoms()]
		distance_buffer = np.max(distances)

		for atom1 in extended_mol.GetAtoms():
			if ' H' in res_info(atom1,'atom_name'):
				continue
			current_res=define_residue(atom1)
			#if str(res_info(atom1,'res_name'))+\
			#   str(res_info(atom1,'res_number')) in include_residues:
			if current_res in include_residues:
				atomic_distance = distance_cutoff
			else:
				coords1 = np.asarray(extended_mol.GetConformer()
									 .GetAtomPosition(atom1.GetIdx()))
				atomic_distance = np.linalg.norm(coords1-centroid_coords)
			if atomic_distance < distance_cutoff+distance_buffer:
				if 'TIP' in current_res[1]:
					print(current_res[2])
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
		for i,atom1 in enumerate(extended_mol.GetAtoms()):
			current_res=define_residue(atom1)
			#if str(res_info(atom1,'res_name'))+\
			#   str(res_info(atom1,'res_number')) in include_residues:
			if current_res in include_residues:
				keep_atom = True
			else:
				keep_atom = False
				coords1 = np.asarray(extended_mol.GetConformer().
									 GetAtomPosition(atom1.GetIdx()))
			for j,atom2 in enumerate(catalytic_center_mol.GetAtoms()):
				if keep_atom == True:
					atomic_distance = distance_cutoff-1
				if keep_atom == False:
					coords2 = np.asarray(catalytic_center_mol.GetConformer().
										 GetAtomPosition(atom2.GetIdx()))
					atomic_distance = np.linalg.norm(coords1-coords2)
				if atomic_distance < distance_cutoff:
					#if str(res_info(atom1,'atom_name'))+\
					#   str(res_info(atom1,'res_number')) in already_added:
					if current_res in already_added:
						continue
					res_name.append(current_res[1])
					res_number.append(current_res[2])
					res_chain.append(current_res[0])
					res_atom.append(res_info(atom1,'atom_name'))
					#already_added.append(str(res_info(atom1,'atom_name'))+\
					#					 str(res_info(atom1,'res_number')))
					already_added.append(current_res)
					if res_info(atom1,'atom_name') in backbone_atoms:
						atom_type.append('Backbone')
					else:
						atom_type.append('Sidechain')
	new_mol = Chem.RWMol(extended_mol)
	
	res_dict = {'Residue Name':res_name,'Residue Number':res_number,
				'Residue Atom':res_atom,'Atom Type':atom_type,
				'Residue Chain':res_chain}
	for atom in reversed(extended_mol.GetAtoms()):
		current_res = define_residue(atom)
		if current_res[2] in res_number:
			if current_res[1] == res_name[res_number.index(current_res[2])]:
				continue
			#if res_info(atom,'chain') == \
			#   res_chain[res_number.index(res_info(atom,'res_number'))]:
			#	continue
			else:
				new_mol.RemoveAtom(atom.GetIdx())
		else:
			new_mol.RemoveAtom(atom.GetIdx())
	if centroid is True:
		print('Initial pass results in {} atoms.'
			  .format(new_mol.GetNumAtoms()))

	if centroid is False:
		print('Final pass results in {} atoms.'.format(new_mol.GetNumAtoms()))
		print("Structure saved as active_site_cutoff_{}.pdb".format(distance_cutoff))
		Chem.MolToPDBFile(new_mol,'active_site_cutoff_{}.pdb'.format(distance_cutoff))
	
	print(delimeter)
	return new_mol, res_dict

###############################################################################
def truncate_new():
	raise Exception("The function 'truncate_new()' is deprecated. Please use "
				 	"'truncate()' instead")

###############################################################################
def add_H(pdb_file=None, output_file=None, remove_files=False):
	'''
	Uses reduce function from AmberTools to add hydrogens.
	'''
	if output_file is None:
		output_file=pdb_file.split('.pdb')[0]+'_reduced.pdb'
	cmd = "reduce -NOFLIP -Quiet {} > {}".format(pdb_file,output_file)
	os.system(cmd)
	new_mol = Chem.MolFromPDBFile(output_file,removeHs=False,sanitize=False)
	if remove_files is True:
		cmd = "rm {} {}".format(pdb_file, output_file)
		os.system(cmd)	

	return new_mol

###############################################################################
def truncate(mol=None, pdb_file=None, scheme='CA_terminal', 
			 skip_residues=['HOH','WAT'], skip_resnumbers=[], 
			 remove_resnumbers=[], remove_atom_ids=[], remove_sidechains=[], 
			 keep_backbones=[], constrain_atoms=[' CA '], add_hydrogens=False):
	'''
	This function is called to prepare truncated enzyme models. 
	It will remove atoms based on the defined scheme and any other 
	specifications when the function is called.
	'''
	
	### Heidi's to do: add CA capping scheme, create capping summary as a 
	# returned item, allow capping summary to be a function input so users
	# can specify the cap they want for each residue. 
	
	if mol is None:
		if pdb_file is None:
			raise ValueError('You must specify at least one of the following:'
							 ' mol={rdkit mol object), or pdb_file={str}.')
		else:
			mol = Chem.MolFromPDBFile(pdb_file,removeHs=False,sanitize=False)	

	if add_hydrogens is True:
		Chem.MolToPDBFile(mol,'temp.pdb')
		mol = add_H('temp.pdb',remove_files=True)

	new_mol = Chem.RWMol(mol)
	proline_count,bb_atom_count = 0,0
	constrain_list,res_num,res_name,N_terminus,C_terminus=[],[],[],[],[]
	previous_res = None
	remove_ids = []

	for atom in reversed(mol.GetAtoms()):
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
				C_atom = mol.GetAtomWithIdx(C_id)
				N_atom = mol.GetAtomWithIdx(N_id)
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
		for atom in reversed(mol.GetAtoms()):
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
	print("Structure saved as truncated_active_site.pdb")
	
	#testing add Hs
	#new_mol = Chem.rdmolops.AddHs(new_mol,addCoords=True)
	Chem.MolToPDBFile(new_mol,'truncated_active_site.pdb')

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
def show_mol(mol):
	mol = Chem.RWMol(mol)
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
