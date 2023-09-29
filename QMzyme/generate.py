###############################################################################
# Code written by Heidi Klem while at
# Colorado State University as a graduate student
# in the Paton and McCullagh groups and at the
# National Institute of Standards and Technology
# as an NRC Postdoc (Fed).
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

'''Generate QM-based enzyme model.'''

import numpy as np
from rdkit import Chem
from urllib.request import urlopen
from rdkit.Chem import rdMolTransforms
import os
from rdkit.Chem import rdDistGeom


delimeter='---------------------------------------------------------'
protein_residues =  ['ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'CYM', 'CYS', 'CYX',
					 'GLH', 'GLN', 'GLU', 'GLY', 'HIS', 'HID', 'HIE', 'HIP',
					 'HYP', 'ILE', 'LEU', 'LYN', 'LYS', 'MET', 'PHE', 'PRO',
					 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HSE', 'HSD', 'HSP' ]

solvent_list=['HOH','WAT','T3P','SOL'] 

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



class generate_model:

	def __init__(self, calculation='QM-only', protein_file=None, pdb_code=None):
		'''
		Initialize QMzyme model. 
		calculation		- string defining the type of calculation this 
					model will be used for. Only option currently is 
					'QM-only'. 'QM/MM' will be supported in the near future. 
		protein_file		- string defining the protein or system .pdb 
					file. This can be preprocessed, or a fresh download 
					from rcsb.org. 
		pdb_code		- string defining the PDB ID code. This option is 
					required if you do not currently have a protein pdb 
					file for your system of interest. Defining the
					4-letter PDB code here will automatically download 
					it from the Protein Data Bank. It is highly recommended 
					that you double check the downloaded file and become 
					acquainted with it. 
		'''
		
		if protein_file is None:
			if pdb_code is None:
				raise ValueError("Either protein_file, or pdb_code must"+
								 " be defined when instantiating model."+
								 " If you currently do not have a pdb file"+
								 " use the pdb_code argument to automatically"+
								 " download the file from the PDB server.")
			else:
				protein_file = self.download(pdb_code)	

		self.protein_file=protein_file
		self.protein_mol = Chem.MolFromPDBFile(protein_file,
											   removeHs=False,
											   sanitize=False)
		self.protein_prefix = protein_file.split('.pdb')[0]
		self.calculation = calculation

###############################################################################
	def download(self, pdb_code):
		base_url = 'http://www.pdb.org/pdb/download/downloadFile.do?fileFormat'+\
		 '=pdb&compression=NO&structureId='
		for structure in pdb_code.split():
			pdb_url = base_url + structure[:4]
			output_file = structure[:4] + '.pdb'
			print(pdb_url)

			with urlopen(pdb_url) as response, open(output_file, 'wb') as outfile:
				data = response.read()
				outfile.write(data)
				print("Downloading {} as {}.".format(structure[:4], output_file))
		print(delimeter)
		return output_file

###############################################################################
	def res_info(self,atom,info='atom_name'):
		'''
		Function to expedite parsing PDB information with rdkit.
		atom			- an atom from an rdkit mol object.
		info			- string defining what information you want. Options 
					are 'chain', 'res_name', 'res_number', 'atom_name'.

		Returns requested atom information.
		'''

		if info == 'res_name':
			return atom.GetPDBResidueInfo().GetResidueName().split()[0]
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
	def define_residue(self,atom):
		'''
		Function to completely define the residue an atom belongs to with rdkit.
		atom			- an atom from an rdkit mol object.

		Returns tuple of format: ('chain','residue name', 'residue number')
		'''

		return (self.res_info(atom,'chain'),self.res_info(atom,'res_name').split()[0],
				self.res_info(atom,'res_number'))


###############################################################################
	def parse_pdb(self,line,info=''):
		'''
		Function to grab specific information from a line in a PDB file.
		line			- string of a line from a PDB file.
		info			- string defining what information you wish to grab 
					from that line. Options correspond to conventional PDB 
					column sections: 'record_type', 'atom_number', 'atom_name',
					'alt_loc', 'residue_name', 'chain_id', 'residue_number', 
					'occupancy', 'temperature_factor', 'segid', 'element_symbol', 
					'charge'.  

		Returns string of requested information.
		'''

		return_info = ''
		if info=='record_type':
			return_info=str(line[0:6]).split()[0]
		if info in ['atom_num','atom_number','atom_serial_number','atom_id']:
			return_info=str(line[6:11])
		if info in ['atom','atom_name']:
			return_info=str(line[12:16])
		if info in ['altloc','alt_loc','alternate_location_indicator']:
			return_info=str(line[16])
		if info in ['residue_name','res','residue','res_name','resname']:
			return_info=str(line[17:20])
		if info in ['chain_identifier','chain_id','chain']:
			return_info=str(line[21])
		if info in ['resid','res_id','residue_id','res_num',
					'residue_number','resnum','residue_sequence_number']:
			return_info=str(line[22:26])
		if info=='occupancy':
			return_info=str(line[54:60])
		if info=='temperature_factor':
			return_info=str(line[60:66])
		if info in ['segid','seg_id','segment_identifier']:
			return_info=str(line[72:76])
		if info in ['element_symbol','element','element_name']:
			return_info=str(line[76:78])
		if info=='charge':
			return_info=str(line[78:80])
		if return_info in ['','\n']:
			raise ValueError("When parsing the pdb for {} nothing was found."\
							 .format(str(info)))
			pass
		return return_info
###############################################################################
	def get_atoms(self):
		return self.protein_mol.GetAtoms()

###############################################################################
	def check_pdb(self,clean=True):
		'''
		Function to assess PDB format, fix any issues that might break 
		the QMzyme code, and gather useful basic information.
		clean			- boolean, default=True. If False, the PDB format will 
					not be checked. Usually not a problem when downloading 
					straight from rcsb.org.

		Prints out basic information of the pdb file, and  if 
		issues are found the PDB file will be fixed, and the original 
		will be copied to a new file with the suffix '_original.pdb'.
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

		with open(self.protein_file,'r') as f:
			data = f.readlines()
		for line in data:
			if 'HET   ' in line:
				non_protein_res_count += 1
				chain, name, number = line[12], line[7:10], line[13:17]
				non_protein_residues['Chain'].append(chain)
				non_protein_residues['Name'].append(name)
				non_protein_residues['Number'].append(number)
			if 'HETNAM' in line:
				non_protein_chemical_name[line[11:14]]=line[15:].split('  ')[0]

		residues_reordered=False
		if clean==True:
			with open(self.protein_file,'r') as f:
				data = f.readlines()
			for i,line in enumerate(data):
				if self.parse_pdb(line,info='record_type') in ['ATOM','HETATOM']:
					atom_type = self.parse_pdb(line,info='atom_name')
					res_num = self.parse_pdb(line,info='res_num')
					res_name = self.parse_pdb(line,info='res_name')
					if res_num != res_num_previous:
						if res_name not in solvent_list:
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
								if self.parse_pdb(data[i+1],info='res_name') == res_name:
									residue_count -= 1
					res_num_previous = res_num
					try:
						element = self.parse_pdb(line,info='element')
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
							print("SUSPECTED ISSUE FOUND ON ATOM{}: Element symbol"\
								  .format(self.parse_pdb(line,info='atom_number'))+
								  " with two letters ({}) must begin in the 13th"\
								  .format(atom_type[1:3])+" columnspace according"+
								  " to PDB format.")

						elif atom_type[2] in elements:
							edits_made=True
							unk_sym = atom_type[:2].split()[0]
							if len(unk_sym) == 2:
								data[i] = line[0:12]+' '+atom_type[2:]+' '+\
										  line[16:76]+' '+atom_type[2]+line[78:]
							if len(unk_sym) == 1:
								data[i] = line[0:12]+' '+atom_type[2:]+' '+\
										  line[16:76]+' '+atom_type[2]+line[78:]
							print("SUSPECTED ISSUE FOUND ON ATOM{}"\
								  .format(self.parse_pdb(line,info='atom_number')) +
								  ": Atom type ({}) ".format(atom_type) +
								  "is preceeded by unknown symbols ({}). These"\
								  .format(atom_type[:2])+" will be removed.")

			if residues_reordered==True:
				print("WARNING: Some residues were renumbered.")
			if edits_made is True:
				new_file = self.protein_prefix+'_original.pdb'
				cmd = 'cp {} {}'.format(self.protein_file,new_file)
				os.system(cmd)
				print("Saved original pdb file as {}.".format(new_file)+
					  " Suspected issue(s)/warning(s) have been resolved in {}."\
					  .format(self.protein_file))
				with open(self.protein_file, 'w+') as g:
					for i,line in enumerate(data):
						g.writelines(data[i])
				print(delimeter)

		for atom in self.protein_mol.GetAtoms():
			if ' H' in self.res_info(atom,'atom_name'):
				h_present = True
			if self.res_info(atom,'chain') == ' ':
				no_chain_info = True
				atom.GetPDBResidueInfo().SetChainId('X')
			current_res = self.define_residue(atom)
			if previous_res == current_res:
				continue
			previous_res = current_res
			if current_res[1] in solvent_list:
				wat_count += 1
				continue
			if atom.GetPDBResidueInfo().GetIsHeteroAtom() is True:
				non_protein_seq.append(current_res)
				continue
			if current_res[1] in protein_residues:
				protein_seq.append(current_res)
				protein_res_count += 1
			else:
				non_protein_seq.append(current_res)
	
		print("Information on PDB file: {}".format(self.protein_file))
		if h_present == True:
			print("Hydrogens are present")
		else:
			print("Hydrogens are not present.")
		if no_chain_info == True:
			print("Chain IDs not defined were set to 'X'.")
		print(delimeter)
		print("Total number of atoms: {}".format(self.protein_mol.GetNumAtoms()))
		print("Water molecules: {}".format(wat_count))
		print("Standard amino acid residues: {}".format(protein_res_count))
		print(delimeter)
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
			
		print(delimeter)
		self.protein_res_count = protein_res_count
		self.h_present = h_present
		self.non_protein_residues = non_protein_residues
		self.non_protein_residue_count = non_protein_res_count
		self.non_protein_chemical_names = non_protein_chemical_name
		
###############################################################################
	def atom_coords(self,mol,atom):
		'''
		mol			- rdkit mol object.
		atom			- atom from rdkit mol object.

		Returns numpy array of the atomic cartesian coorinates.
		'''

		return np.asarray(mol.GetConformer().GetAtomPosition(atom.GetIdx()))

###############################################################################
	def catalytic_center(self, res_name=None, res_number=None, chain=None, 
						 output_file=None, save_file=True):
		'''
		Function to define the center of the QMzyme model. This is 
		typically the ligand/substrate. Currently only supports the 
		definition of 1 residue (i.e., lengths of res_name, 	
		res_number, and chain inputs must not exceed 1. Future 
		versions will adapt for this functionality. If you wish to 
		use multiple residues right now you can bypass this by
		manipulating the PDB file: for example, change the residue
		names to LIG for all residues you want in the catalytic center,
		and set res_name='LIG'.)
		res_name		- string of the three-letter residue name.
		res_number		- integer of the residue number.
		chain			- string of the chain the residue of interest is in.
		output_file		- string to define the name of the output .pdb file.
					Generates the following attributes of the class object: a 
					catalytic_center_mol rdkit mol object, catalytic_center_definition, 
					and creates a .pdb file containing onlt the catalytic center atoms. 
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

		if output_file is None:
			output_file = '{}_catalytic_center{}.pdb'\
						  .format(self.protein_prefix,file_suffix[:-1])
		self.cat_center = 'catalytic_center{}'.format(file_suffix[:-1])
		catalytic_center_mol = Chem.RWMol(self.protein_mol)

		for atom in reversed(self.protein_mol.GetAtoms()):
			current_res = []
			for component in definition:
				current_res.append(self.res_info(atom,component))
			if current_res != catalytic_center:
				catalytic_center_mol.RemoveAtom(atom.GetIdx())
				#remove_ids.append(atom.GetIdx())
			else:
				count += 1
				continue

		if count==0:
			print("WARNING: No atoms found matching"+
                  " catalytic center definition.")
			print(delimeter)
		if count!=0:
			# sanity check
			raise_issue = False
			previous_res = None
			catalytic_center_residues=[]
			for atom in catalytic_center_mol.GetAtoms():
				current_res = self.define_residue(atom)
				if previous_res is None:
					previous_res=current_res
					catalytic_center_residues.append(current_res)
				if current_res!=previous_res:
					catalytic_center_residues.append(current_res)
					previous_res=current_res
					raise_issue = True
			if raise_issue is True:
				print("WARNING: Catalytic center definition is not unique"+
					  " and multiple residues were therefore included: {}."\
					  .format(catalytic_center_residues)+" Please ensure this"+
					  " is what you intended!")

			print("Catalytic center contains {} atoms."
				  .format(catalytic_center_mol.GetNumAtoms()))
			if save_file is True:
				print("Structure saved as {}".format(output_file))
				Chem.MolToPDBFile(catalytic_center_mol,output_file)
			if save_file is False:
				self.catalytic_center_pdb = Chem.MolToPDBFile(catalytic_center_mol,output_file)

			self.catalytic_center_definition = current_res
			self.catalytic_center_mol = catalytic_center_mol
			print(delimeter)
			#return catalytic_center_mol, protein_mol

###############################################################################
	def active_site(self, distance_cutoff=0, include_residues=[],
					output_file=None, save_file=True,
					solvent=solvent_list, intermediate_mol=None, center=None):
		'''
		Function that selects all residues that have at least 
		one atom within the cutoff distance from the predefined catalytic 
		center atoms. If there are residues you want to include that may 
		not be within the distance cutoff you can specify that in 
		include_residues=[].
		distance_cutoff		- integer
		include_residues	- list of a tuple of strings:
					(str(chain),str(res_name),int(res_number)). 
					I.e., include_residues=[('A','GLY',101),('A','ASP',20)].
					It follows the define_residue() format.
		output_file		- string defining the .pdb output file name.
		intermediate_mol	- rdkit mol object to be used instead of the 
			protein mol that was created during model initialization. This is
			useful if, for example, you are creating multiple active sites 
			from the same initial PDB with varying cutoff distances, in 
			which case the intermediate_mol object must have been generated
			with a larger cutoff distance than what is currently specified.
		center		- (Work in progress) string defining the .pdb file 
			of your catalytic center. Default is None. This is helpful to 
			use if you want to alter the ligand in any way, for example, 
			add hydrogens to certain atoms. Please not the coordinates of 
			the molecule need to be consistent with the coordinates of 
			the protein.

		Generates active_site_mol attribute, and creates .pdb file. 
		'''

		if distance_cutoff==0:
			raise ValueError("Please specify a distance cutoff by"
                             " 'distance_cutoff={int}'.")

		res_name, res_number, res_chain = [], [], []
		add_residue=[]
		if intermediate_mol is not None:
			mol = intermediate_mol
		else:
			mol = self.protein_mol

		temp_mol = Chem.RWMol(mol)
		new_mol = Chem.RWMol(mol)
		if center is not None:
			self.catalytic_center_mol = Chem.MolFromPDBFile(center)

		centroid_coords = np.asarray(Chem.rdMolTransforms.\
			ComputeCentroid((self.catalytic_center_mol.GetConformer())))
		distances = [np.linalg.norm(np.asarray(self.catalytic_center_mol.\
            GetConformer().GetAtomPosition(atom.GetIdx()))-centroid_coords)\
            for atom in self.catalytic_center_mol.GetAtoms()]
		distance_buffer = np.max(distances)

		for atom in reversed(mol.GetAtoms()):
			if ' H' in self.res_info(atom,'atom_name'):
				continue
			current_res=self.define_residue(atom)
			if current_res in include_residues:
				continue
			if current_res in add_residue:
				continue
			else:
				coords = self.atom_coords(mol,atom)
				atomic_distance = np.linalg.norm(coords-centroid_coords)
				if atomic_distance < distance_cutoff+distance_buffer:
					add_residue.append(current_res)
				else:
					temp_mol.RemoveAtom(atom.GetIdx())

		# second pass: goes over all atoms in catalytic center
		keep_residue=[] 
		for atom1 in reversed(temp_mol.GetAtoms()):
			current_res=self.define_residue(atom1)
			if current_res in include_residues:
				keep_residue.append(current_res)
				continue
			if current_res in keep_residue:
				continue
			else:
				coords1 = self.atom_coords(temp_mol,atom1)
			for atom2 in self.catalytic_center_mol.GetAtoms():
				if current_res in keep_residue:
					continue
				coords2 = self.atom_coords(self.catalytic_center_mol,atom2)
				#coords2 = np.asarray(self.catalytic_center_mol.GetConformer().
				#                    GetAtomPosition(atom2.GetIdx()))
				atomic_distance = np.linalg.norm(coords1-coords2)
				if atomic_distance < distance_cutoff:
					res_chain.append(current_res[0])
					res_name.append(current_res[1])
					res_number.append(current_res[2])
					keep_residue.append(current_res)

		for atom in reversed(mol.GetAtoms()):
			current_res=self.define_residue(atom)
			if current_res not in keep_residue:
				new_mol.RemoveAtom(atom.GetIdx())

		res_dict = {'Residue Chain':res_chain,
					'Residue Name':res_name,
					'Residue Number':res_number}

		print("Active site contains {} atoms.".format(new_mol.GetNumAtoms()))
		if save_file is True:
			if output_file is None:
				output_file = '{}_{}_active_site_distance_cutoff_{}.pdb'\
						  	  .format(self.protein_prefix,self.cat_center,distance_cutoff)
			print("Structure saved as {}".format(output_file))
			Chem.MolToPDBFile(new_mol,output_file)
		self.distance_cutoff=distance_cutoff	
		self.active_site_mol = new_mol
		self.active_site_residues = res_dict
		print(delimeter)
		#return new_mol, res_dict

###############################################################################
	def add_H(self, pdb_file=None, output_file=None):
		'''
		Function that calls the reduce function from AmberTools 
		to add hydrogens.
		pdb_file		- string defining the .pdb file to be reduced.
		output_file		- string defining the reduced .pdb file name.

		Returns new rdkit mol object that has been reduced and saves 
		the reduced .pdb. 
		'''
		
		if pdb_file is None:
			raise ValueError("PDB file must be defined.")
		
		if output_file is None:
			output_file=pdb_file.split('.pdb')[0]+'_reduced.pdb'
		cmd = "reduce -NOFLIP -Quiet {} > {}".format(pdb_file,output_file)
		os.system(cmd)
		new_mol = Chem.MolFromPDBFile(output_file,
                                      removeHs=False,
                                      sanitize=False)

		return new_mol

###############################################################################
	def protonate_solvent(self,mol=None,solvent=solvent_list):
		'''
		Function to add hydrogens to water molecules because reduce does not
		recognize all water residue names. I found that for some weird 
		reason you cannot do rdkit's Chem.rdmolops.AddHs and include 
		PDB residue information when Hs are already present, so we can't
		simply add Hs to the relevant water oxygen atoms because then
		the PDB information cannot be included if any other Hs are present
		in the structure. To avoid that situation, this function separates
		water molecules from the mol object, protonates them, adds coordinate
		information back, and then combines with final mol object. Note: if 
		any solvent molecules already have bound H, they will be skipped. 
		
		mol : rdkit mol object, must be defined.
		solvent : List of strings containing names of solvent residues in PDB.
			default is solvent=['HOH','WAT','T3P','SOL'] 
		
		Returns mol object. 
		'''
		solvent_mol = Chem.RWMol()
		no_solvent_mol = Chem.RWMol(mol)
		pos = []
		for atom in reversed(mol.GetAtoms()):
			if self.define_residue(atom)[1] in solvent:
				if self.res_info(atom,info='atom_name') == ' O  ':
					if len(atom.GetNeighbors()) != 2:
						pos.append(self.atom_coords(mol,atom))
						no_solvent_mol.RemoveAtom(atom.GetIdx())
						solvent_mol.AddAtom(atom)
		if solvent_mol.GetNumAtoms() > 0:
			rdDistGeom.EmbedMolecule(solvent_mol)
			for i,atom in enumerate(reversed(solvent_mol.GetAtoms())):
				solvent_mol.GetConformer().SetAtomPosition(atom.GetIdx(),pos[i])
				atom.SetHybridization(Chem.HybridizationType.SP3)
		
			solvent_mol = Chem.rdmolops.AddHs(solvent_mol,
											  addCoords=True,
											  addResidueInfo=True)

			combined_mol = Chem.CombineMols(no_solvent_mol,solvent_mol)
			return combined_mol
		else:
			return mol

###############################################################################
	def truncate(self, scheme='CA_terminal', output_file=None,
				 skip_residues=solvent_list, skip_resnumbers=[],
				 remove_resnumbers=[], remove_atom_ids=[], 
				 remove_sidechains=[], keep_backbones=[], 
				 constrain_atoms=[' CA '], add_hydrogens=False, 
				 exclude_solvent=False):
		'''
		Function to prepare truncated QMzyme model.
		scheme			- string to define what truncation scheme to use.
					Currently only option is 'CA_terminal'. This will remove 
					backbone atoms of residues that are not bondedto other 
					residues. I.e., if you have the following residues in 
					the active site (ALA120, GLY121, ASP200), the N-terminus 
					of ALA120 will be removed, the C-terminus of GLY121 will 
					be removed, and both termini will be removed from ASP 
					200 like a typical methyl-capping scheme. 
		output_file 		- string to define the output .pdb file name. 
		skip_residues		- list containing strings of residue names you 
					do not want to truncate at all.
		skip_resnumbers 	- list containing integers of residue numbers you
					do not want to truncate.
		remove_resnumbers	- list containing integers of residue numbers you
					want completely removed from the model.
		remove_atom_ids		- list containing integers of atom IDs you want
					completed removed from the model.
		remove_sidechains	- list containing strings of residue names you
					want the sidechains removed of. This can be useful if its 
					a bulky residue pointing away from the active site, and  
                    only its backbone is important to consider in the model. 
		keep_backbones		- list containing strings of residue names that
					you want to keep backbone atoms of, regardless of the 
                    capping scheme.
		constrain_atoms		- list containing strings of the atom names 
					(according to PDB format) that you want to constrain in the 
					calculation. 
		add_hydrogens		- boolean, default=False. If True, hydrogens will
					be added to your model with the reduce package. 
                    Recommdended if you started with a .pdb file that did not 
                    contain hydrogens, but you should still evaluate the 
                    final structure.
		exclude_solvent		- boolean, default=False. If True, solvent 
					molecules will be removed from the model.

		Generates the truncated_active_site_mol and constrain_atom_list
		attributes, and saves new .pdb file. 
		'''

		### Heidi's to do: add CA capping scheme, create capping summary as a
		# returned item, allow capping summary to be a function input so users
		# can specify the cap they want for each residue.

		if add_hydrogens is True:
			Chem.MolToPDBFile(self.active_site_mol,'temp.pdb')
			self.active_site_mol = self.add_H(pdb_file='temp.pdb')
			self.active_site_mol = self.protonate_solvent(mol=self.active_site_mol)

		new_mol = Chem.RWMol(self.active_site_mol)
		proline_count,bb_atom_count = 0,0
		constrain_list, N_terminus, C_terminus=[],[],[]
		previous_res = None
		remove_ids = []
		residues = []
		fix_N = []

		for atom in reversed(self.active_site_mol.GetAtoms()):
			current_res = self.define_residue(atom)
			if current_res==self.catalytic_center_definition:
				continue
			if atom.GetIdx() in remove_atom_ids:
				remove_ids.append(atom.GetIdx())
				continue
			if current_res[1] in skip_residues:
				continue
			if current_res[2] in skip_resnumbers:
				continue
			if current_res[2] in remove_resnumbers:
				remove_ids.append(atom.GetIdx())
				continue
			if current_res != previous_res:
				atom_name = self.res_info(atom,'atom_name')
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
					C_atom = self.active_site_mol.GetAtomWithIdx(C_id)
					N_atom = self.active_site_mol.GetAtomWithIdx(N_id)
					C_bonds = [self.res_info(x,'atom_name') \
							   for x in C_atom.GetNeighbors()]
					C_bonds_atoms = [x for x in C_atom.GetNeighbors()]
					N_bonds = [self.res_info(x,'atom_name') \
							   for x in N_atom.GetNeighbors()]
					N_bonds_atoms = [x for x in N_atom.GetNeighbors()]
				if bb_atom_count !=3:
					continue

			# CA_all capping scheme
			if scheme == 'CA_all':
				print("The C-alpha only capping scheme is currently"+\
					  " under development.")

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
							if 'H' in self.res_info(N_bonds_atoms[a],\
                                          'atom_name'):
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
						if 'H' in self.res_info(C_bonds_atoms[a],'atom_name'):
							remove_ids.append(C_bonds_atoms[a].GetIdx())
				bb_atom_count = 0
				previous_res = current_res

#TO DO: automate code to remove unecessary side chains
		if len(remove_sidechains) > 0:
			for atom in reversed(self.active_site_mol.GetAtoms()):
				current_res = self.define_residue(atom)
				if current_res[2] in remove_sidechains:
					atom_name = self.res_info(atom, 'atom_name')
					if atom_name in (' CA ',' C  ',' O  ',' N  '):
						continue
					if atom_name == ' CB ':
						cap_atom = new_mol.GetAtomWithIdx(atom.GetIdx())
						cap_atom.SetAtomicNum(1)
						cap_atom.GetPDBResidueInfo().SetName(' H* ')
						continue
					if atom_name == ' H  ':
						continue
					if ' CA ' not in [self.res_info(x,'atom_name') for x in \
					   atom.GetNeighbors()]:
						remove_ids.append(atom.GetIdx())
				else:
					continue

		for atom in reversed(new_mol.GetAtoms()):
			if self.res_info(atom,'atom_name') == ' CA ':
				bound_atoms = []
				lone_methyl = True
				for x in atom.GetNeighbors():
					bound_atoms.append(x.GetIdx())
					if 'H' not in self.res_info(x,'atom_name'):
						lone_methyl = False
				if lone_methyl is True:
					remove_ids.append(atom.GetIdx())
					for a in bound_atoms:
						remove_ids.append(a)
			if exclude_solvent is True:
				if self.define_residue(atom)[1] in solvent_list:
					remove_ids.append(atom.GetIdx())

		# Now actually create the truncated and capped mol object
		for a in reversed(np.unique(remove_ids)):
			new_mol.RemoveAtom(int(a))

		new_mol.UpdatePropertyCache(strict=False)
		Chem.SanitizeMol(new_mol,Chem.SanitizeFlags.SANITIZE_FINDRADICALS\
                         |Chem.SanitizeFlags.SANITIZE_KEKULIZE\
                         |Chem.SanitizeFlags.SANITIZE_SETAROMATICITY\
                         |Chem.SanitizeFlags.SANITIZE_SETCONJUGATION\
                         |Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION\
                         |Chem.SanitizeFlags.SANITIZE_SYMMRINGS,\
                         catchErrors=True)
		for atom in new_mol.GetAtoms():
			current_res = self.define_residue(atom)
			if current_res not in residues:
				residues.append(current_res)
			if ' H* ' in self.res_info(atom,info='atom_name'):
				if rdMolTransforms.GetBondLength(new_mol.GetConformer(),
				   atom.GetNeighbors()[0].GetIdx(), atom.GetIdx()) > 1.01:
					# This changes cap H bond length to something more physical
					rdMolTransforms.SetBondLength(new_mol.GetConformer(),
					atom.GetNeighbors()[0].GetIdx(),atom.GetIdx(), 0.970)
			if current_res[1] == 'PRO':
				if self.res_info(atom,info='atom_name') == ' N  ':
					atom.SetHybridization(Chem.HybridizationType.SP2)
					fix_N.append(atom.GetIdx())
			# Record atom ids to add to constrain list
			if self.res_info(atom,'atom_name') in constrain_atoms:
				constrain_list.append(atom.GetIdx()+1)

		if len(fix_N) > 0:
			new_mol = Chem.rdmolops.AddHs(new_mol,addCoords=True,onlyOnAtoms=fix_N)
		print("Final active site model contains {} atoms."
			  .format(new_mol.GetNumAtoms()))

		if output_file is None:
			output_file = '{}_{}_truncated_active_site_distance_cutoff_{}.pdb'\
						  .format(self.protein_prefix,self.cat_center,self.distance_cutoff)
		print("Structure saved as {}".format(output_file))

		Chem.MolToPDBFile(new_mol,output_file)

		if proline_count > 0:
			print("WARNING: Active site model contains {}".\
                  format(proline_count)+" proline(s). The N backbone"
                  " atom was automatically kept.")

		self.model_atom_count = new_mol.GetNumAtoms()
		self.truncated_active_site_mol = new_mol
		self.constrain_atom_list = constrain_list
		self.model_residues = residues
		print(delimeter)
		#return new_mol, constrain_list

###############################################################################
	def show_mol(self,mol=None):
		if mol is None:
			mol = self.catalytic_center_mol
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
