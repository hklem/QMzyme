from QMzyme.utils import protein_residues, set_bond_length
from QMzyme.RegionBuilder import RegionBuilder
import copy

truncation_schemes = {
    'CA_terminal': lambda region: CA_terminal(region),
    'CA_all': lambda region: CA_all(region),
}


#Fix bug: code currently overwrites the QMzyme region.
def CA_terminal(og_region): 
    # u = region.get_AtomGroup().universe
    # residues = [MDAnalysisWrapper.get_parallel_residue(res, u) for res in region.get_AtomGroup().residues]
    region = copy.copy(og_region)
    remove_atoms = []
    cap_bond = [] # list of tuples of atoms. first atom in tuple is where cap goes, second atom is reference pt for bond length calc
    
    #need to make sure region is in increasing resid order?

    for resid in region.resids:
        resname = region.get_resname(resid)
        if resname not in protein_residues:
            continue
        Natom = copy.copy(region.get_residue_atom(resid, 'N'))
        CAatom = copy.copy(region.get_residue_atom(resid, 'CA'))
        Catom = copy.copy(region.get_residue_atom(resid, 'C'))
        Oatom = copy.copy(region.get_residue_atom(resid, 'O'))
        preceding_Catom = copy.copy(get_preceding_Catom(region, resid))
        following_Natom = copy.copy(get_following_Natom(region, resid))
        if resname != 'PRO':
            Hatom = copy.copy(region.get_residue_atom(resid, 'H'))
        # Remove N terminus of first residue in region
        if preceding_Catom.id not in region.ids:
            if resname != 'PRO':
                cap_bond.append((Natom, CAatom))
                #remove_atoms.append(Natom)
                remove_atoms.append(Hatom)
            else:
                cap_bond.append((preceding_Catom, Natom))
        if resid == region.resids[0]:
            if resname != 'PRO':
                cap_bond.append((Natom, CAatom))
                #remove_atoms.append(Natom)
                remove_atoms.append(Hatom)
            if resname == 'PRO':
                cap_bond.append((preceding_Catom, Natom))
            continue
        if following_Natom.id not in region.ids:
            cap_bond.append((Catom, CAatom))
            #remove_atoms.append(Catom)
            remove_atoms.append(Oatom)
        
        for bond in cap_bond:
            cap_H(bond[0], bond[1])
        
        new_region = RegionBuilder(f'{region.name}_truncated')
        for atom in region.atoms:
            if atom not in remove_atoms:
                new_region.init_atom(atom)

        return new_region.region
    

def CA_all(region):
    raise UserWarning("This method is currently under development.")        
        

def cap_H(replace_atom, fixed_atom, bond_length=1.09):
    """
    :param replace_atom: The Atom to be converted to hydrogen.
    :param fixed_atom: The Atom bound to replace_atom that serves as reference point for bond vector calculation.
    """
    new_position = set_bond_length(replace_atom.position, fixed_atom.position, bond_length)
    new_name = name_cap_H(fixed_atom.region, fixed_atom.resid)
    new_atom_dict = {
            'element': 'H',
            'type': 'H',  
            'name': new_name, 
            'position': new_position, 
            'mass': 1.00794,
        }
    new_atom = alter_atom(replace_atom, new_atom_dict)
    if new_atom not in fixed_atom.region.atoms:
        fixed_atom.region.add_atom(new_atom)


def alter_atom(atom, new_atom_dict):
    for key, val in new_atom_dict.items():
        if hasattr(atom, key):
            setattr(atom, key, val)
    return atom


def name_cap_H(region, resid, name='H1'):
    residue_atoms = region.get_residue_atoms(resid)
    i = 1
    while name in [atom.name for atom in residue_atoms]:
        i += 1
        name = f'H{i}'
    return name


def get_preceding_Catom(region, resid):
    mda_atom = region.get_atom_group().universe.select_atoms(f'resid {resid-1} and name C')
    atom = RegionBuilder('temp', mda_atom).get_region().atoms[0]
    return atom


def get_following_Natom(region, resid):
    mda_atom = region.get_atom_group().universe.select_atoms(f'resid {resid+1} and name N')
    atom = RegionBuilder('temp', mda_atom).get_region().atoms[0]
    return atom


def has_Nterm_neighbor(atom):
    return atom.resid-1 in atom.region.resids


def has_Cterm_neighbor(atom):
    return atom.resid+1 in atom.region.resids



