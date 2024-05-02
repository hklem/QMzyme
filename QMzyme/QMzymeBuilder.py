import warnings
import numpy as np
import copy

class QMzymeModel():
    def __repr__(self):
        return f"<QMzymeModel>"

    def set_region(self, region_name, AtomGroup, overwrite=False):
        region = QMzymeRegion(region_name, AtomGroup)
        #setattr(region, 'name', region_name)
        if hasattr(self, region_name):
            raise UserWarning(f"{self} already has a region with the name {region_name}.")
        setattr(self, region_name, region)

    def get_region(self, region_name):
        return getattr(self, region_name)


class QMzymeRegion:
    def __init__(self, name=None, AtomGroup=None):
        if name is not None:
            self.name = name

        if AtomGroup is not None:
            self.atoms = []
            self.residues = []
            #self._AtomGroup = AtomGroup
            self.__AtomGroup = AtomGroup
            for atom in AtomGroup.atoms:
                self.add_atom(atom)

    def __repr__(self):
        #return f"<QMzymeRegion {self.name} with {self.n_atoms} atoms>"
        return f"<QMzymeRegion {self.name} contains {self.n_atoms} atoms and {self.n_residues} residues>"
    
    @property
    def n_atoms(self):
        return len(self.atoms)
    
    @property
    def n_residues(self):
        self.residues = list(set([atom.resid for atom in self.atoms]))
        return len(self.residues)
    

    def init_atom_dict(self, atom):
        atom_attr_dict = {}
        for attr in dir(atom):
            if attr.startswith('_') or attr.startswith('get'):
                continue
            try:
                atom_attr_dict[attr] = getattr(atom, attr)
            except:
                pass
        return atom_attr_dict

    def add_atom(self, atom):
        """
        :param atom: The atom you want to add to the QMzymeRegion. 
        :type atom: QMzymeAtom or MDanalysis.Atom
        """
        warnings.filterwarnings('ignore')

        atom_attr_dict = self.init_atom_dict(atom)

        if str(type(atom)) == "<class 'MDAnalysis.core.groups.Atom'>":
            atom_attr_dict['chain'] = atom_attr_dict['chainID']
            del atom_attr_dict['chainID']

        qmz_atom = QMzymeAtom(**atom_attr_dict)

        # To make sure every atom has unique id
        if atom_attr_dict['id'] in self.ids:
            if self._is_unique(qmz_atom):
                new_id = atom_attr_dict['id']
                while atom_attr_dict['id'] in self.ids:
                    new_id += 1
                qmz_atom.id = new_id
                raise UserWarning(f"A QMzymeAtom with id {atom_attr_dict['id']} already exists in "+
                                f"QMzymeRegion {self.name}. Value has been automatically modified to "+
                                f"next available integer value: {new_id}.")
            else:
                raise UserWarning(f"{atom} already exists in QMzymeRegion {self.name}.")
            
        self.atoms.append(qmz_atom)

    @property
    def ids(self):
        return [atom.id for atom in self.atoms]
    
    def sort(self, key='id'):
        sorted_atoms = []
        original_atoms = copy.copy(self.atoms)
        for i in range(self.n_atoms):
            key_list = [getattr(atom, key) for atom in original_atoms]
            min = np.min(key_list)
            min_index = key_list.index(min)
            sorted_atoms.append(original_atoms[min_index])
            del original_atoms[min_index]
        self.atoms = sorted_atoms


    def _is_unique(self, atom):
        if atom in self.atoms:
            return False
        else:
            return True
        
    def add_attr(self, attr_name, attr_value):
        setattr(self, attr_name, attr_value)
    
    def get_AtomGroup(self):
        return self.__AtomGroup

class QMzymeAtom:
    """
    Required Parameters
    --------------------
    :param id: id starts from 1 NOT 0. 
    :param name:
    :param element:
    :param position:
    :param resid:
    :param resname:

    Optional Parameters
    --------------------
    :param type: optional
    :param is_fixed: optional
    :param force: optional
    :param velocity: optional
    :param id: optional
    :param bfactor: optional
    :param resid: 
    :param resname:
    :param segid: optional
    :param chain: optional
    :param record_type: optional
    :param charge: optional
    """

    def __init__(self, name, element, position, resid, resname, id=1, **kwargs):
        self.name = name
        self.element = element
        self.position = position
        self.resid = resid
        self.resname = resname
        self.id = id

        for attribute, value in kwargs.items():
            setattr(self, attribute, value)

    def __repr__(self):
        return f"<QMzymeAtom {self.id}: {self.name} of resname {self.resname}, resid {self.resid}>"
    
    def is_neighbor(self):
        setattr(self, 'is_neighbor', True)

