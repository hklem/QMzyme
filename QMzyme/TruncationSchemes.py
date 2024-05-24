###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
Module containing functions to truncate a QMzymeRegion based on some logic/workflow.
"""

from QMzyme.data import protein_residues
from QMzyme.RegionBuilder import RegionBuilder
from QMzyme.truncation_utils import *
import abc

class TruncationScheme(abc.ABC):
    def __init__(self, region, name):
        self.region = region
        if name == None:
            name = f'{self.region.name}_truncated'
        self.name = name
        self.truncate()
        self.return_region()

    @abc.abstractmethod
    def truncate(self):
        ...

    def return_region(self):
        self.region.rename(self.name)
        return self.region

class CA_terminal(TruncationScheme):
    def __init__(self, region, name):
        super().__init__(region, name)

    def truncate(self):
        remove_atoms = []
        region_builder = RegionBuilder(name=self.name)
        # does the region need to be in increasing resid order?
        for res in self.region.residues:
            resname = res.resname
            if resname not in protein_residues:
                continue
            # Define necessary backbone atoms
            Natom = res.get_atom('N')
            CAatom = res.get_atom('CA')
            Catom = res.get_atom('C')
            Oatom = res.get_atom('O')
            preceding_Catom = get_preceding_Catom(self.region, res.resid)
            following_Natom = get_following_Natom(self.region, res.resid)
            if resname != 'PRO':
                Hatom = res.get_atom('H')
            if preceding_Catom is not None and preceding_Catom.id not in self.region.ids:
                if resname != 'PRO':
                    region_builder.init_atom(cap_H(Natom, CAatom))
                    remove_atoms.append(Natom)
                    remove_atoms.append(Hatom)
                if resname == 'PRO':
                    region_builder.init_atom(cap_H(preceding_Catom, Natom))
            if following_Natom is not None and following_Natom.id not in self.region.ids:
                region_builder.init_atom(cap_H(Catom, CAatom))
                remove_atoms.append(Catom)
                remove_atoms.append(Oatom)
        for atom in self.region.atoms:
            if atom not in remove_atoms:
                region_builder.init_atom(atom)

        self.region = region_builder.get_region()

class CA_all(TruncationScheme):
    """
    Function to truncate a QMzymeRegion accoring to the CA_all scheme. 
    See Documentation for scheme descriptions.
    *I still need to write up some documentation about how various schemes works."
    """
    def __init__(self, region, name):
        super().__init__(region, name)

    def truncate(self):
        raise UserWarning("This method is currently under development.")  

# def CA_terminal(region): 
#     """
#     Function to truncate a QMzymeRegion accoring to the CA_terminal scheme. 
#     See Documentation for scheme descriptions.
#     *I still need to write up some documentation about how various schemes works."
#     """
#     # u = region.get_AtomGroup().universe
#     # residues = [MDAnalysisWrapper.get_parallel_residue(res, u) for res in region.get_AtomGroup().residues]
#     remove_atoms = []
#     add_atoms = []
#     cap_bond = [] # list of tuples of atoms. first atom in tuple is where cap goes, second atom is reference pt for bond length calc
    
#     region_builder = RegionBuilder(f'{region.name}_truncated')

#     #need to make sure region is in increasing resid order?
#     for res in region.residues:
#         resname = res.resname
#         if resname not in protein_residues:
#             continue
        
#         # Define necessary backbone atoms
#         Natom = res.get_atom('N')
#         CAatom = res.get_atom('CA')
#         Catom = res.get_atom('C')
#         Oatom = res.get_atom('O')
#         preceding_Catom = get_preceding_Catom(region, res.resid)
#         following_Natom = get_following_Natom(region, res.resid)
#         if resname != 'PRO':
#             Hatom = res.get_atom('H')

#         if preceding_Catom is not None and preceding_Catom.id not in region.ids:
#             if resname != 'PRO':
#                 #add_atoms.append(cap_H(Natom, CAatom))
#                 region_builder.init_atom(cap_H(Natom, CAatom))
#                 remove_atoms.append(Natom)
#                 remove_atoms.append(Hatom)
#             if resname == 'PRO':
#                 #add_atoms.append(cap_H(preceding_Catom, Natom))
#                 region_builder.init_atom(cap_H(preceding_Catom, Natom))

#         if following_Natom is not None and following_Natom.id not in region.ids:
#             #add_atoms.append(cap_H(Catom, CAatom))
#             region_builder.init_atom(cap_H(Catom, CAatom))
#             remove_atoms.append(Catom)
#             remove_atoms.append(Oatom)

#     for atom in region.atoms:
#         if atom not in remove_atoms:
#             #new_region.init_atom(atom)
#             region_builder.init_atom(atom)
#     # for atom in add_atoms:
#     #     new_region.init_atom(atom)
#     truncated_region = region_builder.get_region()

#     return truncated_region
    

# def CA_all(region):
#     """
#     Function to truncate a QMzymeRegion accoring to the CA_all scheme. 
#     See Documentation for scheme descriptions.
#     *I still need to write up some documentation about how various schemes works."
#     """
#     raise UserWarning("This method is currently under development.")        
        