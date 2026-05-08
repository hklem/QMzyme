###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
Module containing functions to truncate a QMzymeRegion based on some logic/scheme.
"""

from QMzyme.data import protein_residues, backbone_atoms
from QMzyme.QMzymeRegion import QMzymeRegion
from QMzyme.truncation_utils import *
import QMzyme.MDAnalysisWrapper as MDAwrapper
import abc


class TruncationScheme(abc.ABC):
    def __init__(self, region, name):
        self.region = region
        self.truncated_region = None
        if name == None:
            name = f'{self.region.name}_truncated'
        self.name = name
        self.truncate()
        if hasattr(self.region.atoms[0], "charge"):
            balance_charge(self.region, self.truncated_region)
        elif hasattr(self.region, "charge"):
            self.truncated_region.set_charge(self.region.charge)
        if getattr(self.region, "method") != None:
            self.truncated_region.set_method(self.region.method)
        self.return_region()

    @abc.abstractmethod
    def truncate(self):
        ...

    def return_region(self):
        self.truncated_region.rename(self.name)
        return self.truncated_region

class TerminalAlphaCarbon(TruncationScheme):
    """
    The TerminalAlphaCarbon scheme will 1) remove N-terminal backbone atoms 
    (N and H) if the preceding sequence residue (resid-1) is not included
    in the region and add a hydrogen atom along the CA–N backbone
    bond vector; and 2) remove C-terminal backbone atoms (C and O) if the 
    following sequence residue (resid+1) is not included in the region and
    add a hydrogen atom along the CA–C backbone bond vector. In
    the case of Proline, if the preceding sequence residue is not present 
    the Proline N atom is kept and a hydrogen is added along the N–(resid-1)C
    backbone bond vector.

    .. image:: ../../docs/Images/terminal_alpha_carbon.png
        :width: 250

    Image modified from Klem, H., McCullagh, M. & Paton, R. S. Top Catal. 
    65, 165–186 (2022). 
    """
    def __init__(self, region, name):
        super().__init__(region, name)

    def truncate(self):
        remove_atoms = []
        r = QMzymeRegion(name = self.name, atoms = self.region.atoms, universe = self.region._universe)
        for res in self.region.residues:
            resname = res.resname
            if resname not in protein_residues:
                continue
            # Define necessary backbone atoms
            Natom = res.get_atom(backbone_atoms['N'])
            CAatom = res.get_atom(backbone_atoms['CA'])
            Catom = res.get_atom(backbone_atoms['C'])
            Oatom = res.get_atom(backbone_atoms['O'])
            preceding_Catom = get_preceding_Catom(self.region, res.resid)
            following_Natom = get_following_Natom(self.region, res.resid)
            if resname != 'PRO':
                Hatom = res.get_atom(backbone_atoms['H'])
            if preceding_Catom is not None and preceding_Catom.id not in self.region.ids:
                if resname != 'PRO':
                    cap_atom = cap_H(Natom, CAatom)
                    r.remove_atom(r.get_atom(id=Natom.id))
                    r.remove_atom(r.get_atom(id=Hatom.id))
                    r.add_atom(cap_atom)
                if resname == 'PRO':
                    cap_atom = cap_H(preceding_Catom, Natom)
                    setattr(cap_atom, "id", cap_atom.id-1)
                    r.add_atom(cap_atom)
            if following_Natom is not None and following_Natom.id not in self.region.ids:
                cap_atom = cap_H(Catom, CAatom)
                r.remove_atom(r.get_atom(id=Catom.id))
                r.remove_atom(r.get_atom(id=Oatom.id))
                r.add_atom(cap_atom)
        # if hasattr(self.region, "charge"):
        #     r.set_charge(self.region.charge)
        # if getattr(self.region, "method") != None:
        #     r.set_method(self.region.method)
        self.truncated_region = r
        
class AlphaCarbon(TruncationScheme):
    """
    Function to truncate a QMzymeRegion accoring to the AlphaCarbon scheme. 
    This method is still under development. 

    .. image:: ../../docs/Images/all_alpha_carbon.png
        :width: 250

    Image modified from Klem, H., McCullagh, M. & Paton, R. S. Top Catal. 
    65, 165–186 (2022). 
    """
    def __init__(self, region, name):
        super().__init__(region, name)

    def truncate(self):
        remove_atoms = []
        r = QMzymeRegion(name=self.name, atoms=self.region.atoms, universe=self.region._universe)
        for res in self.region.residues:
            resname = res.resname
            if resname not in protein_residues:
                continue
            # Define necessary backbone atoms
            Natom = res.get_atom(backbone_atoms['N'])
            CAatom = res.get_atom(backbone_atoms['CA'])
            Catom = res.get_atom(backbone_atoms['C'])
            Oatom = res.get_atom(backbone_atoms['O'])
            preceding_Catom = get_preceding_Catom(self.region, res.resid)
            following_Natom = get_following_Natom(self.region, res.resid)
            if preceding_Catom is not None: # fixes issues if this is the very first res in sequence
                if resname != 'PRO':
                    Hatom = res.get_atom(backbone_atoms['H'])
                    cap_atom = cap_H(Natom, CAatom)
                    r.remove_atom(r.get_atom(id=Natom.id))
                    r.remove_atom(r.get_atom(id=Hatom.id))
                    r.add_atom(cap_atom)
                if resname == 'PRO':
                    cap_atom = cap_H(preceding_Catom, Natom)
                    setattr(cap_atom, "id", cap_atom.id-1)
                    r.add_atom(cap_atom)
            if following_Natom is not None: # fixes issues if this is the very last res in sequence
                cap_atom = cap_H(Catom, CAatom)
                r.remove_atom(r.get_atom(id=Catom.id))
                r.remove_atom(r.get_atom(id=Oatom.id))
                r.add_atom(cap_atom)
        self.truncated_region = r

class BetaCarbon(TruncationScheme):
    """
    The Beta Carbon scheme will 1) select for atoms that are within 2Å from CB;
    2) remove all non-backbone atoms that are outside of 2Å distance, and
    3) remove non-hydrogen and non-backbone atoms and replace it with hydrogen
    along the CB-X vector. In the case of Proline and Glycine, it skips and returns
    a warning message.
    """
    def __init__(self, region, ala_atom_group, model, name):
        self.ala_atom_group = ala_atom_group
        self.region = region
        self.model = model
        super().__init__(region, name)

    def truncate(self):
        remove_atoms = []
        r = self.region

        # Since the resid that is needed to be truncated is written in MDA, we need to digest it
        residues_to_truncate = set(a.resid for a in self.ala_atom_group.residues)

        # Iterating over the whole residues
        for res in r.residues:
            # This is needed to select for specific residues that are in ala_atom_group
            if res.resid not in residues_to_truncate:
                continue

            # Raise warning if it contains Gly and Pro within alanine_mutation
            if res.resname == "GLY" or res.resname == "PRO":
                UserWarning("Pro and Gly exists within alanine_mutation. Please remove the residue.")
                continue
            
            # Get residue selection string in the universe
            sel_str = f"resid {res.resid}"

            # Make AtomGroups in the original MDAnalysis universe
            res_atoms = self.model.universe.select_atoms(sel_str)

            # Define necessary backbone atoms and CB
            CBatom = res.get_atom('CB')

            # MDAwrapper can only digest universe, so we convert it
            CB_sel = self.region._universe.select_atoms(f"resid {CBatom.resid} and name {CBatom.name}")

            # Getting the neighbor atoms
            neighbors = MDAwrapper.get_neighbors(res_atoms,CB_sel, 2)

            # Selecting the non neighbors to remove all of them
            non_neighbors = res_atoms.atoms - neighbors.atoms

            keep_names = set(backbone_atoms.values())
            keep_names.update(["HB1", "HB2", "HB3", CBatom.name])
            
            # Remove all atoms that are not neighbors and are not backbone
            for mda_atom in non_neighbors.atoms:
                for qm_atom in res.atoms:
                    # Since we are now going from MD universe to QMzyme region, we need to check
                    if qm_atom.name != mda_atom.name or qm_atom.resid != mda_atom.resid:
                        continue
                    if qm_atom.name in keep_names:
                        continue
                    res.remove_atom(qm_atom)

            # Adding H atoms and removing connected atoms
            for mda_atom in neighbors.atoms:
                for qm_atom in res.atoms:
                    if qm_atom.name != mda_atom.name or qm_atom.resid != mda_atom.resid:
                        continue
                    if qm_atom.name in keep_names:
                        continue

                    # Replacing neighbor atoms with H
                    cap_atom = cap_H(qm_atom,CBatom)
                    r.remove_atom(qm_atom)
                    r.add_atom(cap_atom)

        self.truncated_region = r
