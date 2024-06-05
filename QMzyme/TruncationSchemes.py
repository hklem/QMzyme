###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
Module containing functions to truncate a QMzymeRegion based on some logic/workflow.
"""

from QMzyme.data import protein_residues
from QMzyme.QMzymeRegion import QMzymeRegion
from QMzyme.truncation_utils import *
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
        if hasattr(self.region, "charge"):
            r.set_charge(self.region.charge)
        if getattr(self.region, "method") != None:
            r.set_method(self.region.method)
        self.truncated_region = r
        
class AlphaCarbon(TruncationScheme):
    """
    Function to truncate a QMzymeRegion accoring to the AlphaCarbon scheme. 
    See Documentation for scheme descriptions.
    *I still need to write up some documentation about how various schemes works."
    """
    def __init__(self, region, name):
        super().__init__(region, name)

    def truncate(self):
        raise UserWarning("This method is currently under development.")  
