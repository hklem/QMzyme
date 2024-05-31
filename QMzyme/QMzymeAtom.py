###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

import inspect

class QMzymeAtom:
    """
    Required Parameters
    --------------------
    :param id: Atom id. If the atom is generated from an MDAnalysis AtomGroup (probably) the id will be the universe atom ID. 
    :param name: Atom name: ex., 'C1'
    :param element: Element one-letter name: ex., 'C'
    :param position: Array of cartesian coordinates.
    :param resid: Integer residue number.
    :param resname: Three letter residue name: ex., 'VAL'
    """
    def __init__(self, name, element, position, resid, resname, id=1, region=None, **kwargs):
        self.name = name
        self.element = element
        self.position = position
        self.resid = resid
        self.resname = resname
        self.id = id
        self.__region = region
        self.is_fixed = False
        self.is_point_charge = False
        
        for attr, val in kwargs. items():
            setattr(self, attr, val)

    def __repr__(self):
        return f"<QMzymeAtom {self.id}: {self.name} of resname {self.resname}, resid {self.resid}>"

    @property
    def region(self):
        """
        Returns the QMzymeRegion this atom belongs to.
        """
        return self.__region
    
    @region.setter
    def region(self, value):
        fname = inspect.currentframe().f_code.co_name
        raise AttributeError(f"This attribute is protected. If you truly wish to change its value use self.set_{fname}({value}).")

    def _set_region(self, value):
        self.__region = value

    def set_neighbor(self, value: bool=True):
        """
        Method to set ``is_neighbor=True`` for QMzymeAtom instance."
        """
        self.is_neighbor = value
    
    def set_fixed(self, value: bool=True):
        """
        Method to set ``is_fixed=True`` for QMzymeAtom instance."
        """
        self.is_fixed = value

    def set_point_charge(self, value: bool=True):
        """
        Method to set ``is_point_charge=True`` for QMzymeAtom instance.
        Note: only works if atom has charge attribute.
        """
        if not hasattr(self, "charge"):
            raise UserWarning(f"Cannot set atom {self} as point_charge because no charge information found.")
        self.is_point_charge = value

    def get_chain(self):
        """
        Searches possible attribute names chain can live under 
        (chainID, chain_ID, chain,id) and returns value if found.
        """
        chain = None
        for name in ['chain', 'chainID', 'chain_ID', 'chain_id', 'chainid']:
            if hasattr(self, name):
                chain = getattr(self, name)
        return chain
    
    def is_within(self, region):
        """
        Returns True if the same atom is found in Region. Used to avoid duplication.
        """
        atom = region.get_atom(id=self.id)
        if atom is None:
            return False
        for k in ['name', 'resid', 'resname', 'element']:
            if self.__dict__[k] != atom.__dict__[k]:
                return False
        return True
