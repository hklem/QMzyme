###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

import inspect

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
    :param force: optional
    :param velocity: optional
    :param id: optional
    :param bfactor: optional
    :param segid: optional
    :param chain: optional
    :param record_type: optional
    :param charge: optional
    :param is_fixed: optional
    :param is_neighbor: optional
    :param is_point_charge: optional
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

    def set_name(self, value):
        self.name = value

    def set_element(self, value):
        self.element = value

    @property
    def region(self):
        return self.__region
    
    @region.setter
    def region(self, value):
        fname = inspect.currentframe().f_code.co_name
        raise AttributeError(f"This attribute is protected. If you truly wish to change its value use self.set_{fname}({value}).")

    def set_region(self, value):
        self.__region = value

    def set_neighbor(self, value: bool=True):
        """
        Method to add ``is_neighbor=True`` attribute to QMzymeAtom."
        """
        self.is_neighbor = value
    
    def set_fixed(self, value: bool=True):
        """
        Method to add ``is_fixed=True`` attribute to QMzymeAtom."
        """
        self.is_fixed = value

    def set_point_charge(self, value: bool=True):
        """
        Method to set ``is_point_charge=True`` to QMzymeAtom.
        Note: only works if atom has charge attribute.
        """
        if not hasattr(self, "charge"):
            raise UserWarning(f"Cannot set atom {self} as point_charge because no charge information found.")
        self.is_point_charge = value

    def get_chain(self):
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
