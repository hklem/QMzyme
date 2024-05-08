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
    :param is_neighbor: optional
    """
    def __init__(self, name, element, position, resid, resname, id=1, region=None, **kwargs):
        self.__name = name
        self.__element = element
        self.position = position
        self.resid = resid
        self.resname = resname
        self.id = id
        self.__region = region
        
        for attr, val in kwargs. items():
            setattr(self, attr, val)

    def __repr__(self):
        return f"<QMzymeAtom {self.id}: {self.name} of resname {self.resname}, resid {self.resid}>"

    @property 
    def name(self):
        return self.__name

    @name.setter
    def name(self, value):
        fname = inspect.currentframe().f_code.co_name
        raise AttributeError(f"This attribute is protected. If you truly wish to change its value use self.set_{fname}(value).")
        
    def set_name(self, value):
        self.__name = value
    
    @property
    def element(self):
        return self.__element
    
    @element.setter
    def element(self, value):
        fname = inspect.currentframe().f_code.co_name
        raise AttributeError(f"This attribute is protected. If you truly wish to change its value use self.set_{fname}(value).")
        
    def set_element(self, value):
        self.__element = value

    def set_neighbor(self, value: bool):
        self.neighbor = value

    @property
    def region(self):
        return self.__region
    
    @name.setter
    def region(self, value):
        fname = inspect.currentframe().f_code.co_name
        raise AttributeError(f"This attribute is protected. If you truly wish to change its value use self.set_{fname}(value).")

    def set_region(self, value):
        self.__region = value
    
    def set_fixed(self, value: bool=True):
        """
        Method to add ``is_fixed=True`` attribute to QMzymeAtom."
        """
        self.is_fixed = value

    def set_dummy(self, value: bool=True):
        """
        Method to add ``is_dummy=True`` attribute to QMzymeAtom."
        """
        self.is_dummy = value

    def get_chain(self):
        chain = None
        for name in ['chain', 'chainID', 'chain_ID', 'chain_id']:
            if hasattr(self, name):
                chain = getattr(self, name)
        return chain