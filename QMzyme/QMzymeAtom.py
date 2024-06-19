###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

import inspect

class QMzymeAtom:
    """
    Required Parameters
    --------------------
    
    :param name: Atom name: ex., 'C1'
    :type name: str
    
    :param element: Element one-letter name: ex., 'C'
    :type element: str
    
    :param position: Array of cartesian coordinates.
    :type position: List[float] or np.array
    
    :param resid: Integer residue number.
    :type resid: int
    
    :param resname: Three letter residue name: ex., 'VAL'
    :type resname: str

    Pararmeters with defaults
    ---------------------------
    
    :param id: Atom id. If the atom is generated from an MDAnalysis AtomGroup (probably the case) 
        the id will be the universe atom ID. 
    :type id: int, default=1

    :param region: Region the atom exists within (if any).
    :type region: :class:`~QMzyme.QMzymeRegion.QMzymeRegion`, default=None

    .. note:: User may add any number of additional attributes as keyword arguments (kwargs).
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
        :returns: Region this atom belongs to. If it does not belong to a region, returns None.
        :rtype: :class`~QMzyme.QMzymeRegion.QMzymeRegion`
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
        Sets ``is_neighbor=True`` for QMzymeAtom instance, unless 'value=False' is passed.
        
        .. note:: The :class:`~QMzyme.SelectionSchemes.DistanceCutoff` class calls this method. This might
            be useful if you want to further evaluate why what residues were included in a QMzymeRegion that
            was selected based on the distance cutoff scheme.
        """
        self.is_neighbor = value
    
    def set_fixed(self, value: bool=True):
        """
        Sets ``is_fixed=True`` for QMzymeAtom instance, unless 'value=False' is passed.

        .. note:: The :module:`~QMzyme.CalculateModel` module will read what atoms have ``set_fixed=True`` and 
            use that information to communicate to :module:`~QMzyme.Writers` what atoms to constrain in the
            calculation input file.
        """
        self.is_fixed = value

    def set_point_charge(self, value: bool=True):
        """
        Sets ``is_point_charge=True`` for QMzymeAtom instance, unless 'value=False' is passed. This will eventually
        be used for calculations with charge embedding.

        :raises: UserWarning if the atom does not have the attribute 'charge'.
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
        :param region: Region to search for atom in.
        :type region: :class:`~QMzyme.QMzymeRegion.QMzymeRegion`, required
        :returns: True if the same atom is found in region. Used to avoid duplication.
        :rtype: bool
        """
        atom = region.get_atom(id=self.id)
        if atom is None:
            return False
        for k in ['name', 'resid', 'resname', 'element']:
            if self.__dict__[k] != atom.__dict__[k]:
                return False
        return True
