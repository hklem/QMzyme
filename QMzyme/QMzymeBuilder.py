from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.core.groups import Atom

required_atom_parameters = {
    'name': None,
    'element': None,
    'coordinates': None,
    'resid': None,
    'resname': None,
}

class QMzymeAtom(Atom):
    """
    Attributes
    ----------
    name:
    element:
    type:
    is_fixed: optional
    coordinates:
    force: optional
    velocity: optional
    ix:
    bfactor: optional
    resid: 
    resname:
    segid: optional
    chainID:
    record_type: optional
    charge: optional
    """

    # def __init__(self, kwargs):
    #     for key in required_atom_parameters.keys():
    #         try:
    #             setattr(self, key, kwargs[key])
    #             del kwargs[key]
    #         except:
    #             raise AttributeError(f"QMzymeAtom requires {key} parameter to be set.")
    #     for key in kwargs[key]:
    #         setattr(self, key, kwargs[key])

    @property
    def is_fixed(self):
        return self.is_fixed
    
    @is_fixed.setter
    def is_fixed(self, val):
        if type(val) != bool:
            raise UserWarning("`is_fixed` must be True or False.")
        setattr(self, 'is_fixed', val)


class QMzymeAtomGroup(AtomGroup):
    """
    An ordered array of atoms.

    ag = QMzymeAtromGroup([QMzymeAtom1, QMzymeAtom2, ...])

    """
