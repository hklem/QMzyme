
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
        self.name = name
        self.element = element
        self.position = position
        self.resid = resid
        self.resname = resname
        self.id = id
        self.region = region

        for attribute, value in kwargs.items():
            setattr(self, attribute, value)

    def __repr__(self):
        return f"<QMzymeAtom {self.id}: {self.name} of resname {self.resname}, resid {self.resid}>"

    def _add_prop(self, name, val, overwrite=False):
        if overwrite is False and hasattr(self, name):
            raise UserWarning(f"{self} already has property {name} set with a value of {val}. "+
                              "To overwrite this property set `overwrite=True`.")
        setattr(self, name, val)

    def set_neighbor(self, value: bool):
        self.is_neighbor = value

    # def is_neighbor(self):
    #     return self.is_neighbor


