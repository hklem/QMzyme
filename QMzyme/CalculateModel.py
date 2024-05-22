###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
Module in charge of creating input files for QM-only or QM/MM calculations. This 
module integrates the `AQME QPREP <https://aqme.readthedocs.io/en/latest/API/aqme.qprep.html>`_ 
workflow.

Notes
...............
    *   Currently optimized to generate QM-only Gaussian input files.
"""

import os
from QMzyme.aqme.qprep import qprep
import numpy as np


class CalculateModel:
    calculation = {}
    def _add(self, calc, region):
        if calc in CalculateModel.calculation:
            if calc == 'QM':
                calc = 'QM2'
                print(f"WARNING: A QM region already defined. Setting {region} as QM2.")
        CalculateModel.calculation[calc] = region

class QM:
    """
    Class to prepare a QMzymeRegion for QM treatment.
    
    Required Parameters
    ====================
    region: QMzymeRegion to apply method to
    basis_set: str, defines basis set to use for calculation
    functional: str, defines functional to use for calculation

    Optional Parameters 
    ====================
    charge: int, charge of the region. If not provided in parameters, charge will be guessed.
    mult: int, multiplicity of the region. default = 1.
    qm_input: str, default = ""
        Keywords to include in the input file route line to 
        declare any details beyond the basis set and functional.
        E.g. "EmpiricalDispersion=GD3BJ opt freq". Not including anything
        here means the calculation will be a single-point energy calculation.
    qm_end: str, default = ""
        Final line(s) in the input file

    """
    def __init__(self, region, basis_set, functional, charge=None, mult=1, qm_input="", qm_end="", program='orca'):
        """
        :param region: QMzyme region to treat at the QM level.
        :type region: QMzymeRegion
        """

        self.qm_input = qm_input
        self.basis_set = basis_set
        self.functional = functional
        if charge is None:
            if not hasattr(region, "charge"):
                region.guess_charge()
            self.charge = region.charge
        self.charge = charge
        self.mult = mult
        self.qm_input = qm_input
        self.qm_end = qm_end
        self.freeze_atoms = []
        self.program = program
        self.files = region.write()
    
        self._set_qm_input()
        self._set_constraints(region)
        #return self
        region.set_method(self.__dict__, _type="QM")
        CalculateModel()._add(calc='QM', region=region)
    
    def _set_constraints(self, region):
        self.freeze_atoms = region.get_indices('is_fixed', True)

    def _set_qm_input(self):
        for info in [self.functional, self.basis_set]:
            if info not in self.qm_input:
                self.qm_input = f"{info} {self.qm_input}"

class XTB:
    def __init__(self, region):
        CalculateModel()._add(calc='XTB', region=region)

class ChargeField:
    def __init__(self, region):
        CalculateModel()._add(calc='ChargeField', region=region)



        
class xTB:
    """
    Class to prepare a QMzymeRegion for xTB treatment.
    """
    def __init__(self, region, charge, mult):
        """
        :param region: QMzyme region to treat at the xTB level.
        :type region: QMzymeRegion
        :param charge: Charge of the region.
        :type charge: int
        :param mult: Multiplicity of the region.
        :type mult: int
        """
        self.charge = charge
        self.multiplicity = mult
        self.freeze_atoms = []
        
        self._set_constraints(self, region)
        region._set_method(self.__dict__, type="xTB")


class ChargeField:
    pass


    
