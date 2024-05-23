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
                if 'QM2' in CalculateModel.calculation:
                    raise UserWarning("Two QM regions have already been defined."
                                      "QMzyme currently only supports a maximum "
                                      "of two QM regions in a single calculation.")
                calc = 'QM2'
                region.method["type"] = 'QM2'
                print(f"WARNING: A QM region already defined. Setting QMzymeRegion {region.name} as QM2.")
        
        CalculateModel.calculation[calc] = region
    
        if len(CalculateModel.calculation) > 1:
            self._check_overlap(region)
    
    def _check_overlap(self, region):
        if 'QM' in CalculateModel.calculation:    
            high_region = CalculateModel.calculation['QM']
            low_region = region
        elif region.method["type"] == 'QM':
            high_region = region
            try:
                low_region = CalculateModel.calculation['XTB']
            except:
                low_region = CalculateModel.calculation['ChargeField']

        common_atoms = high_region.get_overlap(low_region)
        if len(common_atoms) != 0:
            #residues = [atom.resname+str(atom.resid) for atom in common_atoms]
            residues = [atom.resid for atom in common_atoms]
            residues = [high_region.get_residue(resid) for resid in list(set(residues))]
            print(f"\nWARNING: Region overlap detected. The following residue(s) were found in both regions: {list(set(residues))}.")
            print(f"Removing duplicate atoms from {low_region.name} and recalculating charge.")
            subtracted = low_region.subtract(high_region)
            subtracted.guess_charge()
            subtracted.method = low_region.method
            subtracted.method["charge"] = subtracted.charge
            subtracted.method["freeze_atoms"] = subtracted.get_indices(attribute='is_fixed', value=True)
            subtracted.name = low_region.name
            CalculateModel.calculation[low_region.method["type"]] = subtracted


    def _reset():
        CalculateModel.calculation = {}
    
class CalculationBase:
    def assign_to_region(self, region, charge=None, mult=1):
        self._set_constraints(region)
        self.mult = mult
        region.set_method(self.__dict__)
        region = CalculateModel()._add(calc=self.type, region=region)
        self._set_charge(CalculateModel.calculation[self.type], charge)

    def _set_charge(self, region, charge):
        if charge is None:
            if not hasattr(region, "charge"):
                region.guess_charge()
            self.charge = region.charge
        else:
            self.charge = charge

    def _set_constraints(self, region):
        self.freeze_atoms = region.get_indices('is_fixed', True)

class QM_Method(CalculationBase):
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
    def __init__(self, basis_set, functional, qm_input="", qm_end="", program='orca'):
        """
        :param region: QMzyme region to treat at the QM level.
        :type region: QMzymeRegion
        """
        self.type = 'QM'
        self.qm_input = qm_input
        self.basis_set = basis_set
        self.functional = functional
        self.qm_input = qm_input
        self.qm_end = qm_end
        self.program = program
        self._set_qm_input()

    def _set_qm_input(self):
        for info in [self.functional, self.basis_set]:
            if info not in self.qm_input:
                self.qm_input = f"{info} {self.qm_input}"
        self.qm_input.strip()

    # def _set_charge(self, region, charge):
    #     if charge is None:
    #         if not hasattr(region, "charge"):
    #             region.guess_charge()
    #         self.charge = region.charge
    #     else:
    #         self.charge = charge

    # def _set_constraints(self, region):
    #     self.freeze_atoms = region.get_indices('is_fixed', True)

    # def assign_to_region(self, region, charge=None, mult=1):
    #     self._set_constraints(region)
    #     self._set_charge(region, charge)
    #     self.mult = mult
    #     region.set_method(self.__dict__, _type="QM")
    #     CalculateModel()._add(calc='QM', region=region)
    
class XTB_Method(CalculationBase):
    """
    Class to prepare a QMzymeRegion for xTB treatment.
    """
    def __init__(self):
        self.type = 'XTB'

class ChargeField:
    """
    Class to prepare a QMzymeRegion for ChargeField treatment.
    """
    def __init__(self):
        self.type = 'ChargeField'
    