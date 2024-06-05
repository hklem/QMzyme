###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
Module in charge of creating input files for QM-only or QM/MM calculations. This 
module integrates the `AQME QPREP <https://aqme.readthedocs.io/en/latest/API/aqme.qprep.html>`_ 
workflow.
"""

from QMzyme.data import protein_residues

class CalculateModel:
    """
    Class responsible for storing calculation type and associated region. 
    This information is used to determine what type of calculation writer to call 
    in the Writers module.

    This class also determines if there are any model issues, such as overlap between
    two regions, which may lead to incorrect charge estimations.
    """
    calculation = {}
    calc_type = None

    def add(self, type, region):
        if type == 'QM' and 'QM' in CalculateModel.calculation:
            # If region is same, we assume the user is trying to overwrite the old method
            # Otherwise, we will rename the type to allow for multiple QM regions 
            if region != CalculateModel.calculation['QM']:
                count = 2
                while type in CalculateModel.calculation:
                    type = f'QM{count}'
                    count += 1
        CalculateModel.calculation[type] = region
        CalculateModel.calc_type = type

    def combine_regions_and_methods():
        if len(CalculateModel.calculation) == 1:
            return
        # Start with QM region(s)
        # collect QM_regions:
        calc_type = ''
        qm_regions_sorted = sorted([calc for calc in CalculateModel.calculation if 'QM' in calc])
        base_region = CalculateModel.calculation[qm_regions_sorted[0]]
        calc_type += qm_regions_sorted[0]
        if len(qm_regions_sorted) > 1:
             for i in range(1,len(qm_regions_sorted)):
                  combined = base_region.combine(CalculateModel.calculation[qm_regions_sorted[i]])
                  base_region = combined
                  calc_type += qm_regions_sorted[i]
        for calc in ['XTB', 'ChargeField']:
            if calc in CalculateModel.calculation:
                combined = base_region.combine(CalculateModel.calculation[calc])
                calc_type += calc
        CalculateModel.calc_type = calc_type
        CalculationFactory._make_calculation(calc_type)().assign_to_region(region=combined)
        CalculateModel.calculation[calc_type] = combined

    # def _add(self, type, region):
    #     if type != 'QM':
    #         if 'QM' not in CalculateModel.calculation:
    #             raise UserWarning("You must set the high-region QM method first. "
    #                               "See QMzyme.CalculateModel.QM_Method().")
    #     if CalculateModel.calc_type == None:
    #         CalculateModel.calc_type = 'QM'
    #         CalculateModel.calculation[type] = region
    #         return
    #     if type in CalculateModel.calculation:
    #         if type == 'QM':
    #             if CalculateModel.calculation['QM'] == region:
    #                 if CalculateModel.calc_type != 'QM':
    #                     raise UserWarning("Cannot alter QM region after another region "
    #                                       "has been assigned a method. Please reset calculation "
    #                                       "by running QMzyme.CalculateModel.CalculateModel._reset(), then "
    #                                       "you can reassign the methods.")
    #                 del CalculateModel.calculation['QM']
    #             elif 'QM2' in CalculateModel.calculation:
    #                 if CalculateModel.calculation['QM2'] == region:
    #                     del CalculateModel.calculation['QM2']
    #                     CalculateModel.calc_type = CalculateModel.calc_type.remove('QM2')
    #                     type = 'QM2'
    #                     region.method["type"] = 'QM2'
    #                     region.set_atom_segid(region.method["type"])
    #                 else:
    #                     raise UserWarning("Two QM regions have already been defined."
    #                                       "QMzyme currently only supports a maximum "
    #                                       "of two QM methods in a single calculation.")
    #             else:
    #                 type = 'QM2'
    #                 region.method["type"] = 'QM2'
    #                 region.set_atom_segid(region.method["type"])
            
    #         elif CalculateModel.calculation[type] == region:
    #             CalculateModel.calc_type = CalculateModel.calc_type.remove(type)
    #             del CalculateModel.calculation[type]

    #     self.calc_type = type
    #     self.region = region
    #     CalculateModel.calculation[self.calc_type] = region
    
    #     if len(CalculateModel.calculation) > 1:
    #         self.region = region
    #         self._combine_regions()
    #         #self._check_overlap(region)

    # def _combine_regions(self):
    #     #region = self.region
    #     base_region = CalculateModel.calculation['QM']
    #     calc_type = 'QM'
    #     for calc in CalculateModel.calculation:
    #         if calc == 'QM':
    #             continue
    #         calc_type+=calc
    #         #if calc_type in CalculateModel.calculation and calc_type!=self.region.method["type"]:
    #         if calc_type in CalculateModel.calculation and calc_type!=self.calc_type:
    #             # makes 3-scale method possible (QMQMXTB, QMQMChargeField, QMXTBChargeField)
    #             base_region = CalculateModel.calculation[calc_type]
    #         combined_region = base_region+self.region
    #         combined_region.set_method(self.region.method)
    #         combined_region.method["type"] = calc_type
    #         combined_region.method["program"] = base_region.method["program"]
    #         combined_region.name = f'{calc_type}'
    #         #self.calc_type = calc_type

    #     if self.region != combined_region:
    #         print(f"\nCombining regions {self.region} and {base_region} to create {calc_type} region.")
    #         print(f"Region stored in CalculateModel.calculation dictionary under key '{calc_type}'.")
    #         # If self.region does not contain base_region, just add base_region charge and self.region charge
    #         if self.region-base_region == self.region:
    #             combined_region.set_charge(self.region.charge+base_region.charge)
    #         else:
    #             combined_region.guess_charge(verbose=False)
    #             print(f"\nEstimated charge for {combined_region}, the combination of {self.region} and {base_region}, is {combined_region.charge}.")
    #             print(f"\nPlease double check this value. See QMzymeRegion.set_charge() to modify this value if it is incorrect.")
    #     elif self.region == combined_region:
    #         combined_region.set_charge(self.region.charge)
    #     fixed_ids = [atom.id for atom in base_region.get_atoms(attribute='is_fixed', value=True)]
    #     fixed_ids + [atom.id for atom in self.region.get_atoms(attribute='is_fixed', value=True)]
    #     combined_region.set_fixed_atoms(fixed_ids)
    #     CalculateModel.calculation[calc_type] = combined_region
    #     CalculateModel.calc_type = calc_type
    #     CalculationFactory._make_calculation(calc_type)().assign_to_region(region=combined_region)

    def _reset():
        CalculateModel.calculation = {}
        CalculateModel.calc_type = None
    
class CalculationBase:
    def assign_to_region(self, region, charge=None, mult=1):
        self._set_constraints(region)
        self.mult = mult
        region.set_method(self.__dict__)
        self._set_charge(region, charge)
        region.set_atom_segid(region.method["type"])
        CalculateModel().add(type=self.type, region=region)
        #region.set_atom_segid(region.method["type"])
        #self._set_charge(CalculateModel.calculation[self.type], charge)

    def _set_charge(self, region, charge):
        if charge is None:
            if not hasattr(region, "charge"):
                region.guess_charge()
            self.charge = region.charge
        else:
            self.charge = charge

    def _set_constraints(self, region):
        self.freeze_atoms = region.get_indices('is_fixed', True) # these indices are 0 indexed and in order of increasing atom id

class QM_Method(CalculationBase):
    """
    Class to prepare a QMzymeRegion for QM treatment.
    
    Required Parameters
    ---------------------
        region: QMzymeRegion to apply method to
        basis_set: str, defines basis set to use for calculation
        functional: str, defines functional to use for calculation

    Optional Parameters 
    ---------------------
        charge: int, charge of the region. If not provided in parameters, charge will be guessed.
        mult: int, multiplicity of the region. default = 1.
        qm_input: str, default = ""
            Keywords to include in the input file route line to 
            declare any details beyond the basis set and functional.
            E.g. "EmpiricalDispersion=GD3BJ opt freq". Not including anything
            here means the calculation will be a single-point energy calculation.
        qm_end: str, default = ""
            Final line(s) in the input file.
            
    """
    def __init__(self, basis_set, functional, qm_input="", qm_end="", program='orca'):
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
        self.qm_input = self.qm_input.strip()
    
class XTB_Method(CalculationBase):
    """
    Class to prepare a QMzymeRegion for xTB treatment.
    """
    def __init__(self):
        self.type = 'XTB'

class ChargeField_Method(CalculationBase):
    """
    Under development. Class to prepare a QMzymeRegion for ChargeField treatment.
    """
    def __init__(self):
        self.type = 'ChargeField'

class MultiscaleCalculationBase(CalculationBase):
    def assign_to_region(self, region, charge=None, mult=1):
        self.region = region
        self._set_constraints(region)
        self.mult = mult
        region.set_method(self.__dict__)
        region.method["program"] = 'orca'
        #self._set_charge(CalculateModel.calculation[self.type], charge)
        self._set_charge(self.region, charge)
        self._set_qm_input()
        # fix segids
        for atom in self.region.atoms:
            if atom.id in CalculateModel.calculation['QM'].ids:
                setattr(atom, "segid", "QM")
            else:
                #setattr(atom, "segid", CalculateModel.calc_type[2:])
                setattr(atom, "segid", CalculateModel.calc_type[2:])

    def _set_qm_input(self):
        # QM/QM2 or QM/XTB. QMChargeField duck types this method.
        qm_input = f"QM/{CalculateModel.calc_type.replace('QM', '', 1)}" 
        if qm_input not in CalculateModel.calculation['QM'].method['qm_input']:
            qm_input += f" {CalculateModel.calculation['QM'].method['qm_input']}"
        qmmm_section = f"%QMMM {self.qm_input}\n"
        qmmm_section += f" QMATOMS {self._qm_atoms()} END\n"
        qmmm_section += f" Charge_Total {self.region.charge} END"
        self.qm_input = f"{qm_input}\n{qmmm_section}"
        self.region.method['qm_input'] = self.qm_input
        self.region.method['type'] = self.type

    def _qm_atoms(self):
        import numpy as np
        qm_atoms = self.region.get_ix_array_from_ids(ids=CalculateModel.calculation['QM'].ids)
        # qm_atoms = []
        # for i, atom in enumerate(self.region.atoms):
        #     if atom.segid == 'QM':
        #         qm_atoms.append(i)
        range = ''
        for i in np.arange(1, np.max(qm_atoms)+1):
            if i in qm_atoms:
                if i-1 not in qm_atoms:
                    range += "{"+str(i)
                    if i+1 not in qm_atoms:
                        range += "} "
                elif i+1 not in qm_atoms:
                    range += f":{i}"+"} "
        return range.strip()

class _QMQM2_Method(MultiscaleCalculationBase):
    def __init__(self):
        self.type = 'QMQM2'
        qm2_region = CalculateModel.calculation['QM2']
        qm2_input = qm2_region.method['qm_input']
        self.qm_input = qm2_input.lower().replace('opt','')
        self.qm_input = qm2_input.lower().replace('freq','')
        self.qm_input = f"QM2CUSTOMMETHOD '{qm2_input}'"

class _QMXTB_Method(MultiscaleCalculationBase):
    def __init__(self):
        self.type = 'QMXTB'
        self.qm_input = ''

class _QMChargeField_Method(MultiscaleCalculationBase):
    def __init__(self):
        self.type = 'QMChargeField'
        self.qm_input = ''

class CalculationFactory:
    """
    Factory Class to Register Concrete Calculation Method Classes.
    """
    calc_methods = {}
    def _register_method(calculation_type, class_name):
        CalculationFactory.calc_methods[calculation_type] = class_name

    def _make_calculation(calculation_type):
        calculation = CalculationFactory.calc_methods.get(calculation_type)
        if not calculation:
            raise UserWarning(f"Calculation method {calculation} not found.")
        return calculation

possible_methods = {
     'QM': QM_Method,
     'XTB': XTB_Method,
     'ChargeField': ChargeField_Method,
     'QMQM2':_QMQM2_Method,
     'QMXTB': _QMXTB_Method,
     'QMChargeField': _QMChargeField_Method
}
for key, val in possible_methods.items():
    CalculationFactory._register_method(key, val)
     