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
    def _add(self, type, region):
        if type != 'QM':
            if 'QM' not in CalculateModel.calculation:
                raise UserWarning("You must set the high-region QM method first. "
                                  "See QMzyme.CalculateModel.QM_Method().")
        if CalculateModel.calc_type == None:
            CalculateModel.calc_type = 'QM'
            CalculateModel.calculation[type] = region
            return
        if type in CalculateModel.calculation:
            if type == 'QM':
                if CalculateModel.calculation['QM'] == region:
                    if CalculateModel.calc_type != 'QM':
                        raise UserWarning("Cannot alter QM region after another region "
                                          "has been assigned a method. Please reset calculation "
                                          "by running QMzyme.CalculateModel.CalculateModel._reset(), then "
                                          "you can reassign the methods.")
                    del CalculateModel.calculation['QM']
                elif 'QM2' in CalculateModel.calculation:
                    if CalculateModel.calculation['QM2'] == region:
                        del CalculateModel.calculation['QM2']
                        CalculateModel.calc_type = CalculateModel.calc_type.remove('QM2')
                        type = 'QM2'
                        region.method["type"] = 'QM2'
                        region.set_atom_segid(region.method["type"])
                    else:
                        raise UserWarning("Two QM regions have already been defined."
                                          "QMzyme currently only supports a maximum "
                                          "of two QM methods in a single calculation.")
                else:
                    type = 'QM2'
                    region.method["type"] = 'QM2'
                    region.set_atom_segid(region.method["type"])
            
            elif CalculateModel.calculation[type] == region:
                CalculateModel.calc_type = CalculateModel.calc_type.remove(type)
                del CalculateModel.calculation[type]

        # elif calc == 'QM':
        #     if 'QM2' in CalculateModel.calculation:
        #         if CalculateModel.calculation['QM2'] == region:
        #             del CalculateModel.calculation['QM2']
        #         else:
        #             raise UserWarning("Two QM regions have already been defined."
        #                               "QMzyme currently only supports a maximum "
        #                               "of two QM regions in a single calculation.")
        #     calc = 'QM2'
        #     region.method["type"] = 'QM2'
        #     region.set_atom_segid(region.method["type"])

            # if calc in CalculateModel.calculation:
            #     if 'QM2' in CalculateModel.calculation:
            #         raise UserWarning("Two QM regions have already been defined."
            #                           "QMzyme currently only supports a maximum "
            #                           "of two QM regions in a single calculation.")
            #     calc = 'QM2'
            #     region.method["type"] = 'QM2'
            #     region.set_atom_segid(region.method["type"])
            #     #print(f"WARNING: A QM region already defined. Setting QMzymeRegion {region.name} as QM2.")

        self.calc_type = type
        self.region = region
        CalculateModel.calculation[self.calc_type] = region
    
        if len(CalculateModel.calculation) > 1:
            self.region = region
            self._combine_regions()
            #self._check_overlap(region)

    def _combine_regions(self):
        #region = self.region
        base_region = CalculateModel.calculation['QM']
        calc_type = 'QM'
        for calc in CalculateModel.calculation:
            if calc == 'QM':
                continue
            calc_type+=calc
            #if calc_type in CalculateModel.calculation and calc_type!=self.region.method["type"]:
            if calc_type in CalculateModel.calculation and calc_type!=self.calc_type:
                # makes 3-scale method possible (QMQMXTB, QMQMChargeField, QMXTBChargeField)
                base_region = CalculateModel.calculation[calc_type]
            combined_region = base_region+self.region
            combined_region.set_method(self.region.method)
            combined_region.method["type"] = calc_type
            combined_region.method["program"] = base_region.method["program"]
            combined_region.name = f'{calc_type}'
            #self.calc_type = calc_type

        if self.region != combined_region:
            print(f"\nCombining regions {self.region} and {base_region} to create {calc_type} region.")
            print(f"Region stored in CalculateModel.calculation dictionary under key '{calc_type}'.")
            # If self.region does not contain base_region, just add base_region charge and self.region charge
            if self.region-base_region == self.region:
                combined_region.set_charge(self.region.charge+base_region.charge)
            else:
                combined_region.guess_charge(verbose=False)
                print(f"\nEstimated charge for {combined_region}, the combination of {self.region} and {base_region}, is {combined_region.charge}.")
                print(f"\nPlease double check this value. See QMzymeRegion.set_charge() to modify this value if it is incorrect.")
        elif self.region == combined_region:
            combined_region.set_charge(self.region.charge)
        fixed_ids = [atom.id for atom in base_region.get_atoms(attribute='is_fixed', value=True)]
        fixed_ids + [atom.id for atom in self.region.get_atoms(attribute='is_fixed', value=True)]
        combined_region.set_fixed_atoms(fixed_ids)
        CalculateModel.calculation[calc_type] = combined_region
        CalculateModel.calc_type = calc_type
        CalcMethodsRegistry._get_calc_method(calc_type)().assign_to_region(region=combined_region)

    def _check_overlap(self, region):
        high_region = CalculateModel.calculation['QM']
        low_region = region

        high_region.segid = 'QM'
        common_atoms = high_region.get_overlap(low_region)
        if len(common_atoms) != 0:
            #residues = [atom.resname+str(atom.resid) for atom in common_atoms]
            residues = [atom.resid for atom in common_atoms]
            residues = [high_region.get_residue(resid) for resid in list(set(residues))]
            subtracted = low_region.subtract(high_region)
            subtracted.guess_charge(verbose=False)
            subtracted.method = low_region.method
            subtracted.method["charge"] = subtracted.charge
            subtracted.method["freeze_atoms"] = subtracted.get_indices(attribute='is_fixed', value=True)
            subtracted.name = low_region.name
            print(f"\nWARNING: Region overlap detected. The following residue(s) were found in both regions: {list(set(residues))}.")
            print(f"Removing duplicate atoms and recalculating charge for calculation.")
            print(f"Subtracted region has a charge of {subtracted.charge}. Both regions combine for a total charge of {subtracted.charge+high_region.charge}.")
            print(f"The original region {low_region.name} remains unchanged in the GenerateModel object 'regions' attribute.")
                        
            CalculateModel.calculation[low_region.method["type"]] = subtracted

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
        CalculateModel()._add(type=self.type, region=region)
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
        #self._set_charge(CalculateModel.calculation[self.type], charge)
        #self._set_charge(self.region, charge)
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
        qm_input = f'QM/{self.type[-3:]}' 
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


class CalcMethodsRegistry:
    """
    Factory Class to Register Concrete Calculation Method Classes.
    """
    calc_methods = {}
    def _register_method(calc_type, class_name):
        CalcMethodsRegistry.calc_methods[calc_type] = class_name

    def _get_calc_method(calc_type):
        calc_method = CalcMethodsRegistry.calc_methods.get(calc_type)
        if not calc_method:
            raise UserWarning(f"Calculation method {calc_method} not found.")
        return calc_method
    
CalcMethodsRegistry._register_method('QM', QM_Method)
CalcMethodsRegistry._register_method('XTB', XTB_Method)
# CalcMethodsRegistry._register_method('ChargeField', ChargeField_Method)
CalcMethodsRegistry._register_method('QMQM2', _QMQM2_Method)
CalcMethodsRegistry._register_method('QMXTB', _QMXTB_Method)
# CalcMethodsRegistry._register_method('QMChargeField', _QMChargeField_Method)