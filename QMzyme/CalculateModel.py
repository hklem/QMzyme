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

    This class also determines if there are any model issues, such as guessing region charge if partial 
    overlap between occurs between regions (you do not want to double-count charge from duplicate atoms).

    This class will also alert you if you if the methods you are assigning to regions do not add up to a 
    multi-scale calculation method that exists as a subclass of :class:`~QMzyme.Writers.Writer`, because 
    QMzyme would not know how to write that input file without the corresponding writer class.
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
        """
        This method is triggered in :class:`~QMzyme.GenerateModel.GenerateModel`: when 
        :func:`~QMzyme.GenerateModel.GenerateModel.write_input` is called.

        It ensures the overall multiscale method contains all atoms of the QMzymeModel regions
        that have been assigned a method. 
        """
        
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

    def _reset():
        CalculateModel.calculation = {}
        CalculateModel.calc_type = None
    
class CalculationBase:
    """
    Base class for all single-method classes. Contains the inherited method assign_to_region().
    """
    def assign_to_region(self, region, charge=None, mult=1):
        """
        Connects a calculation method instance to a QMzymeRegion. This method also searches for any
        atoms in the region with attribute ``is_fixed=True`` to store what atoms will be constrained
        in the calculation file.

        :param region: Region to assign method istance to.
        :type region: :class:`~QMzyme.QMzymeRegion.QMzymeRegion`, required

        :param charge: Charge of the region. If not specified and the QMzymeRegion instance does not have
            a `charge` attribute, the QMzymeRegion :func:`~QMzyme.QMzymeRegion.QMzymeRegion.guess_charge` method 
            will be called.
        :type charge: int, default=None
        """
        self._set_constraints(region)
        self.mult = mult
        region.set_method(self.__dict__)
        self._set_charge(region, charge)
        region.set_atom_segid(region.method["type"])
        CalculateModel().add(type=self.type, region=region)

    def _set_charge(self, region, charge):
        if charge is None:
            if not hasattr(region, "charge"):
                region.guess_charge(verbose=False)
            self.charge = region.charge
        else:
            self.charge = charge

    def _set_constraints(self, region):
        self.freeze_atoms = region.get_indices('is_fixed', True) # these indices are 0 indexed and in order of increasing atom id

class QM_Method(CalculationBase):
    """
    Class to set a quantum mechanics method for a QMzymeRegion.
    
    :param region: QMzymeRegion to apply method to
    :type region: :class:`~QMzyme.QMzymeRegion.QMzymeRegion`, required
    
    :param basis_set: Defines the basis set to use for calculation.
    :type param: str, required
    
    :param functional: Defines the functional to use for calculation.
    :type functional: str, required
    
    :param charge: Charge of the region. If not provided in parameters, charge will be guessed.
    :type charge: int

    :param mult: Multiplicity of the region.
    :type mult: int, default=1

    :param qm_input: Keywords to include in the input file route line to declare any details 
        beyond the basis set and functional. E.g. "EmpiricalDispersion=GD3BJ opt freq". Not including 
        anything here means the calculation will be a single-point energy calculation with no frequency analysis.
    :type qm_input: str, default=''
    
    :param qm_end: Final line(s) in the input file.
    :type qm_end: str, default=''
            
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
    Class to set an xTB method for a QMzymeRegion.
    """
    def __init__(self):
        self.type = 'XTB'

class ChargeField_Method(CalculationBase):
    """
    Under development. Class to prepare a QMzymeRegion for ChargeField treatment (aka charge embedding).
    """
    def __init__(self):
        self.type = 'ChargeField'

class MultiscaleCalculationBase(CalculationBase):
    """
    Base class for all multi-scale methods.
    """
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
     
