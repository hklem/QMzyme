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


qm_params = {
        "program": None,
        "qm_input": None,
        "qm_end": None,
        "charge": None,
        "mult": None,
        "chk": None,
        "chk_path": None,
        "mem": None,
        "nprocs": None,
        "gen_atoms": None,
        "bs_gen": None,
        "bs_nogen": None,
        "freeze": None
    }


def write_QM(layer):
    qprep(**layer.method)
    # clean up file name because AQME QPREP adds _conf ending.
    rename_file(layer.method["files"], layer.method["program"])

def write_QMxTB(filename, high_layer, low_layer, total_charge, program='orca', **kwargs):
    """
    !QM/XTB WB97X-D3 DEF2-TZVP
    %QMMM
      QMATOMS {2:3} {6:13} END

    Need to add check for charge and mult. If xTB and QM region have no atom overlap then just add 
    charge and mult. If xTB contains the entire QM region then just just xTB region charge and mult.
    But if there is partial overlap between the two regions, raise warning. 
    """
    qprep(**kwargs)
    # clean up file name because AQME QPREP adds _conf ending.
    filename = rename_file(filename, kwargs["program"])
    qmmm_section = ""
    with open(filename, 'r') as f:
        lines = f.readlines()

    high_level_atoms = get_atom_range(high_level_atoms)
    qmmm_section = "%QMMM\n"
    qmmm_section += f" QMATOMS {high_level_atoms} END\n"
    qmmm_section += f" Charge_Total {total_charge} END\n"

    new_lines = []
    for i, line in enumerate(lines):
        if line.startswith("!"):
            if "QM/XTB" not in line:
                new_lines.append(f"{line[0]} QM/XTB {line[1:]}")
            new_lines.append(qmmm_section)
        else:
            new_lines.append(line)
    with open(filename, "w+") as f:
        f.writelines(new_lines)

def write_QMQM(filename, high_layer, low_layer, program='orca', **kwargs):
    """
    !QM/QM2 WB97X-D3 DEF2-TZVP
    %QMMM QM2CUSTOMMETHOD "PBE D3BJ DEF2-SVP"
      QMATOMS {2:3} {6:13} END
    """
    qprep(**kwargs)
    # clean up file name because AQME QPREP adds _conf ending.
    filename = rename_file(filename, kwargs["program"])
    qmqm2_section = ""
    with open(filename, 'r') as f:
        lines = f.readlines()

    high_level_atoms = get_atom_range(high_level_atoms)
    qmmm_section = f"%QMMM QM2CUSTOMMETHOD '{low_layer['qm_input']}'\n"
    qmmm_section += f" QMATOMS {high_level_atoms} END\n"

    new_lines = []
    for i, line in enumerate(lines):
        if line.startswith("!"):
            if "QM/QM2" not in line:
                new_lines.append(f"{line[0]} QM/QM2 {line[1:]}")
            new_lines.append(qmmm_section)
        else:
            new_lines.append(line)
    with open(filename, "w+") as f:
        f.writelines(new_lines)

def write_QMMM(high_layer, low_layer):
    """
    This method will need a QM method assignment for the high_layer, and then 
    takes the point charges for the low_layer. 
    One prerequisite for this is that charges will need to be defined at the atom level.
    """
    write_QM(high_layer, high_layer.method)
    # remove any qm_atoms from low_layer then proceed
    write_charge_field()
    # add line '% pointcharges "pointcharges.pc"' to QM input file

def write_charge_field(region):
    """
    3 # number of point charges
     -0.834  -1.3130   0.0000  -0.0310 # charge in a.u. x y z in Ang
      0.417  -1.8700   0.7570   0.1651
      0.417  -1.8700  -0.7570   0.1651

    From Orca manual: 
    However, it should be noted that ORCA treats point charges from an external 
    file differently than “Q” atoms. When using an external point charge file, 
    the interaction between the point charges is not included in the nuclear energy. 
    This behavior originates from QM/MM, where the interactions among the point charges 
    is done by the MM program. These programs typically use an external point charge 
    file when generating the ORCA input. To add the interaction of the point charges 
    to the nuclear energy, the DoEQ keyword is used either in the simple input or the 
    %method block as shown below.

    # A non QM/MM pointcharge calculation
    ! DoEQ
    
    %pointcharges "pointcharges.pc"

    %method
        DoEQ true
    end
    """
    for atom in region:
        atom.set_point_charge()
    # write pointcharges.pc


def get_atom_range(atom_indices: list):
    """
    Utility function used for ORCA file input.
    """
    range = ''
    for i in np.arange(1, np.max(atom_indices)+1):
        if i in atom_indices:
            if i-1 not in atom_indices:
                range += "{"+str(i)
                if i+1 not in atom_indices:
                    range += "} "
            elif i+1 not in atom_indices:
                range += f":{i}"+"} "
    return range


def rename_file(filename, program):
    """
    Function that cleans up aqme file naming convention.
    """
    if program == 'orca':
        end = 'inp'
    if program == 'gaussian':
        end = 'com'
    calc_file = f"./QCALC/{filename.split('.pdb')[0]}.{end}"
    try:
        os.rename(f"./QCALC/{filename.split('.pdb')[0]}_conf_1.{end}", calc_file)
        with open(calc_file, "r") as f:
            lines = f.readlines()
        for i,line in enumerate(lines):
            if "_conf_1" in line:
                print(lines[i])
                lines[i] = line.replace("_conf_1", "")
                print(lines[i])
        with open(calc_file, "w+") as f:
            f.writelines(lines)
    except:
        pass
    return calc_file

        
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

class QM:
    """
    Class to prepare a QMzymeRegion for QM treatment.
    
    Required Parameters
    ====================
    region: QMzymeRegion to apply method to
    basis_set: str, defines basis set to use for calculation
    functional: str, defines functional to use for calculation
    charge: int, charge of the region
    mult: int, multiplicity of the region

    Optional Parameters 
    ====================
    qm_input: str, default = ""
        Keywords to include in the input file route line to 
        declare any details beyond the basis set and functional.
        E.g. "EmpiricalDispersion=GD3BJ opt freq". Not including anything
        here means the calculation will be a single-point energy calculation.
    qm_end: str, default = ""
        Final line(s) in the input file

    """
    def __init__(self, region, basis_set, functional, charge, mult, qm_input="", qm_end="", program='orca'):
        """
        :param region: QMzyme region to treat at the QM level.
        :type region: QMzymeRegion
        """
        self.qm_input = qm_input
        self.basis_set = basis_set
        self.functional = functional
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
        region.set_method(self.__dict__, type="QM")

    def _set_constraints(self, region):
        self.freeze_atoms = region.get_indices('is_fixed', True)

    def _set_qm_input(self):
        for info in [self.functional, self.basis_set]:
            if info not in self.qm_input:
                self.qm_input = f"{info} {self.qm_input}"
    
