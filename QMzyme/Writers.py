from QMzyme.aqme.qprep import qprep
import os
import numpy as np
from QMzyme.CalculateModel import CalculateModel


class QMWriter:
    def write(region, memory='24GB', nprocs=12):
        if region is None:
            region = CalculateModel.calculation['QM']
        qprep(**qprep_dict(region.method), mem=memory, nprocs=nprocs)
        # clean up file name because AQME QPREP adds _conf ending.
        rename_file(filename=region.method["files"], program=region.method["program"])
        return 
    
class QMQM2Writer:
    """
    !QM/QM2 WB97X-D3 DEF2-TZVP
    %QMMM QM2CUSTOMMETHOD "PBE D3BJ DEF2-SVP"
    QMATOMS {2:3} {6:13} END
    """
    def write(high_layer=None, low_layer=None, program='orca', **kwargs):
        if high_layer is None:
            high_layer = CalculateModel.calculation['QM']
        if low_layer is None:
            low_layer = CalculateModel.calculation['QM2']

        if filename is None:
            low_layer.method["files"]

        qprep(**qprep_dict(high_layer.method))
        # clean up file name because AQME QPREP adds _conf ending.
        filename = rename_file(filename, low_layer.method["program"])
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


class QMMMWriter:
    def write(high_layer=None, low_layer=None):
        """
        This method will need a QM method assignment for the high_layer, and then 
        takes the point charges for the low_layer. 
        One prerequisite for this is that charges will need to be defined at the atom level.
        """
        if high_layer is None:
            high_layer = CalculateModel.calculation['QM']
        if low_layer is None:
            low_layer = CalculateModel.calculation['ChargeField']
        QMWriter.write(high_layer)
        # remove any qm_atoms from low_layer then proceed
        ChargeFieldWriter.write(low_layer)
        # add line '% pointcharges "pointcharges.pc"' to QM input file

class ChargeFieldWriter:
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
    def write(region):
        for atom in region:
            atom.set_point_charge()
        # write pointcharges.pc

class QMXTBWriter:
    """
    !QM/XTB WB97X-D3 DEF2-TZVP
    %QMMM
    QMATOMS {2:3} {6:13} END

    Need to add check for charge and mult. If xTB and QM region have no atom overlap then just add 
    charge and mult. If xTB contains the entire QM region then just just xTB region charge and mult.
    But if there is partial overlap between the two regions, raise warning. 
    """
    def write(high_layer=None, low_layer=None, total_charge=None, program='orca', **kwargs):
        if high_layer is None:
            high_layer = CalculateModel.calculation['QM']
        if low_layer is None:
            low_layer = CalculateModel.calculation['XTB']
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

def qprep_dict(method_dict):
    d = method_dict
    for key in ["type", "basis_set", "functional"]:
        if key in d:
            del d[key]
    return d

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


