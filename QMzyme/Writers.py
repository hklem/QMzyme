"""
Module responsible for converting one or more QMzymeRegion objects 
"""

from QMzyme.aqme.qprep import qprep
import os
import copy
import numpy as np
from QMzyme.CalculateModel import CalculateModel
from QMzyme.utils import check_filename

writers = {
    'QM': lambda Writer: QMWriter.write(Writer),
    'QMQM2': lambda Writer: QMQM2Writer.write(Writer),
    'QMXTB': lambda Writer: QMXTBWriter.write(Writer),
    'QMChargeField': lambda Writer: QMMMWriter.write(Writer),
}

class Writer:
    """
    Base Writer class that configures common calculation input information and then
    sends information to more specific writer classes to deal with remaining information.

    Parameters
    -----------
    :param model: QMzymeModel containing information needed to write calculation input.
    :type model: :class:`~QMzyme.QMzymeModel.QMzymeModel` (required)
    
    :param filename: Name to be given to calculation input file. 
    Does not need to contain file format suffix.
    :type filename: str (required) Example: filename='1oh0_cutoff3'

    :param writer: Tells the class what format of input file to create. Options
    are entries in the writers dict found in Writers.py
    :type writer: str (required), options: 'QM', 'QMQM2', 'QMXTB', 'QMChargeField'

    :param memory: Memory for the QM calculation 
    (i) Gaussian: total memory; (ii) ORCA: memory per processor.
    :type memory: str (optional, default memory='24GB')

    :param nprocs: Number of processors used in the QM calculation.
    :type nprocs: int (optional, default nprocs=12)

    """
    def __init__(self, filename, writer, memory, nprocs):
        if memory == None:
            memory = '24GB'
        if nprocs == None:
            nprocs = 12
        self.filename = filename
        self.memory = memory
        self.nprocs = nprocs
        writers[writer](self)


def print_details(filename, format):
    filename = check_filename(filename, format)
    pth = os.path.join(os.path.abspath(''), 'QCALC')
    print(f"File {os.path.join(pth, filename)} created.")

class QMWriter:
    def write(Writer, region=None):
        if region is None:
            region = CalculateModel.calculation['QM']
        filename = region.write(Writer.filename)
        region.method["starting_file"] = filename
        qprep(**qprep_dict(region.method), mem=Writer.memory, nprocs=Writer.nprocs)
        if region.method["program"] == 'orca':
            format = '.inp'
        if region.method["program"] == 'gaussian':
            format = '.com'
        print_details(filename, format)
        
    
class QMQM2Writer:
    """
    !QM/QM2 WB97X-D3 DEF2-TZVP
    %QMMM QM2CUSTOMMETHOD "PBE D3BJ DEF2-SVP"
    QMATOMS {2:3} {6:13} END

    Note: The charge and mult above the coordinates section is assigned to the high region (QM atoms).
    """
    def write(Writer, high_region=None, low_region=None, total_charge=None):
        if high_region is None:
            high_region = CalculateModel.calculation['QM']
        if low_region is None:
            low_region = CalculateModel.calculation['QM2']

        if high_region.method["program"] != 'orca':
            raise UserWarning("QM/QM2 calculation only supported for ORCA program.")

        combined = high_region.combine(low_region)
        combined.name = high_region.name+'_'+low_region.name
        filename = combined.write(Writer.filename)
        combined.method = high_region.method
        combined.method["starting_file"] = filename
        combined.method["freeze_atoms"] = combined.get_indices("is_fixed", True)

        if 'QM/QM2' not in combined.method['qm_input']:
            combined.method['qm_input'] = 'QM/QM2 '+combined.method['qm_input']

        qm_atoms = '{'+f'0:{high_region.n_atoms}'+'}'
        qmmm_section = f"%QMMM QM2CUSTOMMETHOD '{low_region.method['qm_input']}'\n"
        qmmm_section += f" QMATOMS {qm_atoms} END\n"
        if total_charge is None:
            total_charge = low_region.charge + high_region.charge
        qmmm_section += f" Charge_Total {total_charge} END"
        combined.method['qm_input'] += f'\n{qmmm_section}'

        qprep(**qprep_dict(high_region.method), mem=Writer.memory, nprocs=Writer.nprocs)
        print_details(filename, 'inp')


class QMMMWriter:
    def write(Writer, high_region=None, low_region=None):
        """
        This method will need a QM method assignment for the high_region, and then 
        takes the point charges for the low_region. 
        One prerequisite for this is that charges will need to be defined at the atom level.
        """
        if high_region is None:
            high_region = CalculateModel.calculation['QM']
        if low_region is None:
            low_region = CalculateModel.calculation['ChargeField']
        QMWriter.write(high_region)
        # remove any qm_atoms from low_region then proceed
        ChargeFieldWriter.write(low_region)
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
    def write(Writer, region):
        for atom in region:
            atom.set_point_charge()
        # write pointcharges.pc

class QMXTBWriter:
    """
    !QM/XTB WB97X-D3 DEF2-TZVP
    %QMMM
    QMATOMS {2:3} {6:13} END

    Note: The charge and mult above the coordinates section is assigned to the high region (QM atoms).

    """
    def write(Writer, high_region=None, low_region=None, total_charge=None):
        if high_region is None:
            high_region = CalculateModel.calculation['QM']
        if low_region is None:
            low_region = CalculateModel.calculation['XTB']
        if high_region.method["program"] != 'orca':
            raise UserWarning("QMXTB calculation only supported for ORCA program.")

        combined = high_region.combine(low_region)
        combined.name = high_region.name+'_'+low_region.name
        filename = combined.write(Writer.filename)
        combined.method = high_region.method
        combined.method["starting_file"] = filename
        combined.method["freeze_atoms"] = combined.get_indices("is_fixed", True)

        if 'QM/XTB' not in combined.method['qm_input']:
            combined.method['qm_input'] = 'QM/XTB '+combined.method['qm_input']

        #high_level_atoms = get_atom_range(high_region.ix_array)
        qm_atoms = '{'+f'0:{high_region.n_atoms}'+'}'
        qmmm_section = "%QMMM\n"
        #qmmm_section += f" QMATOMS {high_level_atoms} END\n"
        qmmm_section += f" QMATOMS {qm_atoms} END\n"
        if total_charge is None:
            total_charge = low_region.charge + high_region.charge
        qmmm_section += f" Charge_Total {total_charge} END"

        if qmmm_section not in combined.method['qm_input']:
            combined.method['qm_input'] += f'\n{qmmm_section}'

        qprep(**qprep_dict(combined.method), mem=Writer.memory, nprocs=Writer.nprocs)
        print_details(filename, 'inp')
        
        # high_level_atoms = get_atom_range(high_region.ix_array)
        # qmmm_section = "%QMMM\n"
        # qmmm_section += f" QMATOMS {high_level_atoms} END\n"
        # if total_charge is None:
        #     total_charge = low_region.charge + high_region.charge
        # qmmm_section += f" Charge_Total {total_charge} END\n"

        # with open(filename, "r") as f:
        #     lines = f.readlines()
        # new_lines = []
        # for i, line in enumerate(lines):
        #     if line.startswith("!"):
        #         if "QM/XTB" not in line:
        #             new_lines.append(f"{line[0]} QM/XTB {line[1:]}")
        #         new_lines.append(qmmm_section)
        #     else:
        #         new_lines.append(line)
        # with open(filename, "w+") as f:
        #     f.writelines(new_lines)


def qprep_dict(method_dict):
    d = copy.copy(method_dict)
    d["files"] = method_dict["starting_file"] # qprep has files keyword
    for key in ["type", "basis_set", "functional", "starting_file"]:
        if key in d:
            del d[key]
    return d

# def get_atom_range(atom_indices: list):
#     """
#     Utility function used for ORCA file input.
#     """
#     range = ''
#     for i in np.arange(1, np.max(atom_indices)+1):
#         if i in atom_indices:
#             if i-1 not in atom_indices:
#                 range += "{"+str(i)
#                 if i+1 not in atom_indices:
#                     range += "} "
#             elif i+1 not in atom_indices:
#                 range += f":{i}"+"} "
#     return range


