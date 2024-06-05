"""
Module responsible for converting one or more QMzymeRegion objects 
"""

from QMzyme.aqme.qprep import qprep
import os
import numpy as np
import copy
from QMzyme.CalculateModel import CalculateModel
from QMzyme.utils import check_filename
import abc


### Writer abstract class ###
class Writer(abc.ABC):
    @abc.abstractmethod
    def __init__(self, filename, memory, nprocs, full_region):
        self.memory = memory
        self.nprocs = nprocs
        self.full_region = full_region
        filename = full_region.write(filename)
        self.filename = filename
        self.full_region.method["starting_file"] = filename
        self.set_constraints()

    @abc.abstractmethod
    def write(self):
        ...

    def set_constraints(self):
        self.full_region.method['freeze_atoms'] = self.full_region.get_indices('is_fixed', True) 
        # These indices are 0 indexed and in order of increasing atom id

### Concrete writer subclasses ###
class QMWriter(Writer):
    """
    Writes a QM input file for ORCA or Gaussian using AQME qprep. 

    :param filename: Name to be given to calculation input file. 
        Does not need to contain file format suffix.
    :type filename: str (required) Example: filename='1oh0_cutoff3'

    :param memory: Memory for the QM calculation 
        (i) Gaussian: total memory; (ii) ORCA: memory per processor.
    :type memory: str (optional, default memory='24GB')

    :param nprocs: Number of processors used in the QM calculation.
    :type nprocs: int (optional, default nprocs=12)

    :returns:
        Saves QM input file and region pdb file to working directory. 
    """
    def __init__(self, filename, memory, nprocs, full_region=None):
        if full_region is None:
            full_region = CalculateModel.calculation['QM']
        super().__init__(filename, memory, nprocs, full_region)
        self.write()

    def write(self):
        qprep(**qprep_dict(self.full_region.method), mem=self.memory, nprocs=self.nprocs)
        if self.full_region.method["program"] == 'orca':
            format = '.inp'
        if self.full_region.method["program"] == 'gaussian':
            format = '.com'
        print_details(self.filename, format)

   
class QMQM2Writer(Writer):
    """
    Writes a QM input file for ORCA or Gaussian using AQME qprep. 

    :param filename: Name to be given to calculation input file. 
        Does not need to contain file format suffix.
    :type filename: str (required) Example: filename='1oh0_cutoff3'

    :param memory: Memory for the QM calculation 
        (i) Gaussian: total memory; (ii) ORCA: memory per processor.
    :type memory: str (optional, default memory='12GB')

    :param nprocs: Number of processors used in the QM calculation.
    :type nprocs: int (optional, default nprocs=12)

    :param high_region: QMzymeRegion with assigned QM method. If not provided, the
    code will search in `CalculateModel.calculation` for the 'QM' entry.
    :type high_region: :class:`~QMzyme.QMzymeRegion.QMzymeRegion`, optional

    :param low_region: QMzymeRegion with assigned QM method. If not provided, the
    code will search in `CalculateModel.calculation` for the 'QM2' entry.
    :type low_region: :class:`~QMzyme.QMzymeRegion.QMzymeRegion`, optional

    :returns:
        Saves QM input file and combined regions pdb file to working directory. 

    :notes:

        Example text added to QM input file:
        .. code-block:: bash

            !QM/QM2 WB97X-D3 DEF2-TZVP
            %QMMM QM2CUSTOMMETHOD "PBE D3BJ DEF2-SVP"
            QMATOMS {2:3} {6:13} END

        The charge and mult above the coordinates section is assigned to the high region (QM atoms).
    """
    def __init__(self, filename, memory, nprocs, full_region=None, total_charge=None, total_mult=None, high_region=None, low_region=None):
        if high_region is None:
            high_region = CalculateModel.calculation['QM']
        if high_region.method["program"] != 'orca':
            raise UserWarning("QM/QM2 calculation only supported for ORCA program.")
        if low_region is None:
            low_region = CalculateModel.calculation['QM2']
        if full_region is None:
            full_region = CalculateModel.calculation['QMQM2']
        super().__init__(filename, memory, nprocs, full_region)
        if total_charge is None:
            total_charge = high_region.method["charge"]
        if total_mult is None:
            total_mult = high_region.method["mult"]
        self.total_charge = total_charge
        self.total_mult = total_mult
        self.write(full_region, high_region)

    def write(self, full_region, high_region):
        method = qprep_dict(full_region.method)
        # charge and mult here are just for the high region qm_atoms
        method['charge'] = self.total_charge
        method['mult'] = self.total_mult

        qprep(**method, mem=self.memory, nprocs=self.nprocs)
        print_details(self.filename, 'inp')


class QMMMWriter(Writer):
    """
    Under development.
    """
    pass
#     def __init__(self, filename, memory, nprocs, high_region=None, low_region=None, total_charge=None):
#         super().__init__(filename, memory, nprocs)
#         self.write(high_region, low_region, total_charge)

#     def write(self, high_region=None, low_region=None, total_charge=None):
#         """
#         This method will need a QM method assignment for the high_region, and then 
#         takes the point charges for the low_region. 
#         One prerequisite for this is that charges will need to be defined at the atom level.
#         """
#         if high_region is None:
#             high_region = CalculateModel.calculation['QM']
#         if low_region is None:
#             low_region = CalculateModel.calculation['ChargeField']

#         QMWriter(self.filename, self.memory, self.nprocs).write(high_region)
#         # remove any qm_atoms from low_region then proceed
#         ChargeFieldWriter.write(low_region)
#         # add line '% pointcharges "pointcharges.pc"' to QM input file

class ChargeFieldWriter:
    """
    Under development.

    .. code-block:: bash

        3 # number of point charges
        -0.834  -1.3130   0.0000  -0.0310 # charge in a.u. x y z in Ang
        0.417  -1.8700   0.7570   0.1651
        0.417  -1.8700  -0.7570   0.1651

    **From Orca manual**: 
    However, it should be noted that ORCA treats point charges from an external 
    file differently than “Q” atoms. When using an external point charge file, 
    the interaction between the point charges is not included in the nuclear energy. 
    This behavior originates from QM/MM, where the interactions among the point charges 
    is done by the MM program. These programs typically use an external point charge 
    file when generating the ORCA input. To add the interaction of the point charges 
    to the nuclear energy, the DoEQ keyword is used either in the simple input or the 
    %method block as shown below.

    # A non QM/MM pointcharge calculation

    .. code-block:: bash

        ! DoEQ
    
        %pointcharges "pointcharges.pc"

        %method
            DoEQ true
        end
    """
    pass
#     def write(Writer, region):
#         for atom in region:
#             atom.set_point_charge()
#         # write pointcharges.pc

class QMXTBWriter(Writer):
    """
    Writes a QM input file for ORCA or Gaussian using AQME qprep. 

    :param filename: Name to be given to calculation input file. 
        Does not need to contain file format suffix.
    :type filename: str (required) Example: filename='1oh0_cutoff3'

    :param memory: Memory for the QM calculation 
        (i) Gaussian: total memory; (ii) ORCA: memory per processor.
    :type memory: str (optional, default memory='12GB')

    :param nprocs: Number of processors used in the QM calculation.
    :type nprocs: int (optional, default nprocs=12)

    :param high_region: QMzymeRegion with assigned QM method. If not provided, the
    code will search in `CalculateModel.calculation` for the 'QM' entry.
    :type high_region: :class:`~QMzyme.QMzymeRegion.QMzymeRegion`, optional

    :param low_region: QMzymeRegion with assigned xTB method. If not provided, the
    code will search in `CalculateModel.calculation` for the 'XTB' entry.
    :type low_region: :class:`~QMzyme.QMzymeRegion.QMzymeRegion`, optional

    :returns:
        Saves QM input file and combined regions pdb file to working directory. 

    :notes:
        Example input:

        .. code-block:: bash

            !QM/XTB WB97X-D3 DEF2-TZVP
            %QMMM
            QMATOMS {2:3} {6:13} END
 
        The charge and mult above the coordinates section is assigned 
        to the high region (QM atoms).
    """
    def __init__(self, filename, memory, nprocs, full_region=None, high_region=None, low_region=None, total_charge=None, total_mult=None):
        if high_region is None:
            high_region = CalculateModel.calculation['QM']
        if high_region.method["program"] != 'orca':
            raise UserWarning("QM/QM2 calculation only supported for ORCA program.")
        if low_region is None:
            low_region = CalculateModel.calculation['XTB']
        if full_region is None:
            full_region = CalculateModel.calculation['QMXTB']
        super().__init__(filename, memory, nprocs, full_region)
        if total_charge is None:
            total_charge = high_region.method["charge"]
        if total_mult is None:
            total_mult = high_region.method["mult"]
        self.total_charge = total_charge
        self.total_mult = total_mult
        self.write(full_region, high_region)

    def write(self, full_region, high_region):
        method = qprep_dict(full_region.method)
        # charge and mult here are just for the high region qm_atoms
        method['charge'] = self.total_charge
        method['mult'] = self.total_mult
        qprep(**method, mem=self.memory, nprocs=self.nprocs)
        print_details(self.filename, 'inp')
        

### Factory Class to Register Concrete Writer Classes ###
class WriterFactory:
    """
    Factory Class to Register Concrete Writer Classes.
    """
    writers = {}

    @staticmethod
    def register_writer(writer_type, writer):
        WriterFactory.writers[writer_type] = writer

    @staticmethod
    def make_writer(writer_type, filename, memory, nprocs, **kwargs):
        """
        Instantiates and returns the specific concrete writer subclass.
        """
        writer = WriterFactory.writers.get(writer_type)
        if not writer:
            raise UserWarning(f"No writer detected for calculation type {writer_type}.")
        return writer(filename, memory, nprocs, **kwargs)
    
WriterFactory.register_writer('QM', QMWriter)
WriterFactory.register_writer('QMQM2', QMQM2Writer)
WriterFactory.register_writer('QMXTB', QMXTBWriter)
WriterFactory.register_writer('QMChargeField', QMMMWriter)


### Main Writer Class that GenerateModel Calls ###

# class Writer:
#     """
#     Base Writer class that configures common calculation input information and then
#     sends information to more specific writer classes to deal with remaining information.

#     Parameters
#     -----------
#     :param filename: Name to be given to calculation input file. 
#         Does not need to contain file format suffix.
#     :type filename: str (required) Example: filename='1oh0_cutoff3'

#     :param writer: Tells the class what format of input file to create. Options
#         are entries in the writers dict found in Writers.py
#     :type writer: str (required), available options are in WritersRegistry.writers.keys()

#     :param memory: Memory for the QM calculation 
#         (i) Gaussian: total memory; (ii) ORCA: memory per processor.
#     :type memory: str (optional, default memory='12GB')

#     :param nprocs: Number of processors used in the QM calculation.
#     :type nprocs: int (optional, default nprocs=12)

#     """
#     def __init__(self, filename, memory, nprocs, writer=None):
#         self.filename = filename
#         self.memory = memory
#         self.nprocs = nprocs
#         if writer is None:
#             writer = WritersRegistry._get_writer(CalculateModel.calc_type)
#         else:
#             writer = WritersRegistry._get_writer(writer)
#         self.writer = writer

#     def write(self):
#         self.writer(self.filename, self.memory, self.nprocs)


### Auxilliary functions ###
def print_details(filename, format):
    filename = check_filename(filename, format)
    pth = os.path.join(os.path.abspath(''), 'QCALC')
    print(f"File {os.path.join(pth, filename)} created.")

def qprep_dict(method_dict):
    d = copy.copy(method_dict)
    d["files"] = method_dict["starting_file"] # qprep has files keyword
    for key in ["type", "basis_set", "functional", "starting_file", "region"]:
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
#     return range.strip()
