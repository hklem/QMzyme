###############################################################################
# Code written by Heidi Klem while at
# Colorado State University as a graduate student
# in the Paton and McCullagh groups and at the
# National Institute of Standards and Technology
# as an NRC Postdoc (Fed).
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

'''
Module in charge of create input files for QM-only or QM/MM calculations. This 
module integrates the `AQME QPREP <https://aqme.readthedocs.io/en/latest/API/aqme.qprep.html>`_ 
workflow.

Required input
...............

    *   QMzyme model object.
    *   Model charge and multiplicity.

Notes
...............

    *   Currently optimized to generate QM-only Gaussian input files.
'''

import os
from QMzyme import BiopythonWrapper
from QMzyme import aqme
from QMzyme.aqme.qprep import qprep


class CalculateModel():

    def __init__(self, model, charge = None, mult = None):
        if charge is None:
            if hasattr(model, 'charge'):
                charge = model.charge
            else:
                print("Charge must be specified by charg=(int).")
        if mult is None:
            if hasattr(model, 'multiplicity'):
                mult = model.multiplicity
            else:
                print("Multiplicity must be specified by 'mult=(int)'.")
        self.charge = model.charge
        self.mult = model.multiplicity
        if not hasattr(model, 'calculations'):
            setattr(model, 'calculations', [])
        model.calculations.append(self)
        if self.suffix == '':
            self.suffix = f'calc_{len(model.calculations)}'
        #self.coords = [a.coord for a in model.list_atoms()]
        #self.elements = [a.element for a in model.list_atoms()]
        #self.atom_ids = [a.serial_number for a in model.list_atoms()]


    # def QM_only(self, model, functional, basis_set, frozen_atoms = None):
    #     '''
    #     Requirements:
    #         - elements
    #         - coords
    #         - functional
    #         - basis set(s)
    #         - constraints
    #         - 
    #     '''
    #     self.type = 'QM_only'
    #     self.functional = functional
    #     self.basis_set = basis_set
    #     self.frozen_atoms = frozen_atoms
    #     self.charge = model.charge
    #     self.mult = model.multiplicity
    #     model.calculations.append(self)

class CalculateQM(CalculateModel):
    '''
    Requirements:
            - elements
            - coords
            - functional
            - basis set
            - charge
            - mult
    '''
    def __init__(self, model, functional, 
                basis_set, opt=True, freq=True, freeze_atoms = [], 
                charge = 0, mult = 1, mem='32GB', 
                nprocs=16, program='gaussian', suffix=''):
    
        self.type = 'QM_only'
        self.functional = functional
        self.basis_set = basis_set
        self.frozen_atoms = BiopythonWrapper.get_atom_idx(model, freeze_atoms)
        self.status = None 
        self.mem = mem
        self.nprocs = nprocs
        self.program = program
        self.pdb_file = BiopythonWrapper.write_pdb(model)
        self.opt = opt
        self.freq = freq
        self.charge = charge
        self.mult = mult
        self.suffix = suffix
        CalculateModel.__init__(self, model, charge, mult)
        self.aqme_qprep()


    def aqme_qprep(self):
        if self.suffix.startswith('_'):
            self.suffix = self.suffix.split('_')[-1]
        qm_input = ''
        if self.opt is True:
            qm_input += 'opt '
        if self.freq is True:
            qm_input += 'freq=noraman '
        qm_input += f'{self.functional} {self.basis_set}'
        qprep(files=self.pdb_file,
              charge=self.charge,
              mult=self.mult,
              freeze=self.frozen_atoms,
              qm_input=qm_input,
              program=self.program,
              mem=self.mem,
              nprocs=self.nprocs,
              suffix=self.suffix)
        # clean up file name because AQME QPREP adds _conf ending.
        calc_file = f"./QCALC/{self.pdb_file.split('.pdb')[0]}_{self.suffix}.com"
        try:
            os.rename(f"./QCALC/{self.pdb_file.split('.pdb')[0]}_conf_1_{self.suffix}.com", calc_file)
        except:
            pass



    #def QM_MM(self):
    #def QM_xTB(self):
