###############################################################################
# Code written by Heidi Klem while at
# Colorado State University as a graduate student
# in the Paton and McCullagh groups and at the
# National Institute of Standards and Technology
# as an NRC Postdoc (Fed).
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

import numpy as np
import os
from aqme.qprep import qprep
from aqme.qcorr import qcorr
from QMzyme.utils import get_outlines


def qm_only(qm_input, program='gaussian', suffix='', verbose=True):
    if type(qm_input) is str:
        qm_input = get_outlines(qm_input)
    file = qm_input["starting structure"]
    qprep(files=file,
              charge=qm_input["charge"],
              mult=qm_input["multiplicity"],
              freeze=qm_input["frozen atoms"],
              qm_input=qm_input["level of theory"],
              program=program,
              mem=qm_input["memory"],
              nprocs=qm_input["# processors"],
              suffix=suffix)
    
    # clean up file name because AQME QPREP adds _conf ending.
    calc_file = './QCALC/'+file.split('.pdb')[0]+suffix+'.com'
    os.rename('./QCALC/'+file.split('.pdb')[0]+suffix+'_conf_1.com', calc_file)
    verbose_str = ''
    if verbose is True:
        verbose_str += "INITIALIZING... AQME.QPREP QM INPUT FILE GENERATION\n"
        verbose_str += "STARTING_STRUCTURE: {}\n".format(file)
        verbose_str += "CALCULATION: {}\n".format(qm_input["level of theory"])
        verbose_str += "CHARGE: {}\n".format(qm_input["charge"])
        verbose_str += "MULTIPLICITY: {}\n".format(qm_input["multiplicity"])
        verbose_str += "MEM: {}\n".format(qm_input["memory"])
        verbose_str += "NPROCS: {}\n".format(qm_input["# processors"])
        verbose_str += "FROZEN_ATOMS: {}\n".format(qm_input["frozen atoms"])
        verbose_str += "CALCULATION_FILE: {}\n".format(calc_file)

    # write json
    info = {
            'Calculation file': calc_file,
            'Calculation': qm_input["level of theory"],
            'Charge': chrg,
            'Multiplicity': mult,
            'Memory': mem,
            'Number of processors': nprocs,
            }
    self.to_dict(section='QM preparation', dict=info)
    
    return verbose_str