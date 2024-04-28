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
import QMzyme
from QMzyme import MDAnalysisWrapper
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
        self.frozen_atoms = MDAnalysisWrapper.get_atom_idx(model, freeze_atoms)
        self.status = None 
        self.mem = mem
        self.nprocs = nprocs
        self.program = program
        self.pdb_file = MDAnalysisWrapper.write_pdb(model, 'truncated_model.pdb')
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

    def amber_gau(basis='6-31G*',
                  method='BLYP',
                  scf_conv=8,
                  num_threads=1,
                  executable=None,
                  use_template=0,
                  ntpr=None,
                  dipole=0,
                  mem='256MB'
                  ):
        """
        Currently under development.

        Parameters as described in the Amber23 Manual. Note the default 
        values match the default Amber settings and are not to be understood
        as best practice recommendations.

        :param basis: Basis set type to be used in the calculation. Any 
        basis set that is natively supported by Gaussian can be used. 
        Examples are the single zeta, split valence or triple zeta Pople 
        type basis sets STO-3G, 3- 21G, 6-31G and 6-311G. The split-valence 
        or triple zeta basis sets can be augmented with diffuse functions 
        on heavy atoms or additionally hydrogen by adding one or two plus 
        signs, respectively, as in 6-31++G. Polarization functions on heavy 
        atoms or additionally hydrogens are used by adding one or two stars, 
        respectively, as in 6-31G**. 
        :type basis: str, default: '6-31G*'

        :param method: Method to be used in the calculation. Can either be 
        one of the WFT models for which Gaussian supports gradients, for 
        example RHF or MP2, or some supported DFT functional. Popular choices 
        are BLYP, PBE and B3LYP. 
        :type method: str, default: 'BLYP'

        :param scf_conv: Threshold upon which to stop the SCF procedure. 
        The tested error is the commutator of the Fock matrix and the density 
        matrix. Convergence is considered to be achieved if the maximum 
        element of the commutator (which is zero for an optimized wave 
        function) is smaller than scf_conv}. Set in the form of 10−N.
        :type scf_conv: int, default: 8

        :param num_threads: Number of threads (and thus CPU cores) for 
        Gaussian to use. Unless num_threads is explicitly specified, 
        Gaussian will only use one thread (run on one core). 
        :type num_threads: int, default: 1

        :param executable: Optional name of the Gaussian executable. 
        Note that if a string is specified then it is a fatal error 
        if that executable is not found. Default searches for 'g16',
        'g09' or 'g03' path executables in that order.
        :type executable: str, default: None

        :param use_template: Determine whether or not to use a 
        user-provided template file for running external programs.
        If use_template = 1, a file named 'gau_job.tpl' must be in
        run directory. The file should only contain the route section 
        of a Gaussian input file. The route section defines the method 
        to be used and SCF convergence criteria. Do not include any 
        information about coordinates or point charge treatment since 
        this will all be handled by sander. Also, do not include any 
        Link 0 Commands (line starting with %) since these are handled 
        by sander. 
        :type use_template: int: 0 or 1, default: 0 (off)

        :param ntpr: Controls frequency of printing for dipole moment 
        to file gau_job.dip Defaults to &cntrl namelist variable ntpr.
        :type ntpr: int, default: None

        :param dipole: Toggles writing of dipole moment to file gau_job.dip.
        :type dipole: int: 0 or 1, default: 0 (off)

        :param mem: String that specifies how much memory Gaussian should 
        be allowed to use. 
        :type mem: str, default: '256MB'

        Example resulting input section 
        --------------------------------

            &gau
                method = 'BP86',
                basis  = '6-31G**',
                num_threads = 8,
                mem='1GB',
            /

        Example gau_job.tpl contents
        -----------------------------

            #P B3LYP/6-31G* SCF=(Conver=8)
        """
        pass

    def amber_orca(basis='SV(P)', 
                   cbasis='NONE', 
                   jbasis='NONE', 
                   method='blyp', 
                   convkey='VERYTIGHTSCF',
                   grid=4,
                   finalgrid=6,
                   maxiter=100,
                   maxcore=1024,
                   num_threads=1,
                   use_template=0,
                   ntpr=None,
                   dipole=0
                   ):
        """
        Currently under development.

        Parameters as described in the Amber23 Manual. Note the default 
        values match the default Amber settings and are not to be understood
        as best practice recommendations.
        
        :param basis: Basis set type to be used in the calculation. 
        Possible choices include svp, 6-31g, etc. See Orca manual for a 
        complete list. (Default: basis = "SV(P)")
        :type basis: str, default: 'SV(P)'

        :param cbasis: Auxiliary basis set for correlation fitting. See 
        Orca manual for a complete list. (Default: basis = "NONE")
        :type cbasis: str, default: 'NONE'
    
        :param jbasis: Auxiliary basis set for Coulomb fitting. See 
        Orca manual for a complete list. (Default: basis = "NONE")
        :type jbasis: str, default: 'NONE'

        :param method: Method to be used in the calculation. Popular 
        choices include hf, pm3, blyp, and mp2. 
        :type method: str, default: 'blyp'

        :param convkey: General SCF convergence setting for simplified 
        Orca input. Can take values 'TIGHTSCF', 'VERYTIGHTSCF', etc. 
        :type convkey: str, default: 'VERYTIGHTSCF'

        :param grid: Grid type used during the SCF for the XC quadrature 
        in DFT. (Default: grid = 4, this corresponds to Intacc = 4.34 
        for the radial grid and an angular Lebedev grid with 302 points. 
        Conservatively chosen together with finalgrid to conserve energy.)
        :type grid: int, default: 4

        :param finalgrid: Grid type used for the energy and gradient 
        calculation after the SCF for the XC quadrature in DFT. 
        (Default: finalgrid = 6, this corresponds to Intacc = 5.34 
        for the radial grid and an angular Lebedev grid with 590 points. 
        Conservatively chosen together with grid to conserve energy.)
        :type finalgrid: int, default: 6

        :param maxiter: Maximum number of SCF iteractions. 
        :type maxiter: int, default: 100

        :param maxcore: Global scratch memory (in MB) used by Orca. You 
        may need to increase this when running larger jobs. See Orca 
        manual for more information. 
        :type maxcore: int, default: 1024

        :param num_threads: Number of threads (and thus CPU cores) for
        Orca to use. Note that Orca only supports OpenMPI. 
        :type num_threads: int, default: 1

        :param use_template: Determine whether or not to use a user-provided 
        template file for running external programs. The template file for 
        Orca should be named orc_job.tpl and must at least contain keywords 
        specifying the method and basis set to be used in the calculation, 
        See example below.
        :type use_template: int, default: 0 (off)

        :param ntpr: Controls frequency of printing for the dipole moment to 
        file orc_job.dip Defaults to &cntrl namelist variable ntpr.
        :type ntpr: int

        :param dipole: Toggles writing of the dipole moment to file 
        orc_job.dip 
        :type dipole: int, default: 0 (off)
                   
        Example resulting input section 
        --------------------------------

            &orc
                method = 'blyp',
                basis  = 'svp',
            /

        Example orc_job.tpl contents
        -----------------------------

            # ORCA input file for BLYP/SVP simulation
            ! BLYP SVP
        """
        pass

    def amber_MM():
        """
        Currently under development.
        """
        pass
        

    #def QM_MM(self):
    #def QM_xTB(self):
