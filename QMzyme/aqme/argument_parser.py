#####################################################.
#      This file contains the argument parser         #
#####################################################.

import os

var_dict = {
    "varfile": None,
    "verbose": True,
    "input": "",
    "output_name": "output",
    "command_line": False,
    "name": None,
    "path": "",
    "output": ".sdf",
    "csearch": False,
    "cmin": False,
    "qprep": False,
    "qcorr": False,
    "smi": None,
    "metal_atoms": [],
    "auto_metal_atoms": True,
    "charge": None,
    "mult": None,
    "complex_coord": [],
    "complex_type": "",
    "metal_idx": [],
    "metal_sym": [],
    "constraints_atoms": [],
    "constraints_dist": [],
    "constraints_angle": [],
    "constraints_dihedral": [],
    "ewin_cmin": 5.0,
    "ewin_csearch": 5.0,
    "opt_fmax": 0.05,
    "opt_steps": 1000,
    "opt_steps_rdkit": 1000,
    "heavyonly": True,
    "degree": 120.0,
    "max_torsions": 0,
    "sample": "auto",
    "auto_sample": 20,
    "ff": "MMFF",
    "seed": 62609,
    "rms_threshold": 0.25,
    "max_matches_rmsd": 1000,
    "energy_threshold": 0.25,
    "initial_energy_threshold": 0.0001,
    "max_mol_wt": 0,
    "ani_method": "ANI2x",
    "stacksize": "1G",
    "xtb_keywords": None,
    "max_workers": 4,
    "ewin_sample_fullmonte": 2.0,
    "ewin_fullmonte": 5.0,
    "nsteps_fullmonte": 100,
    "nrot_fullmonte": 3,
    "ang_fullmonte": 30,
    "cregen": False,
    "cregen_keywords": None,
    "program": None,
    "nprocs": 8,
    "mem": "16GB",
    "mol": None,
    "destination": None,
    "qm_input": "",
    "ts_input": "opt=(calcfc,noeigen,ts,maxstep=5)",
    "qm_end": "",
    "chk_path": "",
    "gen_atoms": [],
    "bs_nogen": "",
    "bs_gen": "",
    "freeze_atoms": [],
    "lowest_only": False,
    "lowest_n": None,
    "e_threshold_qprep": None,
    "chk": False,
    "w_dir_main": os.getcwd(),
    "files": [],
    "atom_types": [],
    "cartesians": [],
    "dup": True,
    "dup_threshold": 0.0001,
    "ro_threshold": 0.1,
    "amplitude_ifreq": 0.2,
    "ifreq_cutoff": 0.0,
    "freq_conv": None,
    "im_freq_input": 'opt=(calcfc,maxstep=5)',
    "s2_threshold": 10.0,
    "isom_type": None,
    "isom_inputs": os.getcwd(),
    "vdwfrac": 0.5,
    "covfrac": 1.1,
    "fullcheck": True,
    "suffix": "",
    "geom": [],
    "bond_thres": 0.2,
    "angle_thres": 30,
    "dihedral_thres": 30,
    "crest_keywords": None,
    "crest_force": 0.5,
    "prefix": "",
    "qdescp": False,
    "qdescp_temp": 300,
    "qdescp_acc": 0.2,
    "qdescp_solvent": None,
    "boltz": True,
    "nmr_atoms": [6, 1],  # [C,H]
    "nmr_slope": [-1.0537, -1.0784],  # [C,H]
    "nmr_intercept": [181.7815,31.8723],  # [C,H]
    "nmr_experim": None,
    "nodup_check": False,
    "qdescp_atoms": [],
    "xtb_opt": True,
    "dbstep_r": 3.5,
    "robert": True,
    "csv_name": None,
    "crest_nrun": 1,
    "crest_nclust": 0.4
}


# part for using the options in a script or jupyter notebook
class options_add:
    pass


def set_options(kwargs):
    # set default options and options provided
    options = options_add()
    # dictionary containing default values for options

    for key in var_dict:
        vars(options)[key] = var_dict[key]
    for key in kwargs:
        if key in var_dict:
            vars(options)[key] = kwargs[key]
        elif key.lower() in var_dict:
            vars(options)[key.lower()] = kwargs[key.lower()]
        else:
            print("Warning! Option: [", key,":",kwargs[key],"] provided but no option exists, try the online documentation to see available options for each module.",)

    return options
