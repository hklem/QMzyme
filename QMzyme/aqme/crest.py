#!/usr/bin/python

#######################################################################
#             CREST functions from the CSEARCH module                 #
#######################################################################

from __future__ import print_function, absolute_import
import os
import glob
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdmolfiles
from rdkit import Geometry
import subprocess
import rdkit
from pathlib import Path
import shutil
from QMzyme.aqme.utils import read_file, run_command
from QMzyme.aqme.filter import geom_filter
from rdkit.Chem import rdMolTransforms


def atompairs(mol, atom1, atom2, constraints):
    active = []
    for x in constraints:
        active.append(x[:2])

    for i, x in enumerate(active):
        active[i] = [int(j) for j in x]

    pairs = []
    bonds = [(x.GetBeginAtomIdx(), x.GetEndAtomIdx()) for x in mol.GetBonds()]
    for [a, b] in bonds:
        if [a + 1, b + 1] not in active and [b + 1, a + 1] not in active:
            dist = round(rdMolTransforms.GetBondLength(mol.GetConformer(), a, b), 3)
            at_a, at_b = mol.GetAtoms()[a].GetSymbol(), mol.GetAtoms()[b].GetSymbol()
            if atom1 == "X" and atom2 == "X":
                pairs.append([float(a + 1), float(b + 1), dist])
            elif atom1 == "X" and atom2 == "H":
                if at_a == "H" or at_b == "H":
                    pairs.append([float(a + 1), float(b + 1), dist])
            else:
                if (at_a == atom1 and at_b == atom2) or (
                    at_a == atom2 and at_b == atom1
                ):
                    pairs.append([float(a + 1), float(b + 1), dist])
    return pairs


def get_constraint(mol, constraints):
    # constrained optimization with xtb
    xx_pairs = atompairs(mol, "X", "X", constraints)
    all_fix = []
    for x in constraints:
        all_fix.append(list(x))
    for x in xx_pairs:
        all_fix.append(x)
    return all_fix


def xyzall_2_xyz(xyzin, name):
    # converting multiple xyz to single
    command_run_1 = ["obabel", xyzin, "-oxyz", "-O" + name + "_conf_.xyz", "-m"]
    subprocess.run(command_run_1, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def xtb_opt_main(
    name,
    dup_data,
    dup_data_idx,
    self,
    charge,
    mult,
    smi,
    constraints_atoms,
    constraints_dist,
    constraints_angle,
    constraints_dihedral,
    method_opt,
    geom,
    complex_ts=False,
    mol=None,
    name_init=None,

):

    """
    Run xTB using subprocess to perform CREST/CREGEN conformer sampling
    """

    name_no_path = os.path.basename(Path(name)).split(".xyz")[0]

    # folder to create the files
    if self.args.destination is None:
        if method_opt == 'crest':
            csearch_dir = Path(self.args.w_dir_main) / "CSEARCH"
        elif method_opt == 'xtb':
            csearch_dir = Path(self.args.w_dir_main) / "CMIN"
    else:
        if self.args.initial_dir.as_posix() in f"{self.args.destination}":
            csearch_dir = Path(self.args.destination)
        else:
            csearch_dir = Path(self.args.initial_dir).joinpath(self.args.destination)

    # create the initial xyz input
    if method_opt == 'crest':
        self.args.log.write(f"\no  Starting xTB pre-optimization before CREST sampling")
        dat_dir = csearch_dir / "crest_xyz"
        xyzin = f"{dat_dir}/{name_no_path}.xyz"
    elif method_opt == 'xtb':
        rdmolfiles.MolToXYZFile(mol, f"{name}.xyz")
        self.args.log.write(f"\no  Starting xTB optimization")
        dat_dir = csearch_dir / "xtb_xyz"
        xyzin = f"{dat_dir}/{name_no_path}_xtb.xyz"
    dat_dir.mkdir(exist_ok=True, parents=True)
    shutil.move(f"{name}.xyz", xyzin)

    os.environ["OMP_STACKSIZE"] = self.args.stacksize
    # to run xTB/CREST with more than 1 processor
    os.environ["OMP_NUM_THREADS"] = str(self.args.nprocs)
    opt_valid = True

    os.chdir(dat_dir)
    
    # for systems that were created from 1D and 2D inputs (i.e. SMILES), this part includes two xTB
    # constrained optimizations to avoid geometry problems in noncovalent complexes and transition states

    # xTB optimization with all bonds frozen
    constrained_opt = False
    if len(constraints_atoms) > 0 or len(constraints_dist) > 0 or len(constraints_angle) > 0 or len(constraints_dihedral) > 0:
        constrained_opt = True
        complex_ts = True

    xyzoutxtb1 = str(dat_dir) + "/" + name_no_path + "_xtb1.xyz"
    if complex_ts:
        all_fix = get_constraint(mol, constraints_dist)

        _ = create_xcontrol(
            self.args,
            constraints_atoms,
            all_fix,
            [],
            [],
            xyzin,
            "constrain1.inp",
        )

        command1 = [
            "xtb",
            xyzin,
            "--opt",
            "--input",
            "constrain1.inp",
            "-c",
            str(charge),
            "--uhf",
            str(int(mult) - 1),
            "-P",
            str(self.args.nprocs),
        ]

        if self.args.xtb_keywords is not None:
            for keyword in self.args.xtb_keywords.split():
                command1.append(keyword)

        xtb_out1 = f'{os.path.dirname(Path(xyzoutxtb1))}/{os.path.basename(Path(xyzoutxtb1)).split(".xyz")[0]}'
        run_command(command1, f"{xtb_out1}.out")
        try:
            os.rename(str(dat_dir) + "/xtbopt.xyz", xyzoutxtb1)
        except FileNotFoundError:
            os.rename(str(dat_dir) + "/xtblast.xyz", xyzoutxtb1)

        # remove files that might interfere in subsequent calculations (i.e. wrong electron readings)
        for file in glob.glob('*') + glob.glob('*.*') + glob.glob('.*'):
            if os.path.exists(file):
                if file.find('_xtb2') == -1 and file.find('_xtb1') == -1 and file.find('.out') == -1:
                    try:
                        os.remove(file)
                    except IsADirectoryError:
                        try:
                            shutil.rmtree(file)
                        except OSError: # this avoids problems when running AQME in HPCs
                            pass

        if constrained_opt:
            xyzoutxtb2 = str(dat_dir) + "/" + name_no_path + "_xtb2.xyz"
            # xTB optimization with the user-defined constraints
            _ = create_xcontrol(
                self.args,
                list(constraints_atoms),
                list(constraints_dist),
                list(constraints_angle),
                list(constraints_dihedral),
                xyzoutxtb1,
                "constrain2.inp",
            )

            command2 = [
                "xtb",
                xyzoutxtb1,
                "--opt",
                "--input",
                "constrain2.inp",
                "-c",
                str(charge),
                "--uhf",
                str(int(mult) - 1),
                "-P",
                str(self.args.nprocs),
            ]

            if self.args.xtb_keywords is not None:
                for keyword in self.args.xtb_keywords.split():
                    command2.append(keyword)

            xtb_out2 = f'{os.path.dirname(Path(xyzoutxtb2))}/{os.path.basename(Path(xyzoutxtb2)).split(".xyz")[0]}'
            run_command(command2, f"{xtb_out2}.out")

            try:
                os.rename(str(dat_dir) + "/xtbopt.xyz", xyzoutxtb2)
            except FileNotFoundError:
                os.rename(str(dat_dir) + "/xtblast.xyz", xyzoutxtb2)
        else:
            xyzoutxtb2 = xyzoutxtb1

    else:
        # Preoptimization with xTB to avoid issues from innacurate starting structures in CREST.
        # If you're dealing with a large system, increase the stack size
        try:
            command = [
                "xtb",
                xyzin,
                "--opt",
                "-c",
                str(charge),
                "--uhf",
                str(int(mult) - 1),
                "-P",
                str(self.args.nprocs), 
            ]

            if self.args.xtb_keywords is not None:
                for keyword in self.args.xtb_keywords.split():
                    command.append(keyword)
            xtb_out1 = f'{os.path.dirname(Path(xyzin))}/{os.path.basename(Path(xyzin)).split(".xyz")[0]}'
            run_command(command, f"{xtb_out1}_xtb1.out")

            os.rename(str(dat_dir) + "/xtbopt.xyz", xyzoutxtb1)
        except FileNotFoundError:
            self.args.log.write(f"\nx  There was an error during the xTB pre-optimization. This error is related to parallelization of xTB jobs and is normally observed when using metal complexes in some operative systems/OpenMP versions. AQME is switching to using one processor (nprocs=1).\n")
            self.args.nprocs = 1
            try:
                xtb_out1 = f'{os.path.dirname(Path(xyzin))}/{os.path.basename(Path(xyzin)).split(".xyz")[0]}'
                if self.args.xtb_keywords is None:
                    comm_xtb = f"export OMP_STACKSIZE={self.args.stacksize} && export OMP_NUM_THREADS={self.args.nprocs},1 \
                    && xtb {xyzin} --opt -c {charge} --uhf {int(mult) - 1} >> {xtb_out1}_xtb1.out"
                else:
                    comm_xtb = f"export OMP_STACKSIZE={self.args.stacksize} && export OMP_NUM_THREADS={self.args.nprocs},1 \
                    && xtb {xyzin} --opt -c {charge} --uhf {int(mult) - 1} {self.args.xtb_keywords} >> {xtb_out1}_xtb1.out"
                subprocess.call(comm_xtb, shell=False)
                os.rename(str(dat_dir) + "/xtbopt.xyz", xyzoutxtb1)
            except FileNotFoundError:
                if self.args.program.lower() == "crest":
                    self.args.log.write(f"\nx  There was another error during the xTB pre-optimization that could not be fixed. Trying CREST directly with no xTB preoptimization.\n")
                else:
                    self.args.log.write(f"\nx  There was another error during the xTB pre-optimization that could not be fixed (this molecule will be skipped).\n")
                opt_valid = False
                mol_rd = None

        xyzoutxtb2 = xyzoutxtb1

    xyzoutall = str(dat_dir) + "/" + name_no_path + "_conformers.xyz"

    # CREST sampling
    if self.args.program.lower() == "crest":
        self.args.log.write(f"\no  Starting CREST sampling")
        if constrained_opt:
            _ = create_xcontrol(
                self.args,
                list(constraints_atoms),
                list(constraints_dist),
                list(constraints_angle),
                list(constraints_dihedral),
                xyzoutxtb2,
                ".xcontrol.sample",
            )

        command = [
            "crest",
            xyzoutxtb2,
            "--chrg",
            str(charge),
            "--uhf",
            str(int(mult) - 1),
            "-T",
            str(self.args.nprocs),
            "--ewin",
            str(self.args.ewin_csearch),
        ]

        if constrained_opt:
            command.append("--cinp")
            command.append(".xcontrol.sample")
            const_command = command.copy()

        if self.args.crest_keywords is not None:
            for keyword in self.args.crest_keywords.split():
                command.append(keyword)
        try:
            run_command(command, f"/{dat_dir}/{name_no_path}.out")
            natoms = open("crest_best.xyz").readlines()[0].strip()
        except FileNotFoundError:
            self.args.log.write(f"\nx  CREST optimization failed! This might be caused by different reasons:\n   1) In metal complexes: using metal complexes without specifying any kind of template in the complex_type option (i.e. squareplanar).\n   2) In TSs: include the \"--noreftopo\" option in CREST with the crest_keywords option (i.e. crest_keywords=\"--noreftopo\").\n   3) In big systems: increase stacksize with the stacksize option (i.e. stacksize=\"4GB\").")
            try:
                self.args.log.write(f"\no  Trying the CREST calculations with stacksize=\"4GB\".")
                os.environ["OMP_STACKSIZE"] = '4GB'
                run_command(command, f"/{dat_dir}/{name_no_path}.out")
                natoms = open("crest_best.xyz").readlines()[0].strip()
            except FileNotFoundError:
                self.args.log.write(f"\nx  CREST optimization failed again even with stacksize=\"4GB\"! Contact the administrators to check this issue in more detail.\n")
                if constrained_opt and "--noreftopo" not in command:
                    try:
                        self.args.log.write(f"\no  Constraints were detected, trying a new CREST run with --noreftopo. WARNING! Check that your geometry doesn't isomerize!\n")
                        if self.args.crest_keywords is not None:
                            for keyword in self.args.crest_keywords.split():
                                const_command.append(keyword)
                        const_command.append('--noreftopo')
                        run_command(const_command, f"/{dat_dir}/{name_no_path}.out")
                        natoms = open("crest_best.xyz").readlines()[0].strip()
                    except FileNotFoundError:
                        self.args.log.write(f"\nx  CREST optimization failed again even with --noreftopo! Contact the administrators to check this issue in more detail.\n")
                        opt_valid = False
                else:
                    opt_valid = False

        # CREGEN sorting
        try:
            if self.args.cregen and int(natoms) != 1 and opt_valid:
                self.args.log.write(f"\no  Starting CREGEN sorting")
                command = ["crest", "crest_best.xyz", "--cregen", "crest_conformers.xyz"]

                if self.args.cregen_keywords is not None:
                    for keyword in self.args.cregen_keywords.split():
                        command.append(keyword)

                run_command(command, f"{dat_dir}/{name_no_path}_cregen.out")
        except UnboundLocalError:
            pass

        try:
            if opt_valid:
                if os.path.exists(str(dat_dir) + "/crest_clustered.xyz"):
                    shutil.copy(str(dat_dir) + "/crest_clustered.xyz", xyzoutall)

                elif os.path.exists(str(dat_dir) + "/crest_ensemble.xyz"):
                    shutil.copy(str(dat_dir) + "/crest_ensemble.xyz", xyzoutall)
                else:
                    shutil.copy(str(dat_dir) + "/crest_conformers.xyz", xyzoutall)
        except FileNotFoundError:
            self.args.log.write("\nx  CREST conformer sampling failed! Please, try other options (i.e. include constrains, change the crest_keywords option, etc.)")
            opt_valid = False

    if opt_valid:
        if self.args.program.lower() == "crest":
            xyzall_2_xyz(xyzoutall, name_no_path)
            xyz_files = glob.glob(name_no_path + "_conf_*.xyz")
        if self.args.program.lower() == "xtb":
            xyz_files = [xyzoutxtb1]
        for _, file in enumerate(xyz_files):
            name_conf = f'{os.path.basename(Path(file)).split(".xyz")[0]}'
            command_xyz = ["obabel", "-ixyz", file, "-osdf", "-O" + name_conf + ".sdf"]
            subprocess.run(command_xyz, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        if self.args.program.lower() == "crest":
            sdwriter = Chem.SDWriter(str(f"{csearch_dir}/{name_no_path}.sdf"))
        sdf_files = glob.glob(name_no_path + "*.sdf")
        for file in sdf_files:
            mol = rdkit.Chem.SDMolSupplier(file, removeHs=False, sanitize=False)
            mol_rd = rdkit.Chem.RWMol(mol[0])
            file_nopath = f'{os.path.basename(Path(file)).split(".sdf")[0]}'
            if self.args.program.lower() == "xtb":
                # convert from hartree (default in xtb) to kcal
                energy_Eh = float(open(f'{file_nopath}.xyz', "r").readlines()[1].split()[1])
                energy_kcal = energy_Eh*627.5
                mol_rd.SetProp("_Name", name_init)
                os.remove(file)
            elif self.args.program.lower() == "crest":
                # convert from hartree (default in xtb) to kcal
                energy_Eh = float(open(file, "r").readlines()[0])
                energy_kcal = str(energy_Eh*627.5)
                mol_rd.SetProp("_Name", name_no_path)
                mol_rd.SetProp("Energy", energy_kcal)
                mol_rd.SetProp("Real charge", str(charge))
                mol_rd.SetProp("Mult", str(int(mult)))
                if smi is not None:
                    mol_rd.SetProp("SMILES", str(smi))
                mol_geom = Chem.Mol(mol_rd)
                passing_geom = geom_filter(self,mol_geom,geom)
                if passing_geom:
                    sdwriter.write(mol_rd)
                os.remove(file)
                os.remove(f'{file_nopath}.xyz')
        if self.args.program.lower() == "crest":
            sdwriter.close()
    else:
        xyz_files = []

    # remove xTB/CREST files to avoid wrong readings of molecular information
    for file in glob.glob('*') + glob.glob('*.*') + glob.glob('.*'):
        if os.path.exists(file):
            if file.find('.out') == -1:
                if self.args.program.lower() == "xtb":
                    try:
                        os.remove(file)
                    except OSError: # this avoids problems when running AQME in HPCs
                        pass
                elif self.args.program.lower() == "crest":
                    if file.find('_xtb2') == -1 and file.find('_xtb1') == -1 and file.find('.out') == -1:
                        try:
                            if file == 'crest_clustered.xyz':
                                os.rename('crest_clustered.xyz', f"{dat_dir}/{name_no_path}_clustered.xyz")
                            else:
                                try:
                                    os.remove(file)
                                except OSError: # this avoids problems when running AQME in HPCs
                                    pass
                        except IsADirectoryError:
                            try:
                                shutil.rmtree(file)
                            except OSError: # this avoids problems when running AQME in HPCs
                                pass

    os.chdir(self.args.w_dir_main)

    if method_opt == 'crest':
        dup_data.at[dup_data_idx, "crest-conformers"] = len(xyz_files)
        return 1

    if method_opt == 'xtb':
        return mol_rd, energy_kcal, opt_valid


def create_xcontrol(
    args,
    constraints_atoms,
    constraints_dist,
    constraints_angle,
    constraints_dihedral,
    xyzin,
    name_constraint,
):
    """
    Function to create the .xcontrol.sample if constraints are defined
    """

    constrained_sampling = False

    unique_atoms = []
    for atom in constraints_atoms:
        if atom not in unique_atoms:
            unique_atoms.append(int(atom))
    for x in constraints_dist:
        for i in x[:2]:
            if i not in unique_atoms:
                unique_atoms.append(int(i))
    for x in constraints_angle:
        for i in x[:3]:
            if i not in unique_atoms:
                unique_atoms.append(int(i))
    for x in constraints_dihedral:
        for i in x[:4]:
            if i not in unique_atoms:
                unique_atoms.append(int(i))

    if len(unique_atoms) > 0:

        constrained_sampling = True

        # call --constrain just fo create a coord.ref file
        subprocess.run(
            ["crest", xyzin, "--constrain", "1"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        # add the constraints part
        edited_xcontrol = "$constrain\n"

        if constraints_atoms != []:
            edited_xcontrol += "atoms: "
            edited_xcontrol += ",".join(str(atom_idx) for atom_idx in constraints_atoms)
            edited_xcontrol += "\n"

        for constraint_type in [
            constraints_dist,
            constraints_angle,
            constraints_dihedral,
        ]:
            if constraint_type != []:
                for const in constraint_type:
                    if constraint_type == constraints_dist:
                        edited_xcontrol += "distance: "
                        edited_xcontrol += (
                            ",".join(str(int(val)) for val in const[:2])
                            + ","
                            + str(const[2])
                        )
                    elif constraint_type == constraints_angle:
                        edited_xcontrol += "angle: "
                        edited_xcontrol += (
                            ",".join(str(int(val)) for val in const[:3])
                            + ","
                            + str(const[3])
                        )
                    elif constraint_type == constraints_dihedral:
                        edited_xcontrol += "dihedral: "
                        edited_xcontrol += (
                            ",".join(str(int(val)) for val in const[:4])
                            + ","
                            + str(const[4])
                        )
                    edited_xcontrol += "\n"

        edited_xcontrol += f"force constant={args.crest_force}\n"
        edited_xcontrol += "reference=coord.ref\n"

        # metadyn part
        if name_constraint == ".xcontrol.sample":
            outlines = read_file(os.getcwd(), os.getcwd(), xyzin)
            n_atoms = int(outlines[0])
            edited_xcontrol += "$metadyn\n"
            edited_xcontrol += "atoms: "
            if n_atoms == 1: # just to avoid bugs when parsing single atoms
                edited_xcontrol += "atoms: 1"
            else:
                # if the list is too long, CREST doesn't read it when called from subprocess() in Python
                # I need to include ranges to shorten the lists of atoms for the $metadyn section
                new_cycle = True
                start = True
                for atom_idx in range(1, n_atoms + 1):
                    if new_cycle:
                        new_cycle = False
                        start_idx = atom_idx
                    if atom_idx not in unique_atoms:
                        if start: # just in case the first atom isn't part of the list
                            start = False
                            start_idx = atom_idx
                        elif atom_idx == n_atoms:
                            if start_idx == (atom_idx):
                                edited_xcontrol += f"{start_idx}"
                            else:
                                edited_xcontrol += f"{start_idx}-{atom_idx}"
                    elif atom_idx in unique_atoms and not start:
                        new_cycle = True
                        if start_idx == (atom_idx-1):
                            edited_xcontrol += f"{start_idx}"
                        else:
                            edited_xcontrol += f"{start_idx}-{atom_idx-1}"
                        if atom_idx != n_atoms:
                            edited_xcontrol += ','

        edited_xcontrol += "\n$end\n"

        # write the file
        xcontrol_file = open(name_constraint, "w")
        xcontrol_file.write(edited_xcontrol)
        xcontrol_file.close()

    return constrained_sampling


def nci_ts_mol(
    smi,
    log,
    seed,
    constraints_atoms,
    constraints_dist,
    constraints_angle,
    constraints_dihedral,
):
    using_const = None
    if constraints_atoms is not None:
        using_const = constraints_atoms
        constraints_atoms = [[float(y) for y in x] for x in constraints_atoms]
        constraints_atoms = np.array(constraints_atoms)
    if constraints_dist is not None:
        using_const = constraints_dist
        constraints_dist = [[float(y) for y in x] for x in constraints_dist]
        constraints_dist = np.array(constraints_dist)
    if constraints_angle is not None:
        using_const = constraints_angle
        constraints_angle = [[float(y) for y in x] for x in constraints_angle]
        constraints_angle = np.array(constraints_angle)
    if constraints_dihedral is not None:
        using_const = constraints_dihedral
        constraints_dihedral = [[float(y) for y in x] for x in constraints_dihedral]
        constraints_dihedral = np.array(constraints_dihedral)

    molsH = []
    mols = []
    for m in smi:
        mols.append(Chem.MolFromSmiles(m))
        molsH.append(Chem.AddHs(Chem.MolFromSmiles(m)))

    for m in molsH:
        Chem.EmbedMultipleConfs(m, numConfs=1, randomSeed=seed)
    for m in mols:
        Chem.EmbedMultipleConfs(m, numConfs=1, randomSeed=seed)

    coord = [0.0, 0.0, 5.0]
    molH = molsH[0]
    for _, fragment in enumerate(molsH[1:]):
        offset_3d = Geometry.Point3D(coord[0], coord[1], coord[2])
        molH = Chem.CombineMols(molH, fragment, offset_3d)
        coord[1] += 5
        Chem.SanitizeMol(molH)

    coord = [0.0, 0.0, 5.0]
    mol = mols[0]
    for _, fragment in enumerate(mols[1:]):
        offset_3d = Geometry.Point3D(coord[0], coord[1], coord[2])
        mol = Chem.CombineMols(mol, fragment, offset_3d)
        coord[1] += 5
        Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)

    try:
        mol = Chem.ConstrainedEmbed(mol, molH, randomseed=seed)
    except ValueError: # removingH in the core (molH) avoids problems related to RDKit embedding
        molH = Chem.RemoveHs(molH)
        try:
            mol = Chem.ConstrainedEmbed(mol, molH, randomseed=seed)
        except ValueError:
            log.write(f"\nx  Constrained optimization failed due to an embedding problem with RDKit that could not be fixed!")
            pass

    atom_map = []
    for atom in mol.GetAtoms():
        atom_map.append(atom.GetAtomMapNum())

    max_map = max(atom_map)
    for a in mol.GetAtoms():
        if a.GetSymbol() == "H":
            max_map += 1
            a.SetAtomMapNum(int(max_map))

    adapted_atoms = []
    adapted_dist = []
    adapted_angle = []
    adapted_dihedral = []
    # assign constraints
    if using_const is not None:
        for smi_part in smi:
            if ':' not in smi_part or '[' not in smi_part: # for SMILES that are not mapped
                log.write(f"\nx  Constraints were specified {using_const} but atoms might not be mapped in the SMILES input!")
                adapted_atoms = constraints_atoms
                adapted_dist = constraints_dist
                adapted_angle = constraints_angle
                adapted_dihedral = constraints_dihedral
                break

        if constraints_atoms is not None:
            for _, ele in enumerate(constraints_atoms):
                if adapted_atoms != []: # for mapped SMILES
                    for atom in mol.GetAtoms():
                        if ele == atom.GetAtomMapNum():
                            adapted_atoms.append(float(atom.GetIdx()) + 1)
                else:
                    adapted_atoms = constraints_atoms
            adapted_atoms = np.array(adapted_atoms)

        constraints = [constraints_dist, constraints_angle, constraints_dihedral]
        adapted_consts = [adapted_dist, adapted_angle, adapted_dihedral]
        n_consts= [2, 3, 4]
        for const,adapted_const,n_const in zip(constraints,adapted_consts,n_consts):
            if adapted_const == []: # for mapped SMILES
                for _, r in enumerate(const):
                    nr = []
                    for _, ele in enumerate(r[:n_const]):
                        if const is not None:
                            for atom in mol.GetAtoms():
                                if ele == atom.GetAtomMapNum():
                                    nr.append(float(atom.GetIdx()) + 1)
                    nr.append(r[-1])
                    adapted_const.append(nr)
                adapted_const = np.array(adapted_const)

    return (
        mol,
        adapted_atoms,
        adapted_dist,
        adapted_angle,
        adapted_dihedral,
    )