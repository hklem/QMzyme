###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@lsu.edu
###############################################################################

"""
Module containing functions utilized by truncation_schemes.py.
"""

import numpy as np
from QMzyme.QMzymeAtom import QMzymeAtom
from QMzyme.converters import *
from QMzyme.data import backbone_atoms
from QMzyme.MDAnalysisWrapper import find
import warnings

def cap_H(replace_atom, fixed_atom, name=None, bond_length=1.09, base_atom=None, residue=None):
    """
    Builds a new hydrogen atom that caps a truncation cut point, replacing
    `replace_atom` and positioned along the bond vector from `fixed_atom` at
    the given "bond_length".

    :param replace_atom: The atom being replaced by the new hydrogen.
    :type replace_atom: :class:`~QMzyme.QMzymeAtom.QMzymeAtom`
    :param fixed_atom: The atom bonded to `replace_atom` that anchors the bond
        vector calculation; the new hydrogen is placed along the line from
        `fixed_atom` to `replace_atom`, at `bond_length` from `fixed_atom`.
    :type fixed_atom: :class:`~QMzyme.QMzymeAtom.QMzymeAtom`
    :param name: Atom name for the new hydrogen atom. If None, defaults to
        `H{replace_atom.element}` (e.g. 'HC' when replacing a carbon).
    :type name: str, optional, default=None
    :param bond_length: Bond length (Angstroms) used to place the new hydrogen
        along the `fixed_atom`-`replace_atom` vector.
    :type bond_length: float, default=1.09
    :param base_atom: Atom whose attributes (resid, resname, segid, charge, etc.)
        are copied onto the new atom via `create_new_atom`. If None, defaults to
        `replace_atom` except when `fixed_atom` is the backbone N of a Proline residue,
        in which case `base_atom` defaults to `fixed_atom` instead.
    :type base_atom: :class:`~QMzyme.QMzymeAtom.QMzymeAtom`, optional, default=None
    :param residue: If provided, records that this residue was capped via the
        'cap_H' scheme by setting `residue.region._residue_capping_scheme[residue.resid] = 'cap_H'`.
    :type residue: QMzymeResidue, optional, default=None
    :return: The newly created hydrogen atom.
    :rtype: :class:`~QMzyme.QMzymeAtom.QMzymeAtom`
    """
    new_position = set_bond_length(replace_atom.position, fixed_atom.position, bond_length)
    if name is not None:
        new_atom_dict = {
                'element': 'H',
                'type': 'H',  
                'name': f'{name}',
                'position': new_position, 
                'mass': 1.00794,
            }
    else:
        new_atom_dict = {
                'element': 'H',
                'type': 'H',  
                'name': f'H{replace_atom.element}',
                'position': new_position, 
                'mass': 1.00794,
            }
    if base_atom is None:
        base_atom = replace_atom
        if fixed_atom.resname == 'PRO' and fixed_atom.name == backbone_atoms['N']:
            new_atom_dict['charge'] = 0.0
            new_atom_dict['name'] = 'HN'
            base_atom = fixed_atom
    new_atom = create_new_atom(base_atom, new_atom_dict) # used fixed atom because sometimes replaced atom comes from a different residue
    
    if residue is not None:
        residue.region._residue_capping_scheme[residue.resid] = 'cap_H'
    return new_atom


def cap_ACE(replace_atom, universe):
    """
    Builds the atoms of an ACE (acetyl) capping group from the residue
    preceding `replace_atom`, capping the N-terminal cut point of
    `replace_atom`'s residue. This is done by copying C, O, CA of the
    preceding residue, and methyl-capping the CA. CA of the preceding
    residue is renamed to CH3, and hydrogens attached to it is renamed
    as methyl hydrogens ('HH31', 'HH32', 'HH33'). All cap atoms have
    designated charge of 0.

    :param replace_atom: An atom of the residue being truncated.
    :type replace_atom: :class:`~QMzyme.QMzymeAtom.QMzymeAtom`
    :param universe: MDAnalysis universe used to select the preceding
        residue's atoms.
    :type universe: MDAnalysis.core.universe.Universe

    :raises ValueError: If no atoms exist at the preceding resid (i.e.
        `replace_atom`'s residue has no preceding residue in the same segid),
        or if the preceding residue is missing backbone N, C, O, or CA.

    :return: New atoms making up the ACE cap (C, O, CH3, and up to 3 methyl
        hydrogens).
    :rtype: list of :class:`~QMzyme.QMzymeAtom.QMzymeAtom`

    .. note:: If the cap ends up with other than 3 methyl hydrogens (unusual
        atom naming in the preceding residue), a `UserWarning` is issued
        rather than an exception being raised.
    """
    preceding_resid = replace_atom.resid - 1
    prev_atoms = universe.select_atoms(f"segid {replace_atom.segid} and resid {preceding_resid}")

    # Examines if preceding residue contains any atoms, and raises ValueError if there are no atoms in preceding residue
    if len(prev_atoms) == 0:
        raise ValueError(
            f"Cannot build ACE cap: no preceding residue found at resid {preceding_resid}."
        )

    # Checks for the presence of backbone atoms in the previous residue
    prevN_raw  = find(prev_atoms, backbone_atoms['N'])
    prevC_raw  = find(prev_atoms, backbone_atoms['C'])
    prevO_raw  = find(prev_atoms, backbone_atoms['O'])
    prevCA_raw = find(prev_atoms, backbone_atoms['CA'])
    if any(atom is None for atom in (prevN_raw, prevC_raw, prevO_raw, prevCA_raw)):
        raise ValueError(
            f"Cannot build ACE cap: preceding residue {preceding_resid} is missing "
            "backbone N, C, O, or CA."
        )
    
    # Checks for the presence of CB. It will return None for glycine
    prevCB_raw   = find(prev_atoms, 'CB')
    prev_HA_raws = [atom for atom in prev_atoms if atom.name in ('HA', 'HA1', 'HA2', 'HA3')]

    cap_atoms = []

    # Set resname of the C/O backbone atoms to ACE
    for mda_atom, out_name in ((prevC_raw, 'C'), (prevO_raw, 'O')):
        atom = mda_atom_to_qmz_atom(mda_atom)
        atom.name = out_name
        atom.resname = "ACE"
        cap_atoms.append(atom)

    # Add CA as CH3 for ACE capping residue
    new_CA = mda_atom_to_qmz_atom(prevCA_raw)
    new_CA.name = 'CH3'
    new_CA.resname = "ACE"
    cap_atoms.append(new_CA)
    prevCA = mda_atom_to_qmz_atom(prevCA_raw)

    methyl_h_index = 0

    # Convert the HA atoms to methyl hydrogens
    for mda_atom in prev_HA_raws:
        methyl_h_index += 1
        atom = mda_atom_to_qmz_atom(mda_atom)
        atom.name = f"HH3{methyl_h_index}"
        atom.resname = "ACE"
        cap_atoms.append(atom)

    # If there is a sidechain, do cap_H
    capping_substituents = [a for a in (prevN_raw, prevCB_raw) if a is not None]
    for substituent_raw in capping_substituents:
        methyl_h_index += 1
        substituent = mda_atom_to_qmz_atom(substituent_raw)
        new_H = cap_H(substituent, prevCA, f"HH3{methyl_h_index}")
        new_H.resid = substituent_raw.resid
        new_H.resname = "ACE"
        cap_atoms.append(new_H)

    # Set atom charges of the cap_atoms as 0
    for atom in cap_atoms:
        atom.charge = 0

    # Raise warning if the number of hydrogens in the methyl group is not 3
    methyl_H_count = sum(1 for a in cap_atoms if a.name.startswith('HH3'))
    if methyl_H_count != 3:
        warnings.warn(
            f"cap_ACE built only {methyl_H_count}/3 methyl H atoms for resid "
            f"{replace_atom.resid}. Check preceding residue atom naming.",
            UserWarning,
        )
    
    return cap_atoms

def cap_NME(replace_atom, universe):
    """
    Builds the atoms of an NME (N-methylamine) capping group from the residue
    following `replace_atom`, capping the C-terminal cut point of
    `replace_atom`'s residue. This is done by copying N, H, CA of the
    following residue, and methyl-capping the CA. CA of the following
    residue is renamed to CH3, and hydrogens attached to it is renamed
    as methyl hydrogens ('HH31', 'HH32', 'HH33'). All cap atoms have
    designated charge of 0.

    :param replace_atom: An atom of the residue being truncated.
    :type replace_atom: :class:`~QMzyme.QMzymeAtom.QMzymeAtom`
    :param universe: MDAnalysis universe used to select the following
        residue's atoms.
    :type universe: MDAnalysis.core.universe.Universe

    :raises ValueError: If no atoms exist at the following resid (i.e.
        `replace_atom`'s residue has no following residue in the same segid),
        if the following residue is missing backbone N, C, or CA, or if the
        following residue is Proline (no amide H to cap).

    :return: New atoms making up the NME cap (N, H, CH3, and up to 3 methyl
        hydrogens).
    :rtype: list of :class:`~QMzyme.QMzymeAtom.QMzymeAtom`

    .. note:: If the cap ends up with other than 3 methyl hydrogens (unusual
        atom naming in the following residue), a `UserWarning` is issued
        rather than an exception being raised.
    """
    following_resid = replace_atom.resid + 1
    next_atoms = universe.select_atoms(f"segid {replace_atom.segid} and resid {following_resid}")

    # Examines if following residue contains any atoms, and raises ValueError if there are no atoms in following residue
    if len(next_atoms) == 0:
        raise ValueError(
            f"Cannot build NME cap: no following residue found at resid {following_resid}."
        )

    # Checks for the presence of backbone atoms in the following residue
    nextN_raw  = find(next_atoms, backbone_atoms['N'])
    nextH_raw  = find(next_atoms, backbone_atoms['H'])  # None for Proline
    nextC_raw  = find(next_atoms, backbone_atoms['C'])
    nextCA_raw = find(next_atoms, backbone_atoms['CA'])

    if nextH_raw is None:
        resname = nextN_raw.residue.resname if nextN_raw is not None else None
        if resname == 'PRO':
            raise ValueError(
                f"Cannot build NME cap: residue {following_resid} is Proline. "
                "Proline's backbone N has no amide H due to it bonded to the ring "
                "CD carbon instead. This results in the side chain constraining the"
                "backbone geometry in a way an NME cap can't represent. Choose a"
                "different cut point or handle this residue some other way."
            )
    
    if any(atom is None for atom in (nextN_raw, nextC_raw, nextCA_raw)):
        raise ValueError(
            f"Cannot build NME cap: following residue {following_resid} is missing "
            "backbone N, C, or CA."
        )
    
    # Checks for the presence of CB. It will return None for glycine
    nextCB_raw   = find(next_atoms, 'CB')
    next_HA_raws = [atom for atom in next_atoms if atom.name in ('HA', 'HA1', 'HA2', 'HA3')]

    cap_atoms = []

    # Set resname of the N/H backbone atoms to NME
    new_N = mda_atom_to_qmz_atom(nextN_raw)
    new_N.name = 'N'
    new_N.resname = "NME"
    cap_atoms.append(new_N)

    if nextH_raw is not None:
        new_amide_H = mda_atom_to_qmz_atom(nextH_raw)
        new_amide_H.name = 'H'
        new_amide_H.resname = "NME"
        cap_atoms.append(new_amide_H)

    # Add CA as CH3 for NME capping residue
    new_CA = mda_atom_to_qmz_atom(nextCA_raw)
    new_CA.name = 'CH3'
    new_CA.resname = "NME"
    cap_atoms.append(new_CA)
    nextCA = mda_atom_to_qmz_atom(nextCA_raw)

    methyl_h_index = 0

    # Convert the HA atoms to methyl hydrogens
    for mda_atom in next_HA_raws:
        methyl_h_index += 1
        atom = mda_atom_to_qmz_atom(mda_atom)
        atom.name = f"HH3{methyl_h_index}"
        atom.resname = "NME"
        cap_atoms.append(atom)

    # If there is a sidechain, do cap_H
    capping_substituents = [a for a in (nextC_raw, nextCB_raw) if a is not None]
    for substituent_raw in capping_substituents:
        methyl_h_index += 1
        substituent = mda_atom_to_qmz_atom(substituent_raw)
        new_H = cap_H(substituent, nextCA, f"HH3{methyl_h_index}")
        new_H.resid = substituent_raw.resid
        new_H.resname = "NME"
        cap_atoms.append(new_H)

    # Set atom charges of the cap_atoms as 0
    for atom in cap_atoms:
        atom.charge = 0
        
    # Raise warning if the number of hydrogens in the methyl group is not 3
    methyl_H_count = sum(1 for a in cap_atoms if a.name.startswith('HH3'))
    if methyl_H_count != 3:
        warnings.warn(
            f"cap_NME built only {methyl_H_count}/3 methyl H atoms for resid "
            f"{replace_atom.resid}. Check preceding residue atom naming.",
            UserWarning,
        )

    return cap_atoms

def add_removed_atoms(residue):
    """
    Revives all atoms currently tracked in `removed_atoms` by instantiating 
    new QMzymeAtom objects using metadata from the MDAnalysis universe,
    then appends them back to both this residue and the parent region's master list.

    :return: A list of the revived QMzymeAtom objects.
    :rtype: list
    """
    removed = list(residue.removed_atoms)
    for atom in removed:
        residue.atoms.append(atom)
        residue.region.add_atom(atom)
    residue.removed_atoms = []

def remove_added_atoms(residue):
    """
    Permanently removes the 'added_atoms' from both this residue's 
    active atoms list and the parent region's master atom list.

    :return: A list of the QMzymeAtom objects that were purged.
    :rtype: list
    """
    added = list(residue.added_atoms)
    added_ids = {atom.id for atom in added}
    residue.atoms = [a for a in residue.atoms if a.id not in added_ids]
    for atom in added:
        residue.region.remove_atom(atom)
    residue.added_atoms = []


def balance_charge(region, truncated_region):
    """
    To be used if original region has atom charge information, 
    so the QMzymeRegion guess_charge() and read_charges() methods do not get messed up.
    Within a residue, any atoms added should distribute the charges of any atoms removed. 
    """
    for res in region.residues:
        truncated_res = truncated_region.get_residue(res.resid)
        removed_atoms = []
        added_atoms = []
        chrg = 0
        for atom in res.atoms:
            if atom.name not in [atom.name for atom in truncated_res.atoms]:
                removed_atoms.append(atom)
        for atom in truncated_res.atoms:
            if atom.name not in [atom.name for atom in res.atoms]:
                added_atoms.append(atom)
        if len(removed_atoms) == 0:
            continue
        for atom in removed_atoms:
            chrg += atom.charge
        fractional_charge = chrg/len(added_atoms)
        for atom in added_atoms:
            atom.charge = fractional_charge


def set_bond_length(mobile_coords, fixed_coords, new_length):
    M = new_length/np.linalg.norm(fixed_coords-mobile_coords)
    new_coords = fixed_coords-(M*(fixed_coords-mobile_coords))
    return new_coords


def create_new_atom(base_atom, new_atom_dict):
    for key, val in base_atom.__dict__.items():
        if key.startswith('_QMzymeAtom__'):
            key = key.split('_QMzymeAtom__')[-1]
        if key not in new_atom_dict:
            new_atom_dict[key] = val
    new_atom = QMzymeAtom(**new_atom_dict)
    return new_atom


def get_preceding_Catom(region, resid):
    if resid == 1:
        return None
    if region._universe != None:
        try:
            mda_atom = region._universe.select_atoms(f'resid {resid-1} and name {backbone_atoms["C"]}').atoms[0]
            atom = mda_atom_to_qmz_atom(mda_atom)
        except: #Possibly a structure that had been previous truncated.. or weird backbone atom names!
            return None
    return atom


def get_following_Natom(region, resid):
    if region._universe != None:
        try:
            mda_atom = region._universe.select_atoms(f'resid {resid+1} and name {backbone_atoms["N"]}').atoms[0]
            atom = mda_atom_to_qmz_atom(mda_atom)
        except: #Possibly a structure that had been previous truncated.. or weird backbone atom names!
            atom = None # covers if res is last protein res in universe
    return atom


def has_Nterm_neighbor(atom):
    return atom.resid-1 in atom.region.resids


def has_Cterm_neighbor(atom):
    return atom.resid+1 in atom.region.resids
