###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@lsu.edu
###############################################################################

"""
Module containing functions to truncate a QMzymeRegion based on some logic/scheme.
"""

from QMzyme.data import protein_residues, backbone_atoms
from QMzyme.QMzymeRegion import QMzymeRegion
from QMzyme.truncation_utils import *
from QMzyme.CalculateModel import CalculateModel
import abc


class TruncationScheme(abc.ABC):
    """
    Base class for defining QM region truncation methodologies.
 
    A truncation scheme defines the specific rule for how protein
    amino acid residues in the QMzymeModel are truncated to create
    the desired QM region.
    
    Before any residue is truncated, three checks are run as methods called during initiation of the TruncationScheme:

    1. `_check_gly_ala`: identifies Gly/Ala residues that would become small organic groups (methane or ethane) once truncated, and
       resolves them using the remove_methane, remove_ethane, or
       extend_gly_ala_backbone arguments.
    2. `_check_override_truncation`: identifies residues that have already
       been truncated by a prior TruncationScheme call, and resolves them according to
       the override_truncation argument.
    3. `_check_override_capping`: identifies residues that have already been
       capped (e.g. cap_H, cap_ACE, cap_NME), and resolves them using the
       override_capping argument.
    
    Each check either raises a ValueError describing how to resolve the issue, appends the relevant residues to skip_resids so
    they are excluded from truncation, or passes silently if no action is needed.

    TruncationScheme subclasses are responsible for determining how the
    residues are truncated (e.g. :class:`~QMzyme.TruncationSchemes.TerminalAlphaCarbon`,
    :class:`~QMzyme.TruncationSchemes.AlphaCarbon`, and
    :class:`~QMzyme.TruncationSchemes.BetaCarbon`).

    .. image:: ../../docs/Images/truncation_scheme.png
        :width: 650
    """

    _gly_ala_check = None

    """
    Class attribute that controls how `_check_gly_ala` identifies isolated
    Gly/Ala residues for a given TruncationScheme subclass.

    If `_gly_ala_check` is "terminal", only Gly/Ala residues with no neighboring
    residue (resid +/- 1) present in the region are flagged as isolated
    (e.g. TerminalAlphaCarbon).

    If `_gly_ala_check` is "all", every Gly/Ala residue within the selection is
    flagged as isolated, regardless of whether neighboring residues are present
    (e.g. AlphaCarbon).

    If `_gly_ala_check` is None, no Gly/Ala isolation check is performed for the
    subclass (e.g. BetaCarbon).
    """
 
    def __init__(self, region, name, selection = None, remove_methane:bool = None, remove_ethane:bool = None, extend_gly_ala_backbone:bool = False, override_truncation:bool = None, override_capping:bool = None):
 
        if name is not None:
            self.region = QMzymeRegion(name=name, atoms=list(region.atoms), universe=region._universe)
 
            self.region._residue_truncation_attr = dict(region._residue_truncation_attr)
            self.region._residue_capping_attr = dict(region._residue_capping_attr)
            self.region._residue_method_attr = dict(getattr(region, '_residue_method_attr', {}))
            if hasattr(region, 'method') and region.method is not None:
                self.region.method = region.method
            if hasattr(region, 'charge'):
                self.region.charge = region.charge
            # Record which scheme was used and (if a new region was created) its parent region
            self.region._selection_attr['parent_region'] = region
        else:
            self.region = region
 
        self.truncated_region = None
        self.remove_methane = remove_methane
        self.remove_ethane = remove_ethane
        self.selection = selection if selection is not None else 'all'
        self.override_truncation = override_truncation
        self.override_capping = override_capping
        self.extend_gly_ala_backbone = extend_gly_ala_backbone

        selection_region = region._universe.select_atoms(self.selection)
        self.residues_to_truncate = set(int(a.resid) for a in selection_region.residues)
        self.skip_resids = set()
 
        self._check_gly_ala()
        self._check_override_truncation()
        self._check_override_capping()
 
        # Part for actual truncation! Residues that are no part of residues_to_truncate, skip_resids, and protein_residues are skipped
        for res in list(self.region.residues):
            if res.resid not in self.residues_to_truncate:
                continue
            elif res.resid in self.skip_resids:
                continue
            elif res.resname not in protein_residues:
                continue

            self.truncate(res)
            res.truncation_params = self.__class__.__name__

        # Checking if all residues are truncated
        protein_resid = [res for res in self.region.residues if res.resname in protein_residues]
        all_processed = all(
            (getattr(res, 'truncation_params', None) is not None or
             getattr(res, 'capping_scheme', None) is not None)
            for res in protein_resid
        ) 

        # Designates truncated attribute to the region if all protein residues are either truncated or capped
        if all_processed:
            self.region.truncated = True

        # Raise warning if the region is not truncated
        if not getattr(self.region, 'truncated', True):
            print("Warning: not all residues were truncated or capped. 'truncated' attribute not set.")
 
    def _check_gly_ala(self):
        """
        Identify Gly/Ala residues whose truncation would create a small organic group,
        methane (Gly) or ethane (Ala) fragment, and resolve them according to the kwargs remove_methane, remove_ethane, and extend_gly_ala_backbone arguments.
 
        If a group is present and its corresponding flag (remove_methane for Gly,
        remove_ethane for Ala) is None and extend_gly_ala_backbone is False,
        the workflow returns a ValueError, containing resolution protocols.
        
        If self.remove_methane and/or self.remove_ethane are True, the corresponding Gly and/or Ala
        residues, respectively, are removed from the region if the self.scheme TruncationScheme class
        would result in them becoming free methane and/or ethane groups.

        If self.remove_methane/remove_ethane is False, the corresponding Gly and/or Ala residues,
        respecitve, are kept within the region, and a warning message is raised.
        
        If extend_gly_ala_backbone is True, all isolated Gly/Ala residues are capped with
        acetyl (ACE) and N-methyl amide (NME) groups, using add_N_terminus_ACE and
        add_C_terminus_NME methods from QMzymeRegion.py.
 
        No exception is raised and no action is taken if no Gly/Ala residues
        in the selection would be isolated by truncation.
        """

        if self._gly_ala_check is None:
            return

        # region_resids is the set of resid values actually present in the region
        region_resids = {int(res.resid) for res in self.region.residues}
        organic_group_map = {"gly": "methane", "ala": "ethane"}

        gly_methane_region = []
        ala_ethane_region = []
        isolated_gly = []
        isolated_ala = []

        # Determine which resids we should consider for Gly/Ala isolation
        for res in self.region.residues:
            resname = res.resname.lower()
            resid = int(res.resid)

            # If the residues are not within residues_to_truncate, skip
            if resid not in self.residues_to_truncate:
                continue
            
            # If the resname is not Gly or Ala, skip
            elif resname not in ("gly", "ala"):
                continue
            
            # Assign appropriate ethane or methane label
            organic_group = organic_group_map[resname]

            # Add additional conditions here for further truncation scheme here in case Ala and Gly are isolated in other ways!

            # For TerminalAlphaCarbon, consider that Gly and Ala are isolated when there is no resid +1 or -1
            if self._gly_ala_check == "terminal":
                is_isolated = (
                    (resid - 1 not in region_resids)
                    and (resid + 1 not in region_resids)
                )

            # For Ca truncation, consider that Gly and Ala are isolated when Gly and Ala exists
            elif self._gly_ala_check == "all":
                is_isolated = True

            # If none of these cases, Gly and Ala are not considered isolated
            else:
                is_isolated = False

            # If there is no islated Gly and ala, skip
            if not is_isolated:
                continue
            
            # If there is isolated Gly and Ala, append them to their own list
            if resname == "gly":
                gly_methane_region.append(res)
                isolated_gly.append((res, organic_group))

            if resname == "ala":
                ala_ethane_region.append(res)
                isolated_ala.append((res, organic_group))

        # If there is no isolated Gly or Ala, return
        if not isolated_gly and not isolated_ala:
            return
 
        # Determine, per organic group, if they need to be treated for potential chemical inaccuracy
        # A group only needs a decision if it is isolated Gly/Ala and remove_methane/remove_ethane is None and extend backbone is False
        needs_decision_gly = bool(gly_methane_region) and self.remove_methane is None and self.extend_gly_ala_backbone is False
        needs_decision_ala = bool(ala_ethane_region) and self.remove_ethane is None and self.extend_gly_ala_backbone is False
 
        # If there is any isolated Gly and Ala that needs user's influence, raise ValueError
        if needs_decision_gly or needs_decision_ala:
            pending = []

            if needs_decision_gly:
                pending.extend(isolated_gly)

            if needs_decision_ala:
                pending.extend(isolated_ala)

            msg = ""
            for res, org_group in pending:
                msg += f"Truncation of Residue {res} would result in a(n) {org_group}\n"
            msg += (
                "Representation of the native residue may not be an appropriate chemical substitute.\n\n"
                "RESOLUTION PROTOCOL\n"
                "To resolve this, either:\n"
                "1) Remove them from your GenerateModel instance, run:\n"
                "      GenerateModel.truncate(remove_methane=True, remove_ethane=True)\n"
                "2) Set either flag to False to keep the small organic group.\n"
                "3) Alternatively, you can also include the neighboring residues using\n"
                "      QMzymeRegion.add_residue(QMzymeResidue)\n"
                "      and apply the TerminalAlphaCarbon scheme.\n"
                "4) Or set extend_gly_ala_backbone=True to add ACE and NME capping to isolated Gly/Ala residue(s).\n\n"
                "Please set remove_methane and/or remove_ethane to True or False,\n"
                "add neighboring residues, or set extend_gly_ala_backbone to True."
            )

            raise ValueError(msg)

        # If remove_methane is True, remove the isolated Gly
        if self.remove_methane is True:
            for res in gly_methane_region:
                self.region.remove_residue(res)
 
        # If remove_ethane is True, remove the isolated Ala
        if self.remove_ethane is True:
            for res in ala_ethane_region:
                self.region.remove_residue(res)
        
        # If remove_methane is False, return a warning message
        if self.remove_methane is False and isolated_gly:
            for res, org_group in isolated_gly: 
                print(f"Truncation of Residue {res} would result in a {org_group}")
            print("These organic groups may not be an appropriate representation of the active site region.\n")

        # If remove_ethane is False, return a warning message
        if self.remove_ethane is False and isolated_ala:
            for res, org_group in isolated_ala:
                print(f"Truncation of Residue {res} would result in an {org_group}")
            print("These organic groups may not be an appropriate representation of the active site region.\n")

        # If extend_gly_ala_backbone is True and there is isolated Gly or Ala, go through the protocol
        if self.extend_gly_ala_backbone is True and (isolated_gly or isolated_ala):

            # If it is not TerminalAlphaCarbon or AlphaCarbon, raise ValueError
            if not isinstance(self, (TerminalAlphaCarbon, AlphaCarbon)):
                raise ValueError("Currently, QMzyme only supports extend_gly_ala_backbone with TerminalAlphaCarbon and AlphaCarbon TruncationScheme subclasses.")
            
            else:
                # Create a list of capped gly ala to be skipped during override_capping
                self._capped_gly_ala = set()

                # Create a list of residues to cap
                to_cap = [int(res.resid) for res, _ in (isolated_gly + isolated_ala) if int(res.resid) in self.residues_to_truncate]

                if not to_cap:
                    return

                if isinstance(self, (AlphaCarbon)):
                    both_neighbors_present = [self.region.get_residue(resid) for resid in to_cap if (resid - 1) in region_resids and (resid + 1) in region_resids]
                    if both_neighbors_present:
                        raise UserWarning(
                            f"extend_gly_ala_backbone=True requested for residue(s) {both_neighbors_present},"
                            "but both neighboring residues are already present in the region. Please remove these"
                            "residue(s) using QMzymeRegion.remove_residue() method or manually extend backbone"
                            "using QMzymeRegion.add_N_terminus_ACE() and QMzymeRegion.add_C_terminus_NME() methods"
                            )
                
                # Acquires method from the original region
                composite_types = {'QMQM2', 'QMXTB', 'QMChargeField'}
                calc_type = CalculateModel.calc_type
                is_multiscale = calc_type in composite_types
                secondary_type = calc_type[2:] if is_multiscale else None
                sync_targets = {'QM', secondary_type} if is_multiscale else set()
                orig_method = getattr(self.region, "method", None)
                self.region.method = None

                try:
                    for resid in to_cap:
                        # Precompute once per residue (orig_res used only to inspect pre-cap flags)
                        orig_res = next((r for r in self.region.residues if int(r.resid) == resid), None)

                        # Get the preexisting attribute from the original resiude
                        inherited_type = getattr(orig_res, 'method_type', None)

                        # Track whether this resid actually gained caps/bridge
                        cap_success = False

                        # Apply ACE and NME caps (attempt both)
                        for cap_call in (self.region.add_N_terminus_ACE, self.region.add_C_terminus_NME):
                            before = {int(r.resid) for r in self.region.residues}
                            try:
                                cap_call(resid)
                            except Exception as e:
                                # If there is an exception, return the error message
                                print(f"     {cap_call.__name__}({resid}): skipped/failed: {e}")
                                continue
                            
                            after = {int(r.resid) for r in self.region.residues}
                            new_resids = after - before

                            # If any new residues were added (ACE/NME/bridge), consider that progress
                            if new_resids:
                                cap_success = True

                            # Build lookup once (after caps were added)
                            res_lookup = {int(r.resid): r for r in self.region.residues}

                            # Inherit method_type (if used)
                            if inherited_type:
                                for nr in new_resids:
                                    if nr in res_lookup:
                                        res_lookup[nr].method_type = inherited_type

                                # Sync composite regions if applicable
                                if inherited_type in sync_targets:
                                    target_region = CalculateModel.calculation[inherited_type]
                                    for nr in new_resids:
                                        if nr in res_lookup:
                                            for atom in res_lookup[nr].atoms:
                                                if atom.id not in target_region.ids:
                                                    target_region.add_atom(atom)

                        # Additional check: if ACE+NME cap atoms are present around resid, treat as success
                        has_ace_before = any(a.resname == 'ACE' for a in self.region.atoms if a.resid == resid - 1)
                        has_nme_after = any(a.resname == 'NME' for a in self.region.atoms if a.resid == resid + 1)
                        if has_ace_before and has_nme_after:
                            cap_success = True

                        # If the backbone of the inbetween residue was inserted, treat as success
                        real_res_present = any(a.resname not in ('ACE', 'NME') for a in self.region.atoms if a.resid == resid)
                        if real_res_present:
                            cap_success = True

                        # Only record resid as capped if we actually made progress
                        if cap_success:
                            self._capped_gly_ala.add(resid)

                finally:
                    self.region.method = orig_method
                
                # Make list of residuesthat are capped
                target_resids = {int(r.resid) for r, _ in (isolated_gly + isolated_ala)}
                actual_capped = {r for r in getattr(self, '_capped_gly_ala', set()) if r in target_resids}

                # Only skip residues that are flanked by ACE (resid-1) and NME (resid+1)
                ace_nme_capped = {
                    r for r in actual_capped
                    if any(a.resname == 'ACE' and a.resid == r - 1 for a in self.region.atoms)
                    and any(a.resname == 'NME' and a.resid == r + 1 for a in self.region.atoms)
                }

                self._capped_gly_ala = actual_capped
                self.skip_resids.update(ace_nme_capped)

                # Set attribute _extend_gly_ala_backbone as True
                for res in self.region.residues:
                    if int(res.resid) in ace_nme_capped:
                        res._extend_gly_ala_backbone = True

    def _check_override_truncation(self):
        """
        Identify residues in the selection that have already been truncated
        (i.e. `QMzymeResidue.trucation_params` is not None and resolve them according
        to the override_truncation argument.
 
        If override_truncation is None and all or part of the region has previously been truncated, raise ValueError with resolution protocol.
        
        If override_truncation is False, the previously truncated residues will be added
        to self.skip_resids and will not be truncated in the current truncation call.

        If override_truncation is True, any previously truncated residues will be
        reverted back to their original state and then the current truncation scheme will be applied.

        Does nothing if no residues were previously truncated.
        """
 
        # Check for residues that are already truncated
        already_truncated = [
            res for res in self.region.residues
            if res.resid in self.residues_to_truncate
            and getattr(res, 'truncation_params', None) is not None
            and res.truncation_params is not None
        ]

        # If nothing is already truncated, return
        if not already_truncated:
            return
    
        # If there is already_truncated residue, and override_truncation is None, raise ValueError
        if self.override_truncation is None:
            msg = ""
            for res in already_truncated:
                msg += f"Residue {res} has already been truncated with {res.truncation_params}\n"
            msg += (
                "These residues will not be re-truncated by default.\n\n"
                "RESOLUTION PROTOCOL\n"
                "To override and re-truncate them, run with:\n"
                "  ... .truncate(scheme=..., override_truncation=True)\n"
                "or set override_truncation=False to silence this warning and skip them.\n\n"
                "Please set override_truncation=True to re-truncate, or False to skip."
            )
            raise ValueError(msg)

        # If override_truncation is False, raise warning print statement
        elif self.override_truncation is False:
            for res in already_truncated:
                print(f"Skipping residue {res}: it has already been truncated with {res.truncation_params}.")
                self.skip_resids.add(res.resid)

        # If override_truncation is True, undo the prior truncation
        elif self.override_truncation is True:
            # Only reset the residues that are already_truncated
            for res in already_truncated:
                remove_added_atoms(res)
                add_removed_atoms(res)
                res.truncation_params = None
                if getattr(res, 'capping_scheme', None) == 'cap_H':
                    res.capping_scheme = None
 
    def _check_override_capping(self):
        """
        Identify residues in the selection that have already been capped
        and resolve them according to the override_capping argument.

        If self.override_capping is None, raise ValueError with resolution protocol.

        If self.override_capping is False, previously capped residues will be retained and ignored during current truncation scheme implementation.

        If self.override_capping is True, previously truncated residues will be
        reverted back to their original state and then the current truncation scheme will be applied.

        Does nothing if no residues in the selection have already been
        capped (excluding any residues already handled by extend_gly_ala_backbone).
        """

        gly_ala_capped = getattr(self, '_capped_gly_ala', set())

        # Precompute cap status for each residue in selection
        already_capped = []
        for res in self.region.residues:
            if res.resid not in self.residues_to_truncate:
                continue
            if res.resid in gly_ala_capped:
                continue
            cap_scheme = getattr(res, 'capping_scheme', None)
            if cap_scheme is not None and cap_scheme != 'cap_H':
                already_capped.append(res)

        if not already_capped:
            return

        # If user hasn't decided, return ValueError
        if self.override_capping is None:
            msg = ""
            for res in already_capped:
                msg += f"Residue {res} has been capped with {res.capping_scheme}\n"
            msg += (
                "These residues will not be re-truncated by default.\n\n"
                "RESOLUTION PROTOCOL\n"
                "To override capping, run with:\n"
                "  ... .truncate(scheme=..., override_capping=True)\n"
                "or set override_capping=False to silence this warning and skip them.\n\n"
                "Please set override_capping=True to re-cap, or False to skip."
            )
            raise ValueError(msg)

        # If override_capping is False, skip only the residues we flagged above
        if self.override_capping is False:
            for res in already_capped:
                has_ace_before = any(a.resname == 'ACE' and a.resid == res.resid - 1 for a in self.region.atoms)
                has_nme_after = any(a.resname == 'NME' and a.resid == res.resid + 1 for a in self.region.atoms)
                fully_capped = (has_ace_before and has_nme_after) or res.capping_scheme == 'cap_H'

                if fully_capped:
                    print(f"Skipping residue {res}: it has already been capped with {res.capping_scheme} or ACE+NME caps present.")
                    self.skip_resids.add(res.resid)
                else:
                    print(f"Residue {res}: only one terminus capped ({res.capping_scheme}). Leaving that side alone; remaining uncapped terminus (if any) will still be truncated.")

        # If override_capping is True, undo the caps and remove ACE/NME cap residues
        if self.override_capping is True:
            cap_resids = set()
            for res in already_capped:
                remove_added_atoms(res)
                add_removed_atoms(res)
                if getattr(res, 'capping_scheme', None) == 'cap_H':
                    res.capping_scheme = None

                cap_resids.update({a.resid for a in self.region.atoms if a.resname in ('ACE', 'NME') and (a.resid == res.resid - 1 or a.resid == res.resid + 1)})

            if cap_resids:
                for res in list(self.region.residues):
                    if res.resid in cap_resids:
                        try:
                            self.region.remove_residue(res)
                        except Exception:
                            self.region.atoms = [a for a in self.region.atoms if a.resid not in cap_resids]
 
    @abc.abstractmethod
    def truncate(self, res):
        ...
 
    def return_region(self):
        return self.region

class TerminalAlphaCarbon(TruncationScheme):
    """
    The TerminalAlphaCarbon scheme will 1) remove N-terminal backbone atoms 
    (N and H) if the preceding sequence residue (resid-1) is not included
    in the region and add a hydrogen atom along the CA-N backbone
    bond vector; and 2) remove C-terminal backbone atoms (C and O) if the 
    following sequence residue (resid+1) is not included in the region and
    add a hydrogen atom along the CA-C backbone bond vector. In
    the case of Proline, if the preceding sequence residue is not present 
    the Proline N atom is kept and a hydrogen is added along the N-(resid-1)C
    backbone bond vector.

    .. image:: ../../docs/Images/terminal_alpha_carbon.png
        :width: 80%
    """

    _gly_ala_check = "terminal"
    def __init__(self, region, name, selection=None, remove_methane:bool = None, remove_ethane:bool = None, extend_gly_ala_backbone:bool = False, override_truncation:bool=None, override_capping:bool = None):
        super().__init__(region, name, selection, remove_methane, remove_ethane, extend_gly_ala_backbone=extend_gly_ala_backbone, override_truncation=override_truncation, override_capping=override_capping)

    def truncate(self, res):
        resname = res.resname
        if resname not in protein_residues:
            return

        Natom = res.get_atom(backbone_atoms['N'])
        CAatom = res.get_atom(backbone_atoms['CA'])
        Catom = res.get_atom(backbone_atoms['C'])
        Oatom = res.get_atom(backbone_atoms['O'])
        preceding_Catom = get_preceding_Catom(self.region, res.resid)
        following_Natom = get_following_Natom(self.region, res.resid)

        preceding_C_present = preceding_Catom is not None and any(
            a.id == preceding_Catom.id and a.name == preceding_Catom.name for a in self.region.atoms
        )
        following_N_present = following_Natom is not None and any(
            a.id == following_Natom.id and a.name == following_Natom.name for a in self.region.atoms
        )

        has_ace_before = getattr(res, 'capping_scheme', None) == 'cap_ACE' and preceding_C_present
        has_nme_after = getattr(res, 'capping_scheme', None) == 'cap_NME' and following_N_present

        if preceding_Catom is not None and not preceding_C_present and not has_ace_before:
            if resname != 'PRO':
                Hatom = res.get_atom(backbone_atoms['H'])
                cap_atom = cap_H(Natom, CAatom, residue=res)
                res.remove_atom(Natom)
                res.remove_atom(Hatom)
                self.region.add_atom(cap_atom)
            else:
                cap_atom = cap_H(preceding_Catom, Natom, residue=res)
                setattr(cap_atom, "id", cap_atom.id - 1)
                self.region.add_atom(cap_atom)

        if following_Natom is not None and not following_N_present and not has_nme_after:
            cap_atom = cap_H(Catom, CAatom, residue=res)
            res.remove_atom(Catom)
            res.remove_atom(Oatom)
            self.region.add_atom(cap_atom)
        
class AlphaCarbon(TruncationScheme):
    """
    Function to truncate a QMzymeRegion according to the AlphaCarbon scheme. 

    .. image:: ../../docs/Images/alpha_carbon.png
        :width: 80%
    """

    _gly_ala_check = "all"
    def __init__(self, region, name, selection=None, remove_methane:bool = None, remove_ethane:bool = None, extend_gly_ala_backbone:bool = False, override_truncation:bool=None, override_capping:bool = None):
        super().__init__(region, name, selection, remove_methane, remove_ethane, extend_gly_ala_backbone=extend_gly_ala_backbone, override_truncation=override_truncation, override_capping=override_capping)

    def truncate(self, res):
        resname = res.resname
        if resname not in protein_residues:
            return

        Natom = res.get_atom(backbone_atoms['N'])
        CAatom = res.get_atom(backbone_atoms['CA'])
        Catom = res.get_atom(backbone_atoms['C'])
        Oatom = res.get_atom(backbone_atoms['O'])
        preceding_Catom = get_preceding_Catom(self.region, res.resid)
        following_Natom = get_following_Natom(self.region, res.resid)

        has_ace_before = any(a.resname == 'ACE' and a.resid == res.resid - 1 for a in self.region.atoms)
        has_nme_after = any(a.resname == 'NME' and a.resid == res.resid + 1 for a in self.region.atoms)

        # N-side truncation: truncate whenever preceding_Catom exists, unless has_ace_before
        if preceding_Catom is not None and not (has_ace_before and getattr(res, 'capping_scheme', None) == 'cap_ACE'):
            if resname != 'PRO':
                Hatom = res.get_atom(backbone_atoms['H'])
                cap_atom = cap_H(Natom, CAatom, residue=res)
                if Natom is not None:
                    res.remove_atom(Natom)
                if Hatom is not None:
                    res.remove_atom(Hatom)
                self.region.add_atom(cap_atom)
            else:
                cap_atom = cap_H(preceding_Catom, Natom, residue=res)
                setattr(cap_atom, "id", cap_atom.id - 1)
                self.region.add_atom(cap_atom)

        # C-side truncation: truncate whenever following_Natom exists, unless has_nme_after
        if following_Natom is not None and (not has_nme_after) and not (has_nme_after and getattr(res, 'capping_scheme', None) == 'cap_NME'):
            cap_atom = cap_H(Catom, CAatom, residue=res)
            if Catom is not None:
                res.remove_atom(Catom)
            if Oatom is not None:
                res.remove_atom(Oatom)
            self.region.add_atom(cap_atom)
            
class BetaCarbon(TruncationScheme):
    """
    The Beta Carbon scheme will keep atoms up to CB, remove all other
    sidechain atoms, and cap the heavy atom(s) directly bonded to CB with
    hydrogen along the CB-X vector. In the case of Proline and Glycine, it
    skips and returns a warning message.

    .. image:: ../../docs/Images/beta_carbon.png
        :width: 80%
    """

    _gly_ala_check = None

    # List of all 
    _CB_neighbors = {
        'SER': ['OG'],
        'CYS': ['SG'],
        'THR': ['OG1', 'CG2'],
        'VAL': ['CG1', 'CG2'],
        'ILE': ['CG1', 'CG2'],
    }

    def __init__(self, region, name, selection=None, remove_methane:bool = None, remove_ethane:bool = None, extend_gly_ala_backbone:bool = False, override_truncation:bool=None, override_capping:bool = None):
        super().__init__(region, name, selection, remove_methane, remove_ethane, extend_gly_ala_backbone=extend_gly_ala_backbone, override_truncation=override_truncation, override_capping=override_capping)

    def truncate(self, res):
        resname = res.resname
        if resname not in protein_residues:
            return

        if resname in ("GLY", "PRO"):
            print(f"{res} will be skipped from the BetaCarbon truncation scheme, due to the residue being {resname}.")
            return

        if resname == "ALA":
            # Already exactly the target fragment: backbone + CB + HB1/HB2/HB3.
            return

        CBatom = res.get_atom('CB')
        heavy_neighbors = self._CB_neighbors.get(resname, ['CG'])
        branched = len(heavy_neighbors) == 2

        keep_names = set(backbone_atoms.values())
        keep_names.add(CBatom.name)
        keep_names.add("HB" if branched else "HB2")
        if not branched:
            keep_names.add("HB3")

        # Remove every sidechain atom except CB, its heavy neighbor(s), and CB's hydrogens
        for atom in list(res.atoms):
            if atom.name in keep_names or atom.name in heavy_neighbors:
                continue
            res.remove_atom(atom)

        # Val/Thr: rename the surviving HB → HB1 before capping
        if branched:
            HBatom = res.get_atom("HB")
            if HBatom is not None:
                HBatom.name = "HB1"

        # Cap each heavy neighbor with HB1/HB2/HB3 as appropriate
        cap_index = 2 if branched else 1
        for name in heavy_neighbors:
            atom = res.get_atom(name)
            if atom is None:
                continue
            cap_atom = cap_H(atom, CBatom, name=f"HB{cap_index}", residue=res)
            res.remove_atom(atom)
            self.region.add_atom(cap_atom)
            cap_index += 1

        # The truncated residue now has exactly an alanine's heavy-atom/H
        # complement (backbone + CB + HB1/HB2/HB3), so relabel it as ALA.
        res.resname = "ALA"
        for atom in res.atoms:
            atom.resname = "ALA"