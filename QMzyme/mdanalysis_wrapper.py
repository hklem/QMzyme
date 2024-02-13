###############################################################################
# Code written by Heidi Klem while at
# Colorado State University as a graduate student
# in the Paton and McCullagh groups and at the
# National Institute of Standards and Technology
# as an NRC Postdoc (Fed).
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

try:
    import MDAnalysis as mda
except:
    pass

def res_selection(pdb=None, sel=''):
    '''
    Parameters
    ----------
    pdb : str, required
        Filename of structure. The default is None.
    sel : str, required
        Selection command used by MDAnalysis. The default is ''.

    Returns
    -------
    res_dict : list of dictionaries
        Dictionary containing the residue name, residue number, and chain.

    '''
    u = mda.Universe(pdb)
    if 'chain' in sel:
        sel = sel.split('chain')
        sel = sel[0]+' segid '+sel[1]
    sel = u.select_atoms(sel)
    res = [x for x in sel.residues]
    res_list = []
    for r in res:
        res_dict = {
        'Residue name': r.resname,
        'Residue number': int(r.resnum),
        'Chain': r.segid
        }
        res_list.append(res_dict)
    # just added this right now because if res_list has multiple residues this will not work in the current form of generate.py... so that would need to be fixed first.
    if len(res_list) == 1:
        res_list = res_dict
    return res_list

def atom_selection(pdb=None, sel=''):
    '''
    Parameters
    ----------
    pdb : str, required
        Filename of structure. The default is None.
    sel : str, required
        Selection command used by MDAnalysis. The default is ''.

    Returns
    -------
    idx : list of atom indices (zero indexed).

    '''
    u = mda.Universe(pdb)
    if 'chain' in sel:
        sel = sel.split('chain')
        sel = sel[0]+' segid '+sel[1]
    sel = u.select_atoms(sel)
    idx = [x-1 for x in sel.ids]

    return idx
