###############################################################################
# Code written by Heidi Klem while at
# Colorado State University as a graduate student
# in the Paton and McCullagh groups and at the
# National Institute of Standards and Technology
# as an NRC Postdoc (Fed).
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

import numpy as np

ef_Vm = 5.14220674763e11 # a.u. = 1 V/m
conv = ef_Vm*1e-6/100 # to MV/cm

def electric_field_efg(log_file, atom1_id, atom2_id):
    '''
    To calculate the electric field from a Gaussian prop=efg calculation.

    Parameters
    ----------
    log_file : str, required
        Path to Gaussian log file.
    atom1_id : int, required
        0 indexed id for starting atom.
    atom2_id : int, required
        0 indexed id for ending atom.

    Notes
    ----- 
    - Atom order will alter sign of electric field. Typical to define atom1 as the less electronegative atom, 
    such that the bond vector will be in the direction of the bond dipole.
    
    Returns
    -------
    Scalar electric field value in units MV/cm averaged over the electric field magnitude at atom1 and atom2.

    '''

    with open(log_file) as f:
        data = f.readlines()
    center_start = None
    ef_sart = None

    for i, line in enumerate(data):
        if center_start == None:
            if 'Atomic Center' in line:
                center_start = i
        if '-------- Electric Field --------' in line:
            ef_start = i+3
            break
    
    coords1 = np.array([float(x) for x in data[center_start+atom1_id].split()[-3:]])
    coords2 = np.array([float(x) for x in data[center_start+atom2_id].split()[-3:]])
    
    unit_vec =  unit_vector(coords1,coords2)

    ef1 = np.array([float(x) for x in data[ef_start+atom1_id].split()[-3:]])
    ef2 = np.array([float(x) for x in data[ef_start+atom2_id].split()[-3:]])

    return (0.5*(np.dot(ef1, unit_vec)+np.dot(ef2, unit_vec)))*conv # in units MV/cm

def unit_vector(coords_start, coords_end):
    '''
    To calculate the normalized vector between two points in cartesian space.

    Parameters
    ----------
    coords_start : array or list of floats, required
        Starting point coordinates
    coords_end : array or list of floats, required
        Ending point coordinates

    Returns 
    -------
    unit vector array
    '''
    if type(coords_start) is list:
        coords_start = np.array(coords_start)
    if type(coords_end) is list:
        coords_end = np.array(coords_end)
    bond_length = np.linalg.norm(coords_end-coords_start)
    bond_vec = coords_end-coords_start
    return bond_vec/bond_length