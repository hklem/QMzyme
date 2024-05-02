import numpy as np


def set_bond_length(mobile_coords, fixed_coords, new_length):
    M = new_length/np.linalg.norm(fixed_coords-mobile_coords)
    new_coords = fixed_coords-(M*(fixed_coords-mobile_coords))
    return new_coords

