###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
Module containing functions to define a QMzymeRegion based on some logic/workflow.
"""

import QMzyme.MDAnalysisWrapper as MDAwrapper
from QMzyme.RegionBuilder import RegionBuilder

selection_schemes = {
    'distance_cutoff': lambda model, cutoff: distance_cutoff(model, cutoff),
}

def distance_cutoff(model, cutoff, include_whole_residues=True):
    """
    Function to select QMzymeRegion based on ``distance_cutoff``. Ad hoc approach
    that is often used to create multiple models that will be evaluated for 
    converge of relevant properties to validate region choice.
    """
    if not model.has_region('catalytic_center'):
        raise UserWarning("You must first define a catalytic_center. See method `set_catalytic_center()`.")
    


    neighbors = MDAwrapper.get_neighbors(
        model.universe.select_atoms('all'),
        model.get_region('catalytic_center').get_atom_group(),
        cutoff)
    
    neighbors_byres = neighbors.residues.sorted_unique.atoms

    #neighbors = model.universe.select_atoms(f"around {distance_cutoff}")

    # Convert MDAnalysis AtomGroup to QMzymeRegion
    qmz_neighbors = RegionBuilder(name=f"cutoff_{cutoff}", atom_group=neighbors).get_region()
    qmz_neighbors_byres = RegionBuilder(name=f"cutoff_{cutoff}", atom_group=neighbors_byres).get_region()

    for atom in qmz_neighbors_byres.atoms:
        if atom.id in qmz_neighbors.ids:
            atom.set_neighbor(True)
        else:
            atom.set_neighbor(False)

    region = qmz_neighbors_byres
    if include_whole_residues is False:
        region = qmz_neighbors

    setattr(region, 'catalytic_center', model.get_region('catalytic_center'))
    return region
        