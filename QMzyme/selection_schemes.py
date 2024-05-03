import QMzyme.MDAnalysisWrapper as MDAwrapper
from QMzyme.QMzymeRegion import QMzymeRegion
from QMzyme.RegionBuilder import RegionBuilder

def within_distance(model, distance_cutoff, include_whole_residues=True):
    if not model.has_region('catalytic_center'):
        raise UserWarning("You must first define a catalytic_center. See method `set_catalytic_center()`.")
    
    neighbors = MDAwrapper.get_neighbors(
        model.starting_structure.select_atoms('all'),
        model.get_region('catalytic_center').get_AtomGroup(),
        distance_cutoff)
    
    neighbors_byres = neighbors.residues.sorted_unique.atoms

    # Convert MDAnalysis AtomGroup to QMzymeRegion
    qmz_neighbors = RegionBuilder(name=f"Atoms_within_cutoff_{distance_cutoff}", atom_group=neighbors).get_region()
    qmz_neighbors_byres = RegionBuilder(name=f"Residues_within_cutoff_{distance_cutoff}", atom_group=neighbors_byres).get_region()

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
        
    



