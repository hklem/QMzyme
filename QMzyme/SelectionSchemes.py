###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
Module containing functions to define a QMzymeRegion based on some logic/workflow.
"""

import QMzyme.MDAnalysisWrapper as MDAwrapper
from QMzyme.RegionBuilder import RegionBuilder
from QMzyme.QMzymeRegion import QMzymeRegion
from QMzyme.QMzymeModel import QMzymeModel
from typing import Type


import abc

class SelectionScheme(abc.ABC):
    """
    SelectionScheme is the abstract class used to prescribe concrete selection scheme
    classes. Fork QMzyme to build your own concrete selection scheme class and submit
    a Pull Request on github to have your scheme added to QMzyme! 

    Below is a template that you can use to design a selection scheme concrete class:

    class InformativeName(SelectionScheme):
        This is an example docstring.

        Include a detailed description of how the scheme works, including references 
        to any relevant literature.

        :Parameters:
            - :model: QMzymeModel to provide starting structure that selection will be performed on.
            - :name: Name of the region generated.

        :Returns:
            :class:`~QMzyme.QMzymeRegion.QMzymeRegion`
        
        :Notes:
            Include any notes you want users to be aware of.
        
        def __init__(self, model, name, {any additional kwargs}):
            Assign any kwargs as attributes to self. Then in your
            select_atoms method you can pull any necessary args from
            self attributes, instead of relying on passing them.

            super().__init__(model, name)
        
        def select_atoms(self):
            Write your code to perform the selection. 

            At the end of your code you should set self.region = {region}. 

            The product of your selection scheme needs to be a QMzymeRegion 
            in order for it to work with GenerateModel().set_region().

        def method_name(self):
            You can add whatever other methods you want in your class, but 
            you should call those methods as necessary in __init__() otherwise
            your scheme will not work in GenerateModel.set_region(). 
    """
    def __init__(self, model, name):
        self.name = name
        self.model: QMzymeModel = model
        self.region: QMzymeRegion
        self.select_atoms()
        self.return_region()

    @abc.abstractmethod
    def select_atoms(self):
        ...

    def return_region(self):
        self.region.rename(self.name)
        return self.region

class DistanceCutoff(SelectionScheme):
    """
    The DistanceCutoff class performs a selection simply based on the distance of 
    atoms from a pre-defined catalytic_center region. Users must first call 
    GenerateModel().set_catalytic_center() in order to then run DistanceCutoff via
    GenerateModel().set_region(selection=DistanceCutoff, name={str}, cutoff={int}). 
    
    This scheme is known to require rather large QM regions to achieve agreement 
    with experiment (Ex., Kulik HJ, Zhang J, Klinman JP, Martínez TJ. How Large 
    Should the QM Region Be in QM/MM Calculations? The Case of Catechol 
    O-Methyltransferase. J Phys Chem B. 2016 Nov 10;120(44):11381-11394. 
    doi: 10.1021/acs.jpcb.6b07814.).


    :param model: QMzymeModel to provide starting structure that selection 
        will be performed on.
    :type model: :class:`~QMzyme.QMzymeModel.QMzymeModel`, required.

    :param name: Name of the region generated.
    :type name: str, required.

    :param cutoff: Numerical value to define cutoff.
    :type cutoff: float, required.

    :param include_whole_residues: Informs code whether or not to only include
        atoms within the cutoff, or include whole residues if they have at least one
        atom within the cutoff. 
    :type include_whole_residues: bool, default=True.


    :returns: :class:`~QMzyme.QMzymeRegion.QMzymeRegion`

    :notes:
        Users are encouraged to evaluate the resulting region. There may be situations where
        a charged residue is within the cutoff distance, however, its charge partner is not. 
        Such situations can drastically alter the chemistry of the model! Maybe someone could
        write up a less generic distance based selection scheme that would take such situations
        into consideration. Or modify the current class to include an argument 
        `include_charge_partners=True`. 
        
    """
    def __init__(self, model, name, cutoff, include_whole_residues=True):
        if name is None:
            name = f'cutoff_{cutoff}'
        self.cutoff = cutoff
        self.include_whole_residues = include_whole_residues
        super().__init__(model, name)

    def select_atoms(self):
        if not self.model.has_region('catalytic_center'):
            raise UserWarning("You must first define a catalytic_center. See method `set_catalytic_center()`.")
        
        neighbors = MDAwrapper.get_neighbors(
        self.model.universe.select_atoms('all'),
        self.model.get_region('catalytic_center')._atom_group, self.cutoff)
    
        neighbors_byres = neighbors.residues.sorted_unique.atoms

        # Convert MDAnalysis AtomGroup to QMzymeRegion
        qmz_neighbors = RegionBuilder(name=f"cutoff_{self.cutoff}", atom_group=neighbors).get_region()
        qmz_neighbors_byres = RegionBuilder(name=f"cutoff_{self.cutoff}", atom_group=neighbors_byres).get_region()

        for atom in qmz_neighbors_byres.atoms:
            if atom.id in qmz_neighbors.ids:
                atom.set_neighbor(True)
            else:
                atom.set_neighbor(False)

        region = qmz_neighbors_byres
        if self.include_whole_residues is False:
            region = qmz_neighbors

        setattr(region, 'catalytic_center', self.model.get_region('catalytic_center'))
        self.region = region

