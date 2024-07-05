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
    SelectionScheme is the abstract base class used to prescribe concrete selection scheme
    sub classes. Fork QMzyme to build your own concrete selection scheme class and submit
    a Pull Request (PR) on github to have your scheme added to QMzyme! You will need to create
    comprehensive tests before the PR is accepted. See the 
    `documentation on contributing <https://qmzyme.readthedocs.io/en/latest/Contributing/index.html>`_ 
    for more information.

    Below is a template that you can use to implement your own selection scheme class:

    class InformativeName(SelectionScheme):
        This is an example docstring.

        Include a detailed description of how the scheme works in the class level doc string. Also
        include any parameters the __init__ method of your class accepts so it will automatically
        documented online (the __init__ method doc string is not documented online by default).

        :param model: QMzymeModel to provide starting structure that selection will be performed on. When
            using the main :class:`~QMzyme.GenerateModel.GenerateModel` class, the QMzyme model is automatically 
            passed as an argument to the selection scheme. It is recommended you use the Universe (universe attribute)
            representing the starting structure to perform the selection on. 
        :type model: :class:`~QMzyme.QMzymeModel.QMzymeModel`
        :param name: Name of the region generated.
        :type name: str, required

        The return should always be a QMzyme region.

        :Returns:
            :class:`~QMzyme.QMzymeRegion.QMzymeRegion`
        
        :Notes:
            Include any notes you want users to be aware of.


        def method_name(self):
            You can add whatever other methods you want in your class, but 
            you should call those methods as necessary in __init__() otherwise
            your scheme will not work in GenerateModel.set_region(). 
    """
    def __init__(self, model, name):
        """
        Assign any key word arguments as attributes to self. Then in your
        select_atoms method you can pull any necessary args from
        self attributes, instead of relying on passing them.

        Every concrete scheme __init__ should include this line at the very end:

        super().__init__(model, name)

        This will automatically run your select_atoms method and return the resulting region.
        """
        self.name = name
        self.model: QMzymeModel = model
        self.region: QMzymeRegion
        self.select_atoms()
        self.reference()
        if self.reference is not None:
            print(f"Use of this selection scheme requires citing the following reference(s): \n \t{self.reference}")
        self.return_region()

    @abc.abstractmethod
    def select_atoms(self):
        """
        Write your code to perform the selection. 

        At the end of your code you should set self.region = {region}. 

        The product of your selection scheme needs to be a QMzymeRegion 
        in order for it to work with GenerateModel().set_region().

        This method is automatically called in the ``super().__init__(model, name)``
        line of your __init__ method.
        """
        ...
    
    def method_name(self):
        """
        You can add whatever other methods you want in your class, but 
        you should call those methods as necessary in __init__() otherwise
        your scheme will be automated in GenerateModel.set_region()
        """
        pass

    def return_region(self):
        """
        This method belongs to the base class and is automatically called in 
        the ``super().__init__(model, name)`` line of your __init__ method. All
        you have to do is make sure you have created a class attribute called `region`.
        """
        self.region.rename(self.name)
        return self.region

    @abc.abstractmethod
    def reference(self):
        """
        This method needs to be included in your class. All it should do is create an
        attribute called `reference` that provides a citable reference of the scheme, to 
        give credit where credit is due. The reference will be automatically printed when
        the class is instantiated. This is taken care of in the the ``super().__init__(model, name)`` 
        line of your __init__ method.  
        
        Example:

        .. code:: python
            self.reference = "1. Alegre‐Requena, J. V., Sowndarya S. V., S., Pérez‐Soto, R., Alturaifi, T. M. & Paton, R. S. AQME: Automated quantum mechanical environments for researchers and educators. WIREs Comput Mol Sci 13, e1663 (2023)."

        In some cases, there might not be a direct reference (see DistanceCutoff class), but 
        there might be relevant work a user might be interested in. Please only refer to the 
        work of interest in the class doc string, not in the reference method. 
        
        If there are no references, please only include the line:

        .. code:: python
            self.reference = None

        """
        ...


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

    def reference(self):
        self.reference = None
