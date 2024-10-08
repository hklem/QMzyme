###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
**GenerateModel** is the primary user-facing module in QMzyme. :class:`~GenerateModel`
is used to load a starting structure, define QMzyme regions, and write calculation input.
To avoid unintended behavior, the initial stucture must be pre-processed. I.e., ensure
hydrogens have been added, and the structure is representative of the system you hope
to study. If atomic charge information is not present in the input file(s), QMzyme
will guess atomic charges by refering to the residue names. Any residue name corresponding
to standard protein residue names, defined `here <https://userguide.mdanalysis.org/1.1.1/standard_selections.html>`_, 
are able to be parsed for their total charge. This library can also be found in :py:mod:`~QMzyme.configuration` under
the dictionary protein_residues. If you have a non-protein residue QMzyme will assume its charge is 0. This
is important if you have a ligand with a non-zero charge that you will include in your calculation. After importing
QMzyme you can update the charge dictionary to add this information by adding to the residue_charges dictionary found 
in :py:mod:`~QMzyme.configuration`. 

The starting structure is loaded in using MDAnalysis and converted to a Universe object.
There are a variety of ways to define QMzyme region(s), and once a region has been set it
can be further modified, i.e., through truncation schemes. Lastly, this module interfaces with 
:class:`~QMzyme.CalculateModel.CalculateModel`, :class:`~QMzyme.Writers.Writer` and 
:class:`~QMzyme.aqme.qprep.qprep` to automate the creation of single- or multi-scale 
calculation input files. 
"""

from QMzyme.QMzymeModel import QMzymeModel
from QMzyme.utils import make_selection
from QMzyme.TruncationSchemes import TerminalAlphaCarbon
from QMzyme.CalculateModel import CalculateModel, CalculationFactory
from QMzyme.Writers import WriterFactory


class GenerateModel(QMzymeModel):
    """
    GenerateModel can be instantiated with an MDAnalysis Universe directly,
    or any combination of parameters that MDAnalysis.core.universe.Universe
    accepts to create a Universe.
    See https://userguide.mdanalysis.org/stable/universe.html for details.

    :param name: Name to give to the QMzymeModel. This is used for default file naming 
        purposes throughout the QMzyme package. If not provided, it will default to
        the base name of the universe filename attribute. 
    :type name: str, optional
    :param universe: MDAnalysis Universe object. If not specified, user will need to provide file(s) that
        MDAnalysis can use to create a Universe object.
    :type universe: `MDAnalysis.Universe <https://userguide.mdanalysis.org/stable/universe.html>`_, default=None
    :param select_atoms: `MDAnalysis selection command <https://docs.mdanalysis.org/stable/documentation_pages/selections.html>`_ 
        to specify which atoms are stored in the universe. 
    :type select_atoms: str, default='all'
    :param frame: If trajectory was provided, specify a frame to extract coordinates from.
    :type frame: int, default=0
    :param pickle_file: Provide name/path+file of previously pickled QMzymeModel object to inialize
    :type pickle_file: str, default=None

    :Usage:

        To instantiate a model from a PDB file called "filename.pdb":

        .. code-block:: python

            model = QMzyme.GenerateModel("filename.pdb")

        If "filename.pdb" contains any components you know you do not want included in your model, you can initialize the
        GenerateModel instance from a subselection of atoms by using the select_atoms argument:

        .. code-block:: python

            model = QMzyme.GenerateModel("filename.pdb", select_atoms="not resname WAT")

        You can also initialize the model from a topology and trajectory file(s) and specify what frame to take coordinates from:

        .. code-block:: python

            model = QMzyme.GenerateModel("filename.pdb", "filename.dcd", frame=100)

    """
    def __init__(self, *args, name=None, universe=None, select_atoms='all', frame=0, pickle_file=None, **kwargs):
        CalculateModel._reset()
        super().__init__(*args, name=name, universe=universe, frame=frame, select_atoms=select_atoms, pickle_file=pickle_file, **kwargs)

    def __repr__(self):
        return f"<QMzymeModel built from {self.universe} contains {self.n_regions} region(s)>"

    def set_catalytic_center(self, selection):
        """
        This method calls the set_region method and forces the region name 
        to be 'catalytic_center'. See :py:meth:`~QMzyme.GenerateModel.GenerateModel.set_region`.
        """
        self.set_region(selection=selection, name='catalytic_center')


    def set_region(self, selection, name=None, **kwargs):
        """
        Method to create a QMzymeRegion and add to the QMzymeModel regions list.

        :param selection: Determines what atoms are included in the region. A 
            variety of input options are accepted:

            * str that can be interpreted by `MDAnalysis selection commands <https://docs.mdanalysis.org/stable/documentation_pages/selections.html>`_

            * an `MDAnalysis AtomGroup`

            * a :class:`~QMzyme.QMzymeRegion.QMzymeRegion`

            * any concrete class of :class:`~QMzyme.SelectionSchemes.SelectionScheme`, i.e., :class:`~QMzyme.SelectionSchemes.DistanceCutoff`. Options can be found in :py:mod:`~QMzyme.SelectionSchemes`.
        
        :type selection: See options below, required
        :param name: Name of the resulting region.
        :type name: str, optional
        :param kwargs: Keyword arguments that might be needed if a :class:`~QMzyme.SelectionSchemes.SelectionScheme`
            is used. For example, the parameter `cutoff` is required to use the :class:`~QMzyme.SelectionSchemes.DistanceCutoff` 
            scheme. 

        :Usage:

            .. code-block:: python

                model.set_region(selection="resid 10 or resid 15", name="two_residues")

            .. code-block:: python

                from QMzyme.SelectionSchemes import DistanceCutoff
                model.set_region(selection=DistanceCutoff, cutoff=5)

        .. note::
            When using a :class:`~QMzyme.SelectionSchemes.SelectionScheme` the scheme class must be imported. 

        """
        region = make_selection(selection, model=self, name=name, **kwargs)
        self.add_region(region)
    

    def truncate(self, scheme=TerminalAlphaCarbon, name=None):
        """
        Method to truncate QMzymeModel. All QMzymeModel regions with assigned methods will be 
        combined and truncated according to the specified scheme. The resulting region will
        be saved as the QMzymeModel `truncated` attribute.

        :param scheme: Specifies the truncation scheme to use. Options can be found
            in :py:mod:`~QMzyme.TruncationSchemes`.
        :type scheme: :py:class:`~QMzyme.TruncationSchemes.TruncationScheme` concrete class, 
            default=:class:`~QMzyme.TruncationSchemes.TerminalAlphaCarbon`
        :param name: Name to give the truncated model. If None, the original
            region name will be the combination of calculation methods and the suffix '_combined_region_truncated'.
        :type name: str, optional
        """
        #combine regions
        if hasattr(self, "truncated"):
            raise UserWarning("Your model has already been truncated.")
        if CalculateModel.calculation == {}:
            raise UserWarning("You must first assign calculation method(s) to the model region(s).")
        if len(CalculateModel.calculation) > 1:
            CalculateModel.combine_regions_and_methods()
        calc_type = CalculateModel.calc_type
        s = scheme(region=CalculateModel.calculation[calc_type], name=name)
        region = s.return_region()
        if calc_type != 'QM':
            CalculationFactory._make_calculation(calc_type)().assign_to_region(region=region)
        CalculateModel.calculation[calc_type] = region
        setattr(self, "truncated", region)
        print(f"\nTruncated model has been created and saved to attribute 'truncated' "+
              "and stored in QMzyme.CalculateModel.calculation under key "+
              f"{calc_type}. This model will be used to write the calculation input.")

    def write_input(self, filename=None, memory='24GB', nprocs=12, reset_calculation=False):
        """
        Method to write calculation file input. The code will automatically
        detect what type of calculation file to prepare based on the 
        calculation methods that have been assigned to the model region(s). Once this is called
        the QMzymeModel object will automatically be serialized using the pickle library and saved
        under the filename {self.name+'.pkl'} in the current working directory.

        :param filename: Name to use for resulting file. If not specified, the 
            file will be named according to the region(s) name. The file format ending
            does not need to be included.
        :type filename: str, optional
        :param memory: Amount of memory to specify in the input file. 
        :type memory: str, optional
        :param nprocs: Number of processors to specify in the input file.
        :type nprocs: int, optional

        .. note::

            A :class:`~QMzyme.CalculateModel.QM_Method` must first be assigned
            to a region. 
        """
        if not hasattr(self, "truncated"):
            print("\nWARNING: model has not been truncated. Resulting model may "+
                  "not be a chemically complete structure (i.e., incomplete atomic "+
                  "valencies due to removed atoms).\n")
            CalculateModel.combine_regions_and_methods()
        
        writer_type = CalculateModel.calc_type
        writer = WriterFactory.make_writer(writer_type, filename, memory, nprocs)
        if reset_calculation == True:
            CalculateModel._reset()
        self.store_pickle()
        #writer.write()
