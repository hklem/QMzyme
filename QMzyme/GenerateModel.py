###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@lsu.edu
###############################################################################

"""
**GenerateModel** is the primary user-facing module in QMzyme. :class:`~GenerateModel`
is used to load a starting structure, define QMzyme regions, and write calculation input.
To avoid unintended behavior, the initial structure must be pre-processed. I.e., ensure
hydrogens have been added, and the structure is representative of the system you hope
to study. If atomic charge information is not present in the input file(s), QMzyme
will guess atomic charges by referring to the residue names. Any residue name corresponding
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
from QMzyme.QMzymeRegion import QMzymeRegion
from QMzyme.utils import make_selection
from QMzyme.TruncationSchemes import TerminalAlphaCarbon
from QMzyme.CalculateModel import CalculateModel, CalculationFactory
from QMzyme.Writers import WriterFactory
from QMzyme.data import protein_residues

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
    :param pickle_file: Provide name/path+file of previously pickled QMzymeModel object to initialize
    :type pickle_file: str, default=None

    :Usage:

        To instantiate a model from a PDB file called "filename.pdb":

        .. code-block:: python

            model = QMzyme.GenerateModel("filename.pdb")

        If "filename.pdb" contains any components you know you do not want included in your model, you can initialize the
        GenerateModel instance from a selection of atoms by using the select_atoms argument:

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

        # A short module that allows for storing the selection parameters used to make region selection
        region._selection_attr.update(kwargs)
        if isinstance(selection, QMzymeRegion):
            region._selection_attr.setdefault('selection_scheme', selection._selection_attr.get('selection_scheme'))
        else:
            region._selection_attr['selection_scheme'] = selection
        
        self.add_region(region)
    

    def truncate(self, scheme = TerminalAlphaCarbon, name = None, remove_methane:bool = None, remove_ethane:bool = None, extend_gly_ala_backbone:bool = False, override_truncation:bool = None, override_capping:bool = None):
        """
        Method to truncate QMzymeModel. This method requires users to have methods assign to
        the region prior to using GenerateModel.truncate. All QMzymeModel regions with assigned
        methods will be combined and truncated according to the specified scheme. The resulting
        region will be saved as '{CalculateModel.calc_type}_region' if name = None.

        :param scheme: Specifies the truncation scheme to use. Options can be found
            in :py:mod:`~QMzyme.TruncationSchemes`.
        :type scheme: :py:class:`~QMzyme.TruncationSchemes.TruncationScheme` concrete class, 
            default=:class:`~QMzyme.TruncationSchemes.TerminalAlphaCarbon`
        :param name: Name to give the new truncated region. If None, the original
            region will be truncated.
        :type name: str, optional, default=None
        :param remove_methane: Controls how isolated Gly residues (which would
            otherwise be reduced to a methane upon truncation) are handled.
            If True, isolated Gly residues are removed from the region. If False,
            they are truncated and kept as a methane. If None and an isolated Gly is
            present (and extend_gly_ala_backbone is False), a ValueError is raised
            prompting the user to decide.
        :type remove_methane: bool, optional, default=None
        :param remove_ethane: Controls how isolated Ala residues (which would
            otherwise be reduced to a ethane upon truncation) are handled. If True,
            isolated Ala residues are removed from the region. If False, they are
            are truncated and kept as a ethane. If None and an isolated Ala is
            present (and extend_gly_ala_backbone is False), a ValueError is raised
            prompting the user to decide.
        :type remove_ethane: bool, optional, default=None
        :param extend_gly_ala_backbone: If True, isolated Gly/Ala residues are
            capped with ACE/NME groups instead of being removed or flagged,
            extending the backbone rather than truncating it down to a small
            organic group. Currently only supported with the
            :class:`~QMzyme.TruncationSchemes.TerminalAlphaCarbon` scheme.
        :type extend_gly_ala_backbone: bool, default=False
        :param override_truncation: Controls behavior for residues in the
            selection that have already been truncated. If None, a ValueError
            is raised prompting the user to decide. If False, already-truncated
            residues are skipped. If True, already-truncated residues are reverted
            to their untruncated state and re-truncated.
        :type override_truncation: bool, optional, default=None
        :param override_capping: Controls behavior for residues in the
            selection that have already been capped. If None, a ValueError
            is raised prompting the user to decide. If False, fully-capped
            residues are skipped; residues capped on only one terminus are
            still processed so the remaining uncapped terminus can be truncated.
            If True, already-capped residues are reverted to their uncapped
            state and re-capped/re-truncated.
        :type override_capping: bool, optional, default=None
        """
        # Combine regions
        if CalculateModel.calculation == {}:
            raise UserWarning("You must first assign calculation method(s) to the model region(s).")
        if len(CalculateModel.calculation) > 1:
            CalculateModel.combine_regions_and_methods()
        calc_type = CalculateModel.calc_type

        # Remember the original region name prior to making truncated region
        source_region = CalculateModel.calculation[calc_type]

        if name is None:
            name = f"{calc_type}_region"

        s = scheme(
            region=source_region, 
            name=name, 
            remove_methane=remove_methane, 
            remove_ethane=remove_ethane, 
            extend_gly_ala_backbone=extend_gly_ala_backbone,
            override_truncation=override_truncation,
            override_capping=override_capping
        )
        truncated_region = s.return_region()

        CalculateModel.calculation[calc_type] = truncated_region
        if calc_type != 'QM':
            CalculationFactory._make_calculation(calc_type)().assign_to_region(region=truncated_region)

        # Creates the truncated region as a whole region
        self.set_region(truncated_region)

        print(f"\nTruncated model has been created and saved as {name} "+
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
        if CalculateModel.calculation == {}:
            raise UserWarning("You must first assign calculation method(s) to the model region(s).")
        if len(CalculateModel.calculation) > 1 and not CalculateModel.combined:
            CalculateModel.combine_regions_and_methods()
            calc_type = CalculateModel.calc_type
            combined_region = CalculateModel.calculation[calc_type]
            combined_region.name = f"{calc_type}_region"
            CalculateModel.calculation[calc_type] = combined_region
            self.set_region(combined_region)

        calc_type = CalculateModel.calc_type
        region = CalculateModel.calculation.get(calc_type)
        if region is None or not getattr(region, "truncated", False):
            already_truncated = [res for res in region.residues if getattr(res, 'truncation_params', None) is not None]
            if not already_truncated:
                print("\nWARNING: model has not been truncated. Resulting model may "+
                    "not be a chemically complete structure (i.e., incomplete atomic "+
                    "valencies due to removed atoms).\n")
            else:
                protein_resid = [res for res in region.residues if res.resname in protein_residues]
                not_truncated = not_truncated = [res for res in protein_resid if getattr(res, 'truncation_params', None) is None]
                print("\nWARNING: model is only partially truncated. Resulting model may "+
                    "not be a chemically complete structure (i.e., incomplete atomic "+
                    "valencies due to removed atoms).\n"
                    f"Please truncate {not_truncated} using GenerateModel.truncate() or QMzymeRegion.truncate().")
        
        writer_type = CalculateModel.calc_type
        writer = WriterFactory.make_writer(writer_type, filename, memory, nprocs)
        if reset_calculation == True:
            CalculateModel._reset()
        self.store_pickle()
        #writer.write()
