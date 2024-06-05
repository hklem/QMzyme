###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
**GeneratModel** is the primary user-facing module in QMzyme. :class:`~GenerateModel`
is used to load a starting structure, define QMzyme regions, and write calculation input.
The starting structure is loaded in using MDAnalysis, converting it to a Universe object.
There are a variety of ways to define the QMzyme regions, and once a region has been set it
can be further modified through truncation schemes. Lastly, this module interfaces with 
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
    accepts to create a Universe i.e., (example.prmtop, example.dcd, dt=5).
    See https://userguide.mdanalysis.org/stable/universe.html for details.

    :param name: Name of QMzymeModel.
    :type name: str, default=None
    :param universe: MDAnalysis Universe object.
    :type universe: `MDAnalysis.Universe <https://userguide.mdanalysis.org/stable/universe.html>`_, default=None
    :param frame: If trajectory was provided, specify a frame to base coordinates on
    :type frame: int, default=0
    """
    def __init__(self, *args, name=None, universe=None, frame=0, **kwargs):
        CalculateModel._reset()
        super().__init__(*args, name=name, universe=universe, frame=frame, **kwargs)

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

            * an :class:`~MDAnalysis.core.groups.AtomGroup`

            * a :class:`~QMzyme.QMzymeRegion.QMzymeRegion`

            * any concrete class of :class:`~QMzyme.SelectionSchemes.SelectionScheme`, i.e., :class:`~QMzyme.SelectionSchemes.DistanceCutoff` (note, the class must be imported)
        
        :type selection: Any, see description above, reqiured
        :param name: Name of the resulting region.
        :type name: str, optional
        :param kwargs: Other parameters might need to be passed if a :class:`~QMzyme.SelectionSchemes.SelectionScheme`
            is used. For example, the parameter `cutoff` is required to use the :class:`~QMzyme.SelectionSchemes.DistanceCutoff` 
            scheme. 

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
            default=:py:class:`~QMzyme.TruncationSchemes.CA_terminal`
        :param name: Name to give the truncated model to. If None, the original
            region name will be the combination of calculation methods '_combined_region_truncated' appended.
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

    def write_input(self, filename=None, memory='24GB', nprocs=12):
        """
        Method to write calculation file input. The code will automatically
        detect what type of calculation file to prepare based on the 
        calculation methods that have been assigned to the model region(s). 

        :param filename: Name to use for resulting file. If not specified, the 
            file will be named according to the region(s) name. The file format ending
            does not need to be included.
        :type filename: str, optional
        :param memory: Amount of memory to specify in the input file. 
        :type memory: str, optional
        :param nprocs: Number of processors to specify in the input file.
        :type nprocs: int, optional

        :notes:
            A :class:`~QMzyme.CalculateModel.QM_Method` must have been assigned
            to a region. 
        """
        if not hasattr(self, "truncated"):
            print("\nWARNING: model has not been truncated. Resulting model may "+
                  "not be a chemically complete structure (i.e., incomplete atomic "+
                  "valencies due to removed atoms).\n")
            CalculateModel.combine_regions_and_methods()
        
        writer_type = CalculateModel.calc_type
        writer = WriterFactory.make_writer(writer_type, filename, memory, nprocs)
        #writer.write()
