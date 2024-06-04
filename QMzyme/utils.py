###############################################################################
# Code written by Heidi Klem.
# e: heidiklem@yahoo.com or heidi.klem@nist.gov
###############################################################################

"""
Module containing broad functionality utilized throughout the package.
"""

def check_filename(filename, format):
    if filename.endswith(format):
        return filename
    if not format.startswith('.'):
            format = '.'+format
    return filename.split('.')[0]+format

from functools import singledispatch
from QMzyme.QMzymeModel import QMzymeModel
from QMzyme.RegionBuilder import RegionBuilder
import QMzyme.MDAnalysisWrapper as MDAwrapper
from MDAnalysis.core.groups import AtomGroup
from QMzyme.SelectionSchemes import SelectionScheme
from QMzyme.SelectionSchemes import DistanceCutoff
from abc import ABCMeta
from QMzyme.QMzymeRegion import QMzymeRegion


@singledispatch
def make_selection(selection, model: QMzymeModel, name=None, **kwargs):
    """
    Method to enable variable input comability: will return an MDA AtomGroup if 
    input was an MDA selection command str, or return the input if it was either 
    an MDA AtomGroup or QMzymeRegion.
    """
    raise UserWarning(f"Invalid selection {selection}.")
    #print('make selection from: ', selection)
    #return selection

@make_selection.register
def MDA_str_selection(selection: str, model: QMzymeModel, name, **kwargs):
    region_builder = RegionBuilder(name)
    selection = MDAwrapper.select_atoms(model.universe, selection)
    region_builder.init_atom_group(selection)
    region = region_builder.get_region()
    return region

@make_selection.register
def MDA_AtomGroup_selection(selection: AtomGroup, model: QMzymeModel, name, **kwargs):
    region_builder = RegionBuilder(name)
    region_builder.init_atom_group(selection)
    region = region_builder.get_region()
    return region

@make_selection.register
def MDA_AtomGroup_selection(selection: QMzymeRegion, model: QMzymeModel, name, **kwargs):
    if name is not None:
         selection.name = name
    return selection

@make_selection.register
def scheme_selection(selection: ABCMeta, model: QMzymeModel, name, **kwargs):
    s = selection(model=model, name=name, **kwargs)
    region = s.return_region()
    return region