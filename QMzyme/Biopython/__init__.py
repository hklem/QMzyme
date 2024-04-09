# Portions of this code have been modified or adapted by 
# Heidi Klem for usage in the QMzyme package. 

# Copyright 1999-2003 by Jeffrey Chang.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""
Collection of modules for dealing with biological data in Python.

The Biopython Project is an international association of developers
of freely available Python tools for computational molecular biology.

https://biopython.org

* This code has been modified from the original Biopython code.
"""

import os
import warnings

__version__ = "1.84.dev0"


class MissingExternalDependencyError(Exception):
    """Missing an external dependency.

    Used for things like missing command line tools. Important for our unit
    tests to allow skipping tests with missing external dependencies.
    """


class MissingPythonDependencyError(MissingExternalDependencyError, ImportError):
    """Missing an external python dependency (subclass of ImportError).

    Used for missing Python modules (rather than just a typical ImportError).
    Important for our unit tests to allow skipping tests with missing external
    python dependencies, while also allowing the exception to be caught as an
    ImportError.
    """


class StreamModeError(ValueError):
    """Incorrect stream mode (text vs binary).

    This error should be raised when a stream (file or file-like object)
    argument is in text mode while the receiving function expects binary mode,
    or vice versa.
    """


class BiopythonWarning(Warning):
    """Biopython warning.

    Biopython should use this warning (or subclasses of it), making it easy to
    silence all our warning messages should you wish to:

    >>> import warnings
    >>> from Bio import BiopythonWarning
    >>> warnings.simplefilter('ignore', BiopythonWarning)

    Consult the warnings module documentation for more details.
    """


class BiopythonParserWarning(BiopythonWarning):
    """Biopython parser warning.

    Some in-valid data files cannot be parsed and will trigger an exception.
    Where a reasonable interpretation is possible, Biopython will issue this
    warning to indicate a potential problem. To silence these warnings, use:

    >>> import warnings
    >>> from Bio import BiopythonParserWarning
    >>> warnings.simplefilter('ignore', BiopythonParserWarning)

    Consult the warnings module documentation for more details.
    """


class BiopythonDeprecationWarning(BiopythonWarning):
    """Biopython deprecation warning.

    Biopython uses this warning instead of the built in DeprecationWarning
    since those are ignored by default since Python 2.7.

    To silence all our deprecation warning messages, use:

    >>> import warnings
    >>> from Bio import BiopythonDeprecationWarning
    >>> warnings.simplefilter('ignore', BiopythonDeprecationWarning)

    Code marked as deprecated is likely to be removed in a future version
    of Biopython. To avoid removal of this code, please contact the Biopython
    developers via the mailing list or GitHub.
    """


class BiopythonExperimentalWarning(BiopythonWarning):
    """Biopython experimental code warning.

    Biopython uses this warning for experimental code ('alpha' or 'beta'
    level code) which is released as part of the standard releases to mark
    sub-modules or functions for early adopters to test & give feedback.

    Code issuing this warning is likely to change (or even be removed) in
    a subsequent release of Biopython. Such code should NOT be used for
    production/stable code. It should only be used if:

    - You are running the latest release of Biopython, or ideally the
      latest code from our repository.
    - You are subscribed to the biopython-dev mailing list to provide
      feedback on this code, and to be alerted of changes to it.

    If all goes well, experimental code would be promoted to stable in
    a subsequent release, and this warning removed from it.
    """


_parent_dir = os.path.dirname(os.path.dirname(__file__))
if os.path.exists(os.path.join(_parent_dir, "setup.py")):
    # Looks like we are running from our source directory,
    # a bad idea except if installed in development mode.
    #
    # See https://setuptools.readthedocs.io/en/latest/userguide/development_mode.html
    # Do we have .../site-packages/biopython.egg-link present?
    #
    # Note "pip install -e ." currently calls setuptools internally
    import site

    _dev_mode = False
    for _p in site.getsitepackages():
        if os.path.isfile(os.path.join(_p, "biopython.egg-link")):
            _dev_mode = True
            break
    # Also check the user specific site packages
    if not _dev_mode and os.path.isfile(
        os.path.join(site.getusersitepackages(), "biopython.egg-link")
    ):
        _dev_mode = True
    if not _dev_mode:
        warnings.warn(
            "You may be importing Biopython from inside the source tree."
            " This is bad practice and might lead to downstream issues."
            " In particular, you might encounter ImportErrors due to"
            " missing compiled C extensions. We recommend that you"
            " try running your code from outside the source tree."
            " If you are outside the source tree then you have a"
            " setup.py file in an unexpected directory: " + _parent_dir,
            BiopythonWarning,
        )

# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Classes that deal with macromolecular crystal structures.

Includes: PDB and mmCIF parsers, a Structure class, a module to keep a local
copy of the PDB up-to-date, selective IO of PDB files, etc.

Original Author: Thomas Hamelryck.
Contributions by:
- Peter Cock
- Joe Greener
- Rob Miller
- Lenna X. Peterson
- Joao Rodrigues
- Kristian Rother
- Eric Talevich
- and many others.
"""

import numpy as np


# Get a Structure object from a PDB file
from .PDBParser import PDBParser

from .MMCIFParser import MMCIFParser
from .MMCIFParser import FastMMCIFParser
from .PDBMLParser import PDBMLParser

# Download from the PDB
from .PDBList import PDBList

# Parse PDB header directly
from .parse_pdb_header import parse_pdb_header

# Find connected polypeptides in a Structure
from .Polypeptide import PPBuilder, CaPPBuilder
from .Polypeptide import is_aa, standard_aa_names, is_nucleic

# IO of PDB files (including flexible selective output)
from .PDBIO import PDBIO, Select
from .MMCIFIO import MMCIFIO

# Some methods to eg. get a list of Residues
# from a list of Atoms.
from . import Selection

# Superimpose atom sets
from .Superimposer import Superimposer

# 3D vector class
from .vectors import Vector, calc_angle, calc_dihedral, refmat, rotmat, rotaxis
from .vectors import vector_to_axis, m2rotaxis, rotaxis2m

# Alignment module
from .StructureAlignment import StructureAlignment

# DSSP handle
# (secondary structure and solvent accessible area calculation)
from .DSSP import DSSP, make_dssp_dict

# Fast atom neighbor search
# Depends on kdtrees C module
try:
    from .NeighborSearch import NeighborSearch
except ImportError:
    pass
