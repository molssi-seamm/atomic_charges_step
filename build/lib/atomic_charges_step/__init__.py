# -*- coding: utf-8 -*-

"""
atomic_charges_step
A SEAMM plug-in for computing atomic (partial) charges from a converged
electron density (DDEC6 via Chargemol, Bader via the Henkelman code).
"""

# Bring up the classes so that they appear to be directly in
# the atomic_charges_step package.

from .atomic_charges import AtomicCharges  # noqa: F401, E501
from .atomic_charges_parameters import AtomicChargesParameters  # noqa: F401, E501
from .atomic_charges_step import AtomicChargesStep  # noqa: F401, E501
from .tk_atomic_charges import TkAtomicCharges  # noqa: F401, E501

from .metadata import metadata  # noqa: F401

# Handle versioneer
from ._version import get_versions

__author__ = "Paul Saxe"
__email__ = "psaxe@molssi.org"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
