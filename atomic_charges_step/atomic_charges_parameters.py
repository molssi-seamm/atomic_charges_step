# -*- coding: utf-8 -*-
"""
Control parameters for the Atomic Charges step in a SEAMM flowchart
"""

import logging
import seamm

logger = logging.getLogger(__name__)


class AtomicChargesParameters(seamm.Parameters):
    """The control parameters for the Atomic Charges step.

    Parameters
    ----------
    parameters : {str: {str: str}}
        A dictionary describing each control parameter, keyed by the parameter
        name. See ``seamm.Parameters`` for the full description of each field
        ("default", "kind", "default_units", "enumeration", "format_string",
        "description", "help_text").

    See Also
    --------
    AtomicCharges, TkAtomicCharges
    """

    parameters = {
        "method": {
            "default": "DDEC6",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "DDEC6",
                "Bader",
                "DDEC6 and Bader",
            ),
            "format_string": "",
            "description": "Charge method:",
            "help_text": (
                "The density-partitioning scheme used to assign atomic charges. "
                "DDEC6 (via Chargemol) gives ESP-faithful, transferable charges; "
                "Bader (via the Henkelman code) gives a rigorous topological "
                "partition. Both work for molecular and periodic systems."
            ),
        },
        "charge label": {
            "default": "<method>",
            "kind": "string",
            "default_units": "",
            "enumeration": ("<method>",),
            "format_string": "",
            "description": "Store charges as:",
            "help_text": (
                "Label for the charge set written onto the atoms. The charges are "
                "stored in the configuration as the attribute 'charges_<label>', "
                "so several schemes can coexist on the same structure. The default "
                "'<method>' uses the method name, e.g. 'charges_DDEC6'."
            ),
        },
        "density source": {
            "default": "from the previous step in the flowchart",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "from the previous step in the flowchart",
                "from files",
            ),
            "format_string": "",
            "description": "Electron density:",
            "help_text": (
                "Where to obtain the electron density. Normally it is taken from "
                "the quantum-chemistry step immediately preceding this one, which "
                "must have been told to write its (all-electron) density."
            ),
        },
        "density files": {
            "default": "",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Density file(s):",
            "help_text": (
                "Explicit path(s) to the density grid/cube file(s) when not taking "
                "them from the previous step. For VASP this is the directory holding "
                "CHGCAR plus AECCAR0/AECCAR2; for molecular codes a .cube or .wfx."
            ),
        },
        "atomic densities directory": {
            "default": "~/SEAMM/atomic_charges/atomic_densities",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "DDEC reference densities:",
            "help_text": (
                "Directory containing the reference atomic densities required by "
                "Chargemol for DDEC6. If this path does not exist, the copy "
                "bundled in the seamm-chargemol conda environment is used "
                "automatically. Ignored for the Bader method."
            ),
        },
        "enforce net charge": {
            "default": "yes",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Enforce the total charge:",
            "help_text": (
                "Shift the charges uniformly so they sum exactly to the known net "
                "charge of the system, removing the small residual left by grid / "
                "numerical integration. The residual is reported. Turn off to keep "
                "the raw charges from the partitioning program."
            ),
        },
        "results": {
            "default": {},
            "kind": "dictionary",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "results",
            "help_text": "The results to save to variables or in tables.",
        },
    }

    def __init__(self, defaults={}, data=None):
        """Initialize the parameters, by default with the parameters defined above

        Parameters
        ----------
        defaults: dict
            A dictionary of parameters to initialize. The parameters above are
            used first and any given will override/add to them.
        data: dict
            A dictionary of keys and a subdictionary with value and units for
            updating the current, default values.
        """
        logger.debug("AtomicChargesParameters.__init__")

        super().__init__(
            defaults={**AtomicChargesParameters.parameters, **defaults}, data=data
        )
