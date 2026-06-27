# -*- coding: utf-8 -*-

"""This file contains metadata describing the results from the Atomic Charges step."""

metadata = {}

"""Description of the computational models for Atomic Charges.

Atomic charges are not a single well-defined quantity, so this step exposes the
partitioning schemes it can run. This mirrors the ``metadata["computational
models"]`` convention used by the QM steps, kept deliberately small here.
"""
metadata["computational models"] = {
    "Charge Partitioning": {
        "models": {
            "DDEC6": {
                "parameterizations": {
                    "DDEC6": {
                        "program": "Chargemol",
                        "periodic": True,
                        "reference": "Manz & Limas, RSC Adv. 6, 47771 (2016)",
                        "description": (
                            "Density Derived Electrostatic and Chemical charges; "
                            "ESP-faithful and chemically transferable."
                        ),
                    },
                },
            },
            "Bader": {
                "parameterizations": {
                    "Bader": {
                        "program": "bader",
                        "periodic": True,
                        "reference": "Henkelman, Arnaldsson & Jonsson, "
                        "Comput. Mater. Sci. 36, 354 (2006)",
                        "description": (
                            "Bader / QTAIM topological partition of the density "
                            "into zero-flux basins."
                        ),
                    },
                },
            },
        },
    },
}

"""Properties that the Atomic Charges step produces.

`metadata["results"]` describes the results that this step can produce. The keys are
the internal names used within this step; the values describe each result. See the
QM steps' metadata for the full set of recognized fields.
"""
metadata["results"] = {
    "atomic charges": {
        "description": "The partial charge on each atom",
        "dimensionality": ["n_atoms"],
        "type": "float",
        "units": "e",
    },
    "net charge": {
        "description": "Sum of the atomic charges (check against the system charge)",
        "dimensionality": "scalar",
        "type": "float",
        "units": "e",
    },
    "method": {
        "description": "The charge-partitioning method used",
        "dimensionality": "scalar",
        "type": "string",
    },
}
