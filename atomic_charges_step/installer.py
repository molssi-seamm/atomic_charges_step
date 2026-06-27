# -*- coding: utf-8 -*-

"""Installer for the Atomic Charges plug-in.

This handles any further installation needed after installing the Python package
`atomic-charges-step`, namely locating (or installing) the external charge-analysis
programs -- Chargemol (for DDEC6) and the Henkelman `bader` code -- and registering
their locations in seamm.ini.
"""

import importlib
import logging
import subprocess

import seamm_installer

logger = logging.getLogger(__name__)


class Installer(seamm_installer.InstallerBase):
    """Handle further installation for atomic-charges-step.

    The Python package should already be installed with `pip` or `conda`. This
    plug-in-specific installer then checks for the Chargemol and `bader`
    executables, installing them if needed, and registers their locations in
    seamm.ini.
    """

    def __init__(self, logger=logger):
        # Call the base class initialization, which sets up the commandline
        # parser, amongst other things.
        super().__init__(logger=logger)

        logger.debug("Initializing the Atomic Charges installer object.")

        # Define this step's details
        self.environment = "seamm-atomic-charges"
        self.section = "atomic-charges-step"
        self.executables = ["Chargemol", "bader"]

        self.resource_path = importlib.resources.files("atomic_charges_step") / "data"

        # The environment.yaml file for Conda installations.
        logger.debug(f"data directory: {self.resource_path}")
        self.environment_file = self.resource_path / "seamm-atomic-charges.yml"

    def exe_version(self, config):
        """Get the version of the charge-analysis executables.

        Parameters
        ----------
        config : dict
            The configuration for this section from seamm.ini.

        Returns
        -------
        (str, str)
            The name and version reported by the executable, or 'unknown'.
        """
        # SCAFFOLD: neither Chargemol nor bader has a stable --version flag; the
        # version is usually printed in the program banner / output header. For now
        # just confirm the executable is reachable.
        code = config.get("code", "Chargemol")
        try:
            subprocess.run(
                f"'{code}' --help",
                stdin=subprocess.DEVNULL,
                capture_output=True,
                text=True,
                shell=True,
                timeout=30,
            )
        except Exception as e:
            logger.debug(f"    Could not query {code}: {e}")
            return "Chargemol/bader", "unknown"
        return "Chargemol/bader", "unknown"
