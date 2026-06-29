# -*- coding: utf-8 -*-

"""Installer for the Atomic Charges plug-in.

This handles any further installation needed after installing the Python package
`atomic-charges-step`, namely providing Chargemol (for DDEC6) -- by default in a
dedicated `seamm-chargemol` conda environment -- and registering it in seamm.ini
so the step runs Chargemol there.
"""

import importlib
import logging
from pathlib import Path
import re
import subprocess

import seamm_installer

logger = logging.getLogger(__name__)


class Installer(seamm_installer.InstallerBase):
    """Handle further installation for atomic-charges-step.

    The Python package should already be installed with `pip` or `conda`. This
    plug-in-specific installer then provides Chargemol (for DDEC6), by default in
    a dedicated `seamm-chargemol` conda environment, and registers it in
    seamm.ini so the step runs Chargemol in that environment.

    (The Henkelman `bader` code is not yet wired up -- Bader charges await
    periodic/VASP support -- so it is not installed here.)
    """

    def __init__(self, logger=logger):
        # Call the base class initialization, which sets up the commandline
        # parser, amongst other things.
        super().__init__(logger=logger)

        logger.debug("Initializing the Atomic Charges installer object.")

        # Define this step's details
        self.environment = "seamm-chargemol"
        self.section = "atomic-charges-step"
        self.executables = ["chargemol"]

        self.resource_path = importlib.resources.files("atomic_charges_step") / "data"

        # The environment.yaml file for Conda installations.
        logger.debug(f"data directory: {self.resource_path}")
        self.environment_file = self.resource_path / "seamm-chargemol.yml"

    def exe_version(self, config):
        """Get the version of the Chargemol executable.

        Parameters
        ----------
        config : dict
            The configuration for this section (installation, conda,
            conda-environment, code).

        Returns
        -------
        (str, str)
            ("Chargemol", <version>). Chargemol has no --version flag; run with
            no input it prints "Starting Chargemol version <X>" before exiting,
            which we parse (else 'unknown').
        """
        code = config.get("code", "chargemol")
        command = f"{code}"
        if config.get("installation") == "conda":
            conda = config["conda"]
            environment = config["conda-environment"]
            if environment[0] == "~" or Path(environment).is_absolute():
                environment = str(Path(environment).expanduser())
                flag = "-p"
            else:
                flag = "-n"
            command = f"'{conda}' run --live-stream {flag} '{environment}' {code}"

        version = "unknown"
        try:
            result = subprocess.run(
                command,
                stdin=subprocess.DEVNULL,
                capture_output=True,
                text=True,
                shell=True,
                timeout=60,
            )
        except Exception as e:
            logger.debug(f"    Could not query {code}: {e}")
        else:
            text = (result.stdout or "") + (result.stderr or "")
            match = re.search(r"Chargemol version\s+(\S+)", text)
            if match:
                version = match.group(1)
        return "Chargemol", version
