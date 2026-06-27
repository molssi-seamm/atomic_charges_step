# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""Handle the installation of the Atomic Charges step."""

from .installer import Installer


def run():
    """Handle the extra installation needed.

    * Find and/or install the Chargemol (DDEC6) and bader executables.
    * Add or update information in the seamm.ini file for this step.
    """
    # Create an installer object
    installer = Installer()
    installer.run()


if __name__ == "__main__":
    run()
