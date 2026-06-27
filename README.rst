====================
Atomic Charges step
====================

A SEAMM plug-in for computing atomic (partial) charges from a converged electron
density.

This plug-in is a post-processing step: it runs after a quantum-chemistry step
(VASP, Gaussian, Psi4, ORCA, …) that has produced an electron density, reads that
density, and uses a real-space density-partitioning program to assign a partial
charge to each atom. The charges are written back into the configuration as a
labeled charge set (e.g. ``charges_DDEC6``) so several schemes can coexist on the
same structure.

Supported methods
-----------------

* **DDEC6** (via the Chargemol program) -- ESP-faithful, chemically transferable
  charges that work for molecular *and* periodic systems, including metals. The
  recommended general-purpose / force-field charge.
* **Bader / QTAIM** (via the Henkelman ``bader`` code) -- a rigorous topological
  partition of the density, the charge solid-state users expect.

Both consume a real-space density grid (or cube), so they apply equally to periodic
and molecular calculations. They need the *all-electron* density from the upstream
step (for VASP, ``CHGCAR`` plus ``AECCAR0``/``AECCAR2`` via ``LAECHG=.TRUE.``).

* Free software: BSD-3-Clause license
* Documentation: https://molssi-seamm.github.io/atomic_charges_step/index.html

Developed by the Molecular Sciences Software Institute (MolSSI_), which receives
funding from the `National Science Foundation`_.

.. _MolSSI: https://www.molssi.org
.. _`National Science Foundation`: https://www.nsf.gov
