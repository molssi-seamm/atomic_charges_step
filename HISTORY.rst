=======
History
=======

2026.6.28 -- Citations, charge normalization, and GUI fixes
    * The DDEC6 (Manz & Limas) and Bader (Henkelman) methodology papers are now
      cited for the charge method that is run.
    * The charges are normalized to the known net charge of the system with a
      small, reported uniform shift; this can be turned off.
    * Accepts an analytic wavefunction (.wfx) from a preceding ORCA step, in
      addition to Gaussian, for molecular DDEC6.
    * Fixed: the Results tab in the GUI was empty and selected results were not
      saved. The atomic charges, net charge, residual, and method are now listed
      and stored.
    * The charge method now offers Bader only when a density grid is available
      (a periodic/VASP density or explicit files); molecular densities offer
      DDEC6. Bader reports a clear message that it awaits periodic support.

2026.6.27 -- Initial release of the Atomic Charges step
    * Computes DDEC6 atomic charges (via Chargemol) from a molecular Gaussian
      wavefunction (.wfx), storing them as a labeled charge set on the structure.
