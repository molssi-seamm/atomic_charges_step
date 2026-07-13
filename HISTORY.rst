=======
History
=======

2026.7.13.1 -- Bugfix: clearer Chargemol failures; skip incomplete densities
    * When Chargemol fails to produce charges, the error now quotes its real log
      (``<input>.output``) instead of the empty ``stdout.txt`` -- so causes like
      "Could not find a suitable reference density" are shown directly. Chargemol
      writes its diagnostics to that log, not to stdout.
    * The atomic-densities directory is now used only if it actually contains the
      DDEC6 reference densities (``c2_*.txt``); an empty or incomplete directory
      is skipped in favor of the complete set bundled in the seamm-chargemol conda
      environment. This avoids a silent Chargemol failure when the configured
      directory exists but is not populated.

2026.7.13 -- Option to set the charges on the structure
    * New **Set as the atomic charges on the structure** option (default off).
      When on, the computed charges are also written to the structure's standard
      per-atom ``charge`` attribute (in addition to the labeled ``charges_<label>``
      column), so they travel with the structure -- e.g. to write them to an
      extended-XYZ (extxyz) file for machine-learning training.

2026.6.29 -- Run Chargemol from a dedicated conda environment
    * The installer now creates a 'seamm-chargemol' conda environment containing
      Chargemol (for DDEC6), and the step runs Chargemol in that environment --
      no hand-built installation or PATH setup is needed. Run
      'atomic-charges-step-installer install' to set it up.
    * The reference atomic densities are found automatically inside the
      seamm-chargemol environment; the 'DDEC reference densities' setting is now
      needed only to point at a different copy.

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
