.. _user-guide:

**********
User Guide
**********

The Atomic Charges step partitions an electron density from a preceding quantum
chemistry step into atomic charges: **DDEC6** with Chargemol, from an analytic
wavefunction file (``.wfx``), and **Bader** charges. A molecular Gaussian or ORCA step
writes the ``.wfx`` when asked to ("Write the wavefunction (wfx) file"); the step finds
it in the previous step's directory, or reads density files given to it.

The charges are stored on the atoms as a labeled charge set, ``charges_<label>``, so
several schemes coexist on one structure, and optionally also as the structure's
atomic charges ("Set as the atomic charges on the structure"), which then travel with
it, e.g. into an extended-XYZ file. By default the small residual from the numerical
integration is removed so that the charges sum to the system's charge, and the
residual is reported.

Which density
=============

The charges are those of the density in the file. For a double hybrid, or MP2, ORCA
can write either the SCF density -- for a double hybrid, that of its hybrid-DFT part
-- or the relaxed MP2 density, the method's own (ORCA's "Density in the wfx file").
ORCA records which beside the file (``wfx_density.json``), and this step then says
which density it partitioned, and stores the charges of a relaxed MP2 density under
their own default label, ``charges_DDEC6_MP2relaxed``, so that the two kinds of
charges are never mixed silently. A label given in "Store charges as" is used as
given.

..
    <remove the dots above and this line and unindent the toctree to expose it>
    Contents:

    .. toctree::
       :glob:
       :maxdepth: 2
       :titlesonly:

       *

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
