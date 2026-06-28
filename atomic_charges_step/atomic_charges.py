# -*- coding: utf-8 -*-

"""Non-graphical part of the Atomic Charges step in a SEAMM flowchart.

This is a post-processing step: a preceding quantum-chemistry step produces an
electron density, and this step partitions that density into per-atom charges
using an external program (Chargemol for DDEC6, the Henkelman code for Bader).
The resulting charges are written back onto the configuration as a labeled charge
set, ``charges_<label>`` -- the same convention used by ``seamm_ff_util`` -- so
several schemes can coexist on one structure.

Status:

* Molecular **DDEC6** via Chargemol from a Gaussian ``.wfx`` is implemented and
  validated end-to-end (g09 -> wfx -> Chargemol -> parse -> store). ORCA feeds the
  same path: its Energy substep writes ``orca.wfx`` via orca_2aim, validated
  end-to-end (ORCA -> wfx -> Chargemol).
* The density handoff (``_locate_density``) reads a ``.wfx`` from the preceding
  step's directory (preferred), recognizing ``.cube`` as well. A uniform cube
  cannot represent the all-electron core, so molecular cube input is rejected;
  periodic CHGCAR/AECCAR support will come with VASP.
* The **Bader** backend raises a clear NotImplementedError: the Henkelman code
  reads a density grid (molecular cube or VASP CHGCAR/AECCAR), never a .wfx, and
  molecular all-electron cubes hit the same core-cusp limit as cube-based DDEC6.
  It will be implemented with periodic (VASP) support.
* Citations: the DDEC6 (Manz & Limas) and Bader (Henkelman group) methodology
  papers are cited per the method actually run.
"""

import csv
import importlib
import logging
import os
from pathlib import Path
import pprint  # noqa: F401
import textwrap

import numpy as np
from tabulate import tabulate

import atomic_charges_step
import molsystem  # noqa: F401
import seamm
import seamm_exec
from seamm_util import ureg, Q_  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("Atomic Charges")

# Add this module's properties to the standard properties
path = importlib.resources.files("atomic_charges_step") / "data"
csv_file = path / "properties.csv"
if path.exists() and csv_file.exists():
    molsystem.add_properties_from_file(csv_file)


class AtomicCharges(seamm.Node):
    """The non-graphical part of an Atomic Charges step in a flowchart.

    See Also
    --------
    TkAtomicCharges, AtomicCharges, AtomicChargesParameters
    """

    def __init__(
        self, flowchart=None, title="Atomic Charges", extension=None, logger=logger
    ):
        """A step for Atomic Charges in a SEAMM flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart that contains this step.
        title: str
            The name displayed in the flowchart.
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes
        """
        logger.debug(f"Creating Atomic Charges {self}")

        super().__init__(
            flowchart=flowchart,
            title="Atomic Charges",
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._metadata = atomic_charges_step.metadata
        self.parameters = atomic_charges_step.AtomicChargesParameters()

    @property
    def version(self):
        """The semantic version of this module."""
        return atomic_charges_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return atomic_charges_step.__git_revision__

    def create_parser(self):
        """Set up the command-line / seamm.ini parser for this step.

        Adds ``--max-atoms-to-print`` (default 20), overridable in the
        ``[atomic-charges-step]`` section of seamm.ini.
        """
        parser_name = self.step_type
        parser = self.flowchart.parser

        parser_exists = parser.exists(parser_name)

        result = super().create_parser(name=parser_name)

        if parser_exists:
            return result

        parser.add_argument(
            parser_name,
            "--max-atoms-to-print",
            default=20,
            help=(
                "Print the per-atom charges to the output for systems with no "
                "more than this many atoms. The charges are always written to a "
                "CSV file regardless."
            ),
        )

        return result

    def description_text(self, P=None):
        """Create the text description of what this step will do.

        Parameters
        ----------
        P: dict
            An optional dictionary of the current values of the control parameters.

        Returns
        -------
        str
            A description of the current step.
        """
        if not P:
            P = self.parameters.values_to_dict()

        method = P["method"]
        label = self._charge_label(P, method if method != "DDEC6 and Bader" else None)
        if method == "DDEC6 and Bader":
            text = (
                "Compute atomic charges using both DDEC6 (Chargemol) and Bader "
                "(Henkelman), storing them as separate labeled charge sets."
            )
        else:
            text = (
                f"Compute atomic charges using the {method} method, storing them "
                f"on the atoms as 'charges_{label}'."
            )

        return self.header + "\n" + __(text, indent=4 * " ").__str__()

    def run(self):
        """Run an Atomic Charges step.

        Returns
        -------
        seamm.Node
            The next node object in the flowchart.
        """
        next_node = super().run(printer)

        # Get the values of the parameters, dereferencing any variables
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # Make the working directory
        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)

        # Get the current system and configuration
        system, configuration = self.get_system_configuration(None)
        n_atoms = configuration.n_atoms

        # Echo what we are doing to the job output
        printer.important(__(self.description_text(P), indent=self.indent))

        # Find the density file(s) the upstream QM step produced
        density = self._locate_density(P)

        # Which methods to run
        if P["method"] == "DDEC6 and Bader":
            methods = ["DDEC6", "Bader"]
        else:
            methods = [P["method"]]

        results = {}
        for method in methods:
            if method == "DDEC6":
                charges = self._run_ddec6(P, density, configuration)
            elif method == "Bader":
                charges = self._run_bader(P, density, configuration)
            else:
                raise ValueError(f"Unknown charge method '{method}'")

            charges = np.asarray(charges, dtype=float)
            if charges.shape[0] != n_atoms:
                raise RuntimeError(
                    f"{method} returned {charges.shape[0]} charges for "
                    f"{n_atoms} atoms."
                )

            # Normalize to the known net charge, removing the small residual from
            # grid / numerical integration (a uniform shift preserves the relative
            # charges). The residual is reported so the shift is not hidden.
            target = float(configuration.charge)
            residual = float(charges.sum()) - target
            if P["enforce net charge"] == "yes" and n_atoms > 0:
                charges = charges - residual / n_atoms

            self._cite_method(method)

            label = self._charge_label(P, method)
            self._store_charges(configuration, label, charges)

            net = float(charges.sum())
            results[method] = {
                "charges": charges,
                "net charge": net,
                "charge residual": residual,
            }
            msg = (
                f"{method}: stored {n_atoms} charges as 'charges_{label}' "
                f"(net charge {net:.4f} e)."
            )
            if P["enforce net charge"] == "yes":
                msg += (
                    f" Corrected a residual of {residual:+.4f} e to match the "
                    f"system charge of {target:+.1f} e."
                )
            else:
                msg += f" Residual vs. system charge: {residual:+.4f} e."
            printer.normal(__(msg, indent=self.indent + 4 * " "))

        self._last_results = results

        # Analysis / saving of results, variables, tables
        self.analyze(P=P, results=results)

        return next_node

    # ------------------------------------------------------------------
    # Citations
    # ------------------------------------------------------------------
    # Per method: (bibtex key, citation level, note). The primary methodology
    # paper is level 1; supporting algorithm papers are level 2.
    _CITATIONS = {
        "DDEC6": [
            ("ddec6_part1", 1, "DDEC6 charge-partitioning theory and methodology."),
            (
                "ddec6_part2",
                2,
                "DDEC6 validation across periodic and molecular systems.",
            ),
        ],
        "Bader": [
            ("bader_tang2009", 1, "Grid-based Bader analysis without lattice bias."),
            (
                "bader_henkelman2006",
                2,
                "The original fast Bader decomposition algorithm.",
            ),
            ("bader_sanville2007", 2, "Improved grid-based Bader charge allocation."),
            (
                "bader_yu2011",
                2,
                "Weighted Bader integration (used with a reference density).",
            ),
        ],
    }

    def _cite_method(self, method):
        """Cite the methodology papers for the charge method actually run."""
        for alias, level, note in self._CITATIONS.get(method, []):
            if alias in self._bibliography:
                try:
                    self.references.cite(
                        raw=self._bibliography[alias],
                        alias=alias,
                        module="atomic_charges_step",
                        level=level,
                        note=note,
                    )
                except Exception as e:  # pragma: no cover
                    logger.warning(f"Could not cite {alias}: {e}")

    # ------------------------------------------------------------------
    # Charge storage
    # ------------------------------------------------------------------
    def _charge_label(self, P, method):
        """Resolve the charge-set label for a given method."""
        label = P["charge label"].strip()
        if label in ("", "<method>"):
            return method
        return label

    def _store_charges(self, configuration, label, charges):
        """Write a labeled charge set onto the atoms.

        Mirrors the ``charges_<label>`` pattern in ``seamm_ff_util/forcefield.py``:
        a configuration-dependent float attribute, one column per scheme, so DDEC6,
        Bader, etc. coexist rather than overwriting a single 'charge' column.
        """
        key = f"charges_{label}"
        atoms = configuration.atoms
        if key not in atoms:
            atoms.add_attribute(key, coltype="float", configuration_dependent=True)
        atoms[key][0:] = list(charges)
        logger.debug(f"Set column '{key}' to the charges")

    # ------------------------------------------------------------------
    # Density handoff (SCAFFOLD -- see module docstring)
    # ------------------------------------------------------------------
    def _locate_density(self, P):
        """Find the electron-density file to partition.

        Returns a dict describing the density, e.g.::

            {"format": "wfx", "path": Path(".../1.wfx")}

        For the molecular Gaussian path, the Gaussian step writes its wavefunction
        to ``<job>/<step_no>.wfx`` (mirroring its ``<step_no>.chk``), so the file
        from the immediately preceding step is ``<job>/<my_step - 1>.wfx``.

        (VASP CHGCAR/AECCAR and cube handling will extend this dict's "format".)
        """
        # Density files we know how to hand to the backend, best first: an
        # analytic Gaussian wavefunction (.wfx) is preferred over a grid cube.
        suffixes = {".wfx": "wfx", ".cube": "cube"}

        if P["density source"] == "from files":
            path = Path(P["density files"]).expanduser()
            if path.is_dir():
                found = [
                    p for p in sorted(path.glob("*")) if p.suffix.lower() in suffixes
                ]
                # wfx ahead of cube
                found.sort(key=lambda p: list(suffixes).index(p.suffix.lower()))
                if not found:
                    raise FileNotFoundError(
                        f"No .wfx or .cube density file found in {path}"
                    )
                path = found[0]
            fmt = suffixes.get(path.suffix.lower())
            if fmt is None:
                raise ValueError(f"Expected a .wfx or .cube density file, got {path}")
            return {"format": fmt, "path": path}

        # From the previous step in the flowchart. A molecular QM step writes its
        # wavefunction inside its own directory tree: the Gaussian step renames
        # its wavefunction to <gaussian_dir>/<substep_no>.wfx; the ORCA step
        # writes orca.wfx (via orca_2aim) in its Energy substep directory. Look in
        # the previous node's directory for a .wfx (or .cube), preferring .wfx.
        prev = self.previous()
        search_dirs = []
        if prev is not None and getattr(prev, "_id", None) is not None:
            search_dirs.append(Path(prev.directory))
        # Fall back to the job root in case a step wrote there directly.
        search_dirs.append(Path(self.directory).parent)

        for d in search_dirs:
            for suffix, fmt in suffixes.items():
                files = sorted(d.rglob(f"*{suffix}"), key=lambda p: p.stat().st_mtime)
                if files:
                    return {"format": fmt, "path": files[-1]}

        raise FileNotFoundError(
            "Could not find an electron-density file (.wfx or .cube) from the "
            f"preceding step (looked in {', '.join(str(d) for d in search_dirs)}). "
            "Enable 'Write a wavefunction (wfx) file' in a preceding Gaussian step, "
            "or 'Write the electron-density cube' in a preceding ORCA step, so this "
            "step has a density to partition."
        )

    # ------------------------------------------------------------------
    # DDEC6 / Chargemol backend
    # ------------------------------------------------------------------
    def _run_ddec6(self, P, density, configuration):
        """Run Chargemol to get DDEC6 charges. Returns a list of charges.

        Uses an analytic Gaussian wavefunction (.wfx). For molecular ORCA the
        preceding step writes a .wfx via orca_2aim, mirroring Gaussian.
        """
        if density["format"] == "cube":
            raise NotImplementedError(
                "DDEC6 from a density cube is not supported: a uniform grid "
                "cannot represent the all-electron core cusp, so Chargemol "
                "rejects molecular all-electron cubes. For molecular ORCA, "
                "enable 'Write the wavefunction (wfx) file' in the preceding "
                "ORCA step (it uses orca_2aim) so this step gets an analytic "
                ".wfx. (Periodic CHGCAR/AECCAR support will come with VASP.)"
            )
        if density["format"] != "wfx":
            raise NotImplementedError(
                f"The DDEC6 backend handles 'wfx' input; got '{density['format']}'."
            )

        atomic_densities = self._atomic_densities_path(P)
        src = Path(density["path"])

        # Chargemol reads the density file named in job_control.txt from the run
        # directory, so place a copy there.
        directory = Path(self.directory)
        (directory / src.name).write_bytes(src.read_bytes())

        # 3-periodic for a periodic configuration, molecular otherwise.
        periodic = getattr(configuration, "periodicity", 0) == 3
        job_control = self._chargemol_job_control(
            input_filename=src.name,
            atomic_densities=atomic_densities,
            periodicity=(periodic, periodic, periodic),
        )

        self._execute(
            section="chargemol",
            cmd=["{code}", ">", "stdout.txt", "2>", "stderr.txt"],
            files={"job_control.txt": job_control},
            return_files=[
                "stdout.txt",
                "stderr.txt",
                "DDEC6_even_tempered_net_atomic_charges.xyz",
            ],
        )
        charge_file = directory / "DDEC6_even_tempered_net_atomic_charges.xyz"
        return self._parse_ddec6_charges(charge_file, configuration.n_atoms)

    def _atomic_densities_path(self, P):
        """Resolve the atomic-densities directory, expanding ~ and env vars.

        Note: a leading ``$`` is avoided in the default because SEAMM treats a
        parameter value starting with ``$`` as a variable/expression to eval.
        """
        raw = os.path.expandvars(P["atomic densities directory"])
        return Path(raw).expanduser()

    def _chargemol_job_control(self, input_filename, atomic_densities, periodicity):
        """Build Chargemol's job_control.txt contents for a (molecular) wfx run.

        For a non-periodic wfx input the net charge and multiplicity are taken
        from the wfx itself, so the control file specifies periodicity, the
        reference-density directory, the input filename, and the charge type.
        """
        # Chargemol requires the densities path to end with a separator.
        densities = str(atomic_densities).rstrip("/") + "/"
        pa, pb, pc = (".true." if x else ".false." for x in periodicity)
        return (
            "<periodicity along A, B, and C vectors>\n"
            f"{pa}\n{pb}\n{pc}\n"
            "</periodicity along A, B, and C vectors>\n\n"
            "<atomic densities directory complete path>\n"
            f"{densities}\n"
            "</atomic densities directory complete path>\n\n"
            "<input filename>\n"
            f"{input_filename}\n"
            "</input filename>\n\n"
            "<charge type>\n"
            "DDEC6\n"
            "</charge type>\n\n"
            "<compute BOs>\n"
            ".false.\n"
            "</compute BOs>\n"
        )

    def _parse_ddec6_charges(self, charge_file, n_atoms):
        """Parse DDEC6 charges from DDEC6_even_tempered_net_atomic_charges.xyz.

        The file is XYZ-like: line 1 is the atom count, line 2 a comment, then
        one line per atom of 'Element x y z charge'; the net atomic charge is the
        fifth field. A Chargemol banner follows the atom block.
        """
        if not charge_file.exists():
            raise RuntimeError(
                f"Chargemol did not produce {charge_file.name}; check stdout.txt "
                "/ stderr.txt in the step directory."
            )
        lines = charge_file.read_text().splitlines()
        count = int(lines[0].split()[0])
        if count != n_atoms:
            raise RuntimeError(f"DDEC6 output has {count} atoms, expected {n_atoms}.")
        charges = []
        for line in lines[2 : 2 + count]:
            charges.append(float(line.split()[4]))
        return charges

    # ------------------------------------------------------------------
    # Bader / Henkelman backend
    # ------------------------------------------------------------------
    def _run_bader(self, P, density, configuration):
        """Run the Henkelman `bader` code. Returns a list of charges.

        Not yet implemented. The Henkelman code reads a density *grid* -- a
        Gaussian cube for molecules, or VASP CHGCAR + AECCAR0/AECCAR2 for periodic
        systems (``bader CHGCAR -ref CHGCAR_sum``) -- and never an analytic .wfx.
        A molecular all-electron cube suffers the same core-cusp sampling problem
        that rules out cube-based DDEC6, so Bader is the natural partner for
        periodic (VASP) densities and will be implemented with that support. The
        parse helper below sketches the ACF.dat handling for that work.
        """
        raise NotImplementedError(
            "Bader charges are not available yet. The Henkelman 'bader' code "
            "needs a density grid (a Gaussian cube for molecules, or VASP "
            "CHGCAR/AECCAR for periodic systems); it cannot use the analytic "
            ".wfx that the molecular path produces. Bader will be implemented "
            "with periodic (VASP) support -- for molecular charges now, use DDEC6."
        )

    def _parse_bader_charges(self, out, configuration):
        """Parse Bader charges from ACF.dat.

        SCAFFOLD: ACF.dat lists, per atom, the integrated electron count in the
        'CHARGE' column. The partial charge is ZVAL[element] - CHARGE, so this
        needs the per-element valence electron count from the upstream code.
        """
        raise NotImplementedError(
            "Bader output parsing not implemented yet (scaffold). Parse 'ACF.dat' "
            "and subtract from the per-element valence electron counts."
        )

    # ------------------------------------------------------------------
    # External-program execution (mirrors mopac_step's executor usage)
    # ------------------------------------------------------------------
    def _execute(self, section, cmd, files, return_files):
        """Run an external program through the flowchart executor.

        Parameters
        ----------
        section : str
            The section name in the program ini file (e.g. 'chargemol', 'bader'),
            also used as the executor configuration section.
        cmd : [str]
            The command template; '{code}' is filled in from the ini config.
        files : {str: str}
            Input files to write into the run directory.
        return_files : [str]
            Files to bring back from the run.

        Returns
        -------
        dict
            The executor result (stdout/stderr/return code and returned files).
        """
        executor = self.flowchart.executor
        config = self._program_config(section)

        seamm_exec.computational_environment()  # set resource limits

        result = executor.run(
            cmd=cmd,
            config=config,
            directory=self.directory,
            files=files,
            return_files=return_files,
            in_situ=True,
            shell=True,
        )
        if not result:
            raise RuntimeError(f"There was an error running {section}.")
        return result

    def _program_config(self, section):
        """Resolve the executable configuration for a program section.

        SCAFFOLD: this should read '<root>/atomic_charges.ini' (created from the
        default in data/ on first use) exactly as mopac_step reads 'mopac.ini',
        returning the dict for the active executor type. For now it locates the
        executable on PATH so the wiring can be exercised.
        """
        import shutil

        exe = {"chargemol": "Chargemol", "bader": "bader"}.get(section, section)
        path = shutil.which(exe) or shutil.which(section)
        if path is None:
            raise RuntimeError(
                f"Could not find the '{exe}' executable for the '{section}' "
                "method. Install it and/or register it in atomic_charges.ini."
            )
        return {"code": path, "installation": "local"}

    # ------------------------------------------------------------------
    def analyze(self, indent="", P=None, results=None, **kwargs):
        """Report the charge results: always write a CSV, and print a per-atom
        table for small systems.

        Parameters
        ----------
        indent: str
            An extra indentation for the output
        """
        if results is None:
            return

        _, configuration = self.get_system_configuration(None)
        symbols = list(configuration.atoms.symbols)
        n_atoms = configuration.n_atoms
        directory = Path(self.directory)

        # Threshold for printing the full table, from seamm.ini (default 20).
        try:
            max_print = int(self.options["max_atoms_to_print"])
        except (KeyError, TypeError, ValueError):
            max_print = 20

        for method, data in results.items():
            charges = data["charges"]

            # Save the results the user selected in the Results tab (variables,
            # tables, JSON). For several methods the generic keys hold the last
            # one; the per-method labeled charge sets on the structure
            # (charges_<label>) keep them distinct regardless.
            try:
                self.store_results(
                    configuration=configuration,
                    data={
                        "atomic charges": [float(q) for q in charges],
                        "net charge": float(data["net charge"]),
                        "charge residual": float(data["charge residual"]),
                        "method": method,
                    },
                )
            except Exception as e:  # pragma: no cover
                logger.warning(f"Could not store results for {method}: {e}")

            # Always write the charges to a CSV file.
            csv_path = directory / f"atomic_charges_{method}.csv"
            with open(csv_path, "w", newline="") as fd:
                writer = csv.writer(fd)
                writer.writerow(["Atom", "Element", "Charge"])
                for i, (symbol, q) in enumerate(zip(symbols, charges), start=1):
                    writer.writerow([i, symbol, f"{q:.6f}"])

            text = (
                f"{method} atomic charges (net charge {data['net charge']:.4f} e). "
                f"Written to {csv_path.name}."
            )
            printer.normal(__(text, indent=self.indent + 4 * " "))

            # Print the per-atom table for small systems.
            if n_atoms <= max_print:
                table = {
                    "Atom": [*range(1, n_atoms + 1)],
                    "Element": symbols,
                    "Charge": [f"{q:.4f}" for q in charges],
                }
                tmp = tabulate(
                    table,
                    headers="keys",
                    tablefmt="psql",
                    colalign=("center", "center", "decimal"),
                )
                printer.normal("")
                printer.normal(textwrap.indent(tmp, self.indent + 7 * " "))
                printer.normal("")
