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
  validated end-to-end (g09 -> wfx -> Chargemol -> parse -> store).
* The density handoff (``_locate_density``) implements the Gaussian wfx case: the
  Gaussian step writes ``<job>/<step_no>.wfx`` and this step reads the preceding
  step's file. VASP CHGCAR/AECCAR and cube handling extend the same dict.
* The **Bader** backend is still a stub (needs VASP-style CHGCAR/AECCAR or a cube
  plus per-element valence-electron counts).
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

            label = self._charge_label(P, method)
            self._store_charges(configuration, label, charges)

            net = float(charges.sum())
            results[method] = {"charges": charges, "net charge": net}
            printer.normal(
                __(
                    f"{method}: stored {n_atoms} charges as 'charges_{label}' "
                    f"(net charge {net:.4f} e).",
                    indent=self.indent + 4 * " ",
                )
            )

        self._last_results = results

        # Analysis / saving of results, variables, tables
        self.analyze(P=P, results=results)

        return next_node

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
        if P["density source"] == "from files":
            path = Path(P["density files"]).expanduser()
            if path.is_dir():
                wfx = sorted(path.glob("*.wfx"))
                if not wfx:
                    raise FileNotFoundError(f"No .wfx file found in {path}")
                path = wfx[0]
            if path.suffix.lower() != ".wfx":
                raise ValueError(
                    f"Expected a .wfx file for the molecular DDEC6 path, got {path}"
                )
            return {"format": "wfx", "path": path}

        # From the previous step in the flowchart. The Gaussian step writes its
        # wfx inside its own directory (its Energy substep renames the file to
        # <gaussian_dir>/<substep_no>.wfx), so look in the previous node's
        # directory for a .wfx (most recent if several).
        prev = self.previous()
        search_dirs = []
        if prev is not None and getattr(prev, "_id", None) is not None:
            search_dirs.append(Path(prev.directory))
        # Fall back to the job root in case a step wrote there directly.
        search_dirs.append(Path(self.directory).parent)

        for d in search_dirs:
            wfx_files = sorted(d.rglob("*.wfx"), key=lambda p: p.stat().st_mtime)
            if wfx_files:
                return {"format": "wfx", "path": wfx_files[-1]}

        raise FileNotFoundError(
            "Could not find a wavefunction (.wfx) file from the preceding step "
            f"(looked in {', '.join(str(d) for d in search_dirs)}). Enable "
            "'Write a wavefunction (wfx) file' in the preceding Gaussian step so "
            "this step has a density to partition."
        )

    # ------------------------------------------------------------------
    # DDEC6 / Chargemol backend
    # ------------------------------------------------------------------
    def _run_ddec6(self, P, density, configuration):
        """Run Chargemol to get DDEC6 charges. Returns a list of charges."""
        if density["format"] != "wfx":
            raise NotImplementedError(
                f"The DDEC6 backend currently handles 'wfx' input; got "
                f"'{density['format']}'. (VASP/cube support is planned.)"
            )

        atomic_densities = self._atomic_densities_path(P)
        wfx = Path(density["path"])

        # Chargemol reads the density file named in job_control.txt from the run
        # directory, so place a copy there.
        directory = Path(self.directory)
        (directory / wfx.name).write_bytes(wfx.read_bytes())

        job_control = self._chargemol_job_control(
            input_filename=wfx.name,
            atomic_densities=atomic_densities,
            periodicity=(False, False, False),
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
        """Run the Henkelman `bader` code. Returns a list of charges."""
        # The bader code reads a grid file (e.g. CHGCAR) and, for VASP, the
        # reference all-electron density to reassign core charge:
        #   bader CHGCAR -ref CHGCAR_sum
        # where CHGCAR_sum = AECCAR0 + AECCAR2 (chgsum.pl). Net charge per atom is
        # Z_valence - ACF column; converting to a partial charge needs the number
        # of valence electrons per element (ZVAL), which is code/pseudopotential
        # dependent -- another reason the density handoff must carry metadata.
        out = self._execute(
            section="bader",
            cmd=[
                "{code}",
                "CHGCAR",
                "-ref",
                "CHGCAR_sum",
                ">",
                "stdout.txt",
                "2>",
                "stderr.txt",
            ],
            files={},  # SCAFFOLD: add located density files
            return_files=["stdout.txt", "stderr.txt", "ACF.dat"],
        )
        return self._parse_bader_charges(out, configuration)

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
