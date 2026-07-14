#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the `atomic_charges_step` package."""

import os
import shutil
import subprocess
from pathlib import Path

import pytest

import atomic_charges_step

DATA = Path(__file__).parent / "data"


def _seamm_chargemol():
    """Locate the seamm-chargemol conda env's chargemol + bundled densities.

    Returns ``(conda_exe, atomic_densities_dir)`` if the environment is present
    (so the end-to-end test can run), else ``None`` (so it skips -- e.g. in CI,
    where the installer has not run).
    """
    conda = os.environ.get("CONDA_EXE") or shutil.which("conda")
    if conda is None:
        return None
    prefix = Path(conda).resolve().parent.parent / "envs" / "seamm-chargemol"
    densities = prefix / "share" / "chargemol" / "atomic_densities"
    if (prefix / "bin" / "chargemol").is_file() and densities.is_dir():
        return conda, densities
    return None


_CHARGEMOL = _seamm_chargemol()


def test_factory_construction():
    """Create the stevedore helper class and check its type."""
    result = atomic_charges_step.AtomicChargesStep()
    assert (
        str(type(result))
        == "<class 'atomic_charges_step.atomic_charges_step.AtomicChargesStep'>"
    )


def test_node_construction():
    """Create the non-graphical node and check its type."""
    node = atomic_charges_step.AtomicCharges()
    assert (
        str(type(node)) == "<class 'atomic_charges_step.atomic_charges.AtomicCharges'>"
    )


def test_description():
    """The factory advertises itself in the 'Analysis' group."""
    desc = atomic_charges_step.AtomicChargesStep().description()
    assert desc["name"] == "Atomic Charges"
    assert desc["group"] == "Analysis"


def test_default_parameters():
    """Defaults are sensible: DDEC6, label tracks the method."""
    P = atomic_charges_step.AtomicChargesParameters()
    assert P["method"].value == "DDEC6"
    assert P["charge label"].value == "<method>"


def test_metadata_methods():
    """Both partitioning methods are described in the metadata."""
    models = atomic_charges_step.metadata["computational models"]
    families = models["Charge Partitioning"]["models"]
    assert "DDEC6" in families
    assert "Bader" in families


def test_charge_label_resolution():
    """'<method>' resolves to the method name; explicit labels pass through."""
    node = atomic_charges_step.AtomicCharges()
    assert node._charge_label({"charge label": "<method>"}, "DDEC6") == "DDEC6"
    assert node._charge_label({"charge label": "pbe0"}, "Bader") == "pbe0"


class _FakeAtoms(dict):
    """Minimal stand-in for molsystem atoms: dict of columns + add_attribute."""

    def __init__(self, n):
        super().__init__()
        self._n = n

    def add_attribute(self, key, coltype=None, configuration_dependent=None):
        self[key] = [0.0] * self._n


class _FakeConfig:
    def __init__(self, n):
        self.atoms = _FakeAtoms(n)


def test_apply_to_structure_default_off():
    """The new option defaults to off."""
    P = atomic_charges_step.AtomicChargesParameters()
    assert P["apply to structure"].value == "no"


def test_store_charges_labeled_only():
    """Without apply_to_structure only the labeled column is written."""
    node = atomic_charges_step.AtomicCharges()
    cfg = _FakeConfig(3)
    node._store_charges(cfg, "DDEC6", [-0.8, 0.4, 0.4])
    assert cfg.atoms["charges_DDEC6"] == [-0.8, 0.4, 0.4]
    assert "charge" not in cfg.atoms


def test_store_charges_apply_to_structure():
    """With apply_to_structure the standard 'charge' column is set too."""
    node = atomic_charges_step.AtomicCharges()
    cfg = _FakeConfig(3)
    node._store_charges(cfg, "DDEC6", [-0.8, 0.4, 0.4], apply_to_structure=True)
    assert cfg.atoms["charges_DDEC6"] == [-0.8, 0.4, 0.4]
    assert cfg.atoms["charge"] == [-0.8, 0.4, 0.4]


def test_has_ddec6_densities(tmp_path):
    """A directory counts as densities only if it holds c2_*.txt files."""
    node = atomic_charges_step.AtomicCharges
    assert node._has_ddec6_densities(tmp_path) is False  # empty
    (tmp_path / "c2_006_006_006_500_100.txt").write_text("x")
    assert node._has_ddec6_densities(tmp_path) is True
    assert node._has_ddec6_densities(tmp_path / "missing") is False


def test_conda_env_prefixes_absolute():
    """An absolute environment path is used directly (no derivation)."""
    node = atomic_charges_step.AtomicCharges()
    assert list(node._conda_env_prefixes("", "/abs/env")) == [Path("/abs/env")]


def test_conda_env_prefixes_guess_first():
    """For a named env the cheap <base>/envs/<name> guess comes first."""
    node = atomic_charges_step.AtomicCharges()
    prefixes = list(
        node._conda_env_prefixes("/opt/mc/condabin/conda", "seamm-chargemol")
    )
    assert prefixes[0] == Path("/opt/mc/envs/seamm-chargemol")


def test_chargemol_failure_message_quotes_log(tmp_path):
    """The failure message quotes Chargemol's real log, not the empty stdout."""
    (tmp_path / "orca.output").write_text(
        "...\nc2_006_006_006_500_100.txt\n"
        "Could not find a suitable reference density. Program will terminate.\n"
    )
    (tmp_path / "stderr.txt").write_text("Note: ... IEEE_UNDERFLOW_FLAG IEEE_DENORMAL")
    msg = atomic_charges_step.AtomicCharges._chargemol_failure_message(
        tmp_path, "orca.output"
    )
    assert "Could not find a suitable reference density" in msg
    assert "orca.output" in msg
    # The benign IEEE floating-point note is not surfaced as an error.
    assert "IEEE_" not in msg


def test_chargemol_failure_message_no_log(tmp_path):
    """When Chargemol wrote no log it is flagged as a likely startup failure."""
    msg = atomic_charges_step.AtomicCharges._chargemol_failure_message(
        tmp_path, "orca.output"
    )
    assert "failed to start" in msg


def test_chargemol_job_control():
    """The molecular wfx job_control has the right tags and a trailing slash."""
    node = atomic_charges_step.AtomicCharges()
    text = node._chargemol_job_control(
        input_filename="gaussian.wfx",
        atomic_densities="/some/where/atomic_densities",
        periodicity=(False, False, False),
    )
    assert "<input filename>\ngaussian.wfx\n</input filename>" in text
    assert "/some/where/atomic_densities/\n" in text  # slash appended
    assert text.count(".false.") >= 3  # periodicity
    assert "DDEC6" in text


def test_parse_ddec6_charges(tmp_path):
    """Parse the 5th column of a DDEC6 net-atomic-charges file."""
    f = tmp_path / "DDEC6_even_tempered_net_atomic_charges.xyz"
    f.write_text(
        "    3\n"
        "Nonperiodic system\n"
        "O      0.0   -0.0    0.117   -0.798055\n"
        "H      0.0    0.757  -0.469    0.399027\n"
        "H      0.0   -0.757  -0.469    0.399027\n"
        " \n Chargemol version 3.5 ...\n"
    )
    node = atomic_charges_step.AtomicCharges()
    charges = node._parse_ddec6_charges(f, 3)
    assert charges == pytest.approx([-0.798055, 0.399027, 0.399027])
    assert sum(charges) == pytest.approx(0.0, abs=1e-3)


def test_locate_density_from_files():
    """'from files' resolves a directory or an explicit .wfx path."""
    node = atomic_charges_step.AtomicCharges()
    P = {"density source": "from files", "density files": str(DATA)}
    density = node._locate_density(P)
    assert density["format"] == "wfx"
    assert density["path"].name == "water_b3lyp.wfx"


@pytest.mark.skipif(
    _CHARGEMOL is None,
    reason="the seamm-chargemol conda environment is not installed",
)
def test_ddec6_water_end_to_end(tmp_path):
    """Full external chain: water wfx -> Chargemol -> parse -> DDEC6 charges.

    Runs Chargemol in the seamm-chargemol conda environment (created by the
    installer), with the densities bundled in that environment. Uses the
    committed B3LYP/6-31G** water wfx fixture, so it is reproducible without
    Gaussian. Mirrors what AtomicCharges._run_ddec6 does.
    """
    conda, densities = _CHARGEMOL
    node = atomic_charges_step.AtomicCharges()
    shutil.copy(DATA / "water_b3lyp.wfx", tmp_path / "gaussian.wfx")
    (tmp_path / "job_control.txt").write_text(
        node._chargemol_job_control(
            input_filename="gaussian.wfx",
            atomic_densities=str(densities),
            periodicity=(False, False, False),
        )
    )
    subprocess.run(
        [conda, "run", "-n", "seamm-chargemol", "chargemol"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=300,
    )
    charge_file = tmp_path / "DDEC6_even_tempered_net_atomic_charges.xyz"
    charges = node._parse_ddec6_charges(charge_file, 3)

    assert sum(charges) == pytest.approx(0.0, abs=1e-3)  # neutral
    assert charges[0] == pytest.approx(-0.8, abs=0.1)  # O
    assert charges[1] == pytest.approx(0.4, abs=0.1)  # H
    assert charges[2] == pytest.approx(0.4, abs=0.1)  # H
