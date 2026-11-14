"""ORCA .hess input: a native ORCA 6.1.1 file (water, HF/def2-SVP, TightOpt Freq) in tests/data."""

from pathlib import Path

import numpy as np
import pytest

import pedro
from pedro.io import ReadError, orca

DATA = Path(__file__).parent / "data" / "water_hf_def2svp_orca611.hess"


def test_orca_hess_reads_and_reproduces_orca():
    job = pedro.load(DATA)
    assert job.kind == "orca_hess" and job.elements == ["O", "H", "H"]
    assert job.hessian.shape == (9, 9) and job.dipole_derivatives.shape == (9, 3)
    r = pedro.analyze(DATA, pedro.Options(anharmonic="off"))
    # ORCA prints frequencies to 0.01 cm-1 and intensities to 1e-8 km/mol
    assert np.max(np.abs(r.modes.freqs - job.printed_freqs)) < 0.005
    assert np.max(np.abs(r.modes.ir_km_mol - job.printed_ir)) < 1e-4
    np.testing.assert_allclose(r.ted.sum(axis=1), 100, atol=1e-8)
    assert [d.split()[0] for d in r.descriptions] == ["δ(H2–O1–H3)", "νₛ(O1–H2,H3)", "νₐₛ(O1–H2,H3)"]


def test_orca_hess_rejects_other_files(tmp_path):
    p = tmp_path / "x.hess"
    p.write_text("$hessian\n3\n")
    with pytest.raises(ReadError):
        orca.load(p)
