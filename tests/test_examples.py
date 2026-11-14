"""The public example calculations in examples/ (run in CI; no private data)."""

from pathlib import Path

import numpy as np

import pedro

EX = Path(__file__).resolve().parents[1] / "examples"


def test_benzaldehyde_gvpt2():
    r = pedro.analyze(EX / "benzaldehyde_gvpt2_b3lyp_gaussian.log")
    a = r.anharmonic
    assert r.checks["reconstruction"]["max_abs_cm1"] < 0.001
    assert a["status"] == "ok" and a["counts"] == {"fundamental": 36, "overtone": 36, "combination": 630}
    assert a["resonance_counts"]["fermi"] == 12
    assert a["mapping"]["max_residual_cm1"] < 0.01
    g = a["symmetry"]
    assert g["group"] == "Cs"
    # Strongly mixed Darling-Dennison 1-1 states (51 %/36 %) must be labeled one-to-one, as Gaussian does;
    # labeling by the largest component alone gave a 78 cm-1 mismatch here.
    assert g["variational_parse_check_max_abs_cm-1"] < 0.05 and g["symmetry_filtered_available"]


def test_ethylene_gvpt2():
    a = pedro.analyze(EX / "ethylene_gvpt2_b3lyp_gaussian.log").anharmonic
    assert a["counts"] == {"fundamental": 12, "overtone": 12, "combination": 66}
    assert a["symmetry"]["group"] == "D2h" and a["symmetry"]["variational_parse_check_max_abs_cm-1"] < 0.05


def test_acetone_water_orca_embedded():
    r = pedro.analyze(EX / "acetone_2h2o_b3lyp_orca.hess", pedro.Options(anharmonic="off"))
    assert r.basis.case == "F" and len(r.basis.subsystem) == 10        # acetone described, waters as environment
    assert r.checks["reconstruction"]["max_abs_cm1"] < 0.001
    np.testing.assert_allclose(r.ted.sum(axis=1), 100, atol=1e-8)
