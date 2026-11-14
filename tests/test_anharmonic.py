"""Resonance/overtone layer.

Public tests use synthetic data. Fixture tests reproduce the dev21 reference results (spec V8) and run only when
PEDRO_GVPT2_DATA points to a folder holding c2v_neutral_gvpt2.out and c2v_anion_gvpt2.log (a planar C2v molecule,
neutral and radical anion, wB97X-D/def2-SVP) and water_dimer_gvpt2.out (B3LYP/6-31+G(d)), all Gaussian 16 GVPT2 logs.
"""

import os
import time
from pathlib import Path

import numpy as np
import pytest

import pedro
from pedro import harmonic
from pedro.anharmonic import map_modes, quanta_text, renumber, run
from pedro.anharmonic import symmetry as sg
from pedro.anharmonic.screen import screen

from conftest import make_job


# ----------------------------------------------------------------------------- helpers
def test_quanta_text_and_renumber():
    assert quanta_text({33: 2}) == "2ν₃₃" and quanta_text({40: 1, 12: 1}) == "ν₁₂ + ν₄₀"
    mp = {1: 47, 5: 43, 6: 42, 14: 33}
    assert renumber("COUPLED_TO_FLAGGED(14,1)", mp) == "COUPLED_TO_FLAGGED(ν34,ν48)"
    assert renumber("PARENT_FLAGGED(5:ANOMALOUS_SHIFT)", mp) == "PARENT_FLAGGED(ν44:ANOMALOUS_SHIFT)"
    assert renumber("QUARTIC_OUTLIER(k(5,5,6,6)=+208)", mp) == "QUARTIC_OUTLIER(k(ν44,ν44,ν43,ν43)=+208)"
    assert renumber("FORBIDDEN_COUPLING(1-1:6(1)=+6.9)", mp) == "FORBIDDEN_COUPLING(1-1:ν43=+6.9)"
    assert renumber("RESONANCE_MIXED(57%)", mp) == "RESONANCE_MIXED(57%)"


def test_mapping_flags_near_degenerate_modes():
    freqs = np.array([500.0, 1000.0, 1000.05, 1500.0])
    mp, res, amb = map_modes({1: 1500.001, 2: 1000.049, 3: 1000.0, 4: 500.0}, freqs)
    assert mp == {1: 3, 2: 2, 3: 1, 4: 0} and max(res.values()) < 0.01
    assert amb == [2, 3]


# ----------------------------------------------------------------------------- point groups and speed (spec V9c)
def test_point_groups_synthetic():
    Z, m = [8, 1, 1], [15.995, 1.008, 1.008]
    xyz = np.array([[0, 0, 0.12], [0.76, 0, -0.47], [-0.76, 0, -0.47]])
    assert sg.point_group(Z, m, xyz)["group"] == "C2v"
    xyz2 = xyz.copy()
    xyz2[1, 0] = 0.80
    assert sg.point_group(Z, m, xyz2)["group"] == "Cs"
    Z4, m4 = [6, 1, 9, 17], [12.0, 1.008, 18.998, 34.969]
    xyz4 = np.array([[0, 0, 0], [0.6, 0.8, 0.3], [-0.9, 0.2, 0.7], [0.2, -1.1, -0.6]])
    assert sg.point_group(Z4, m4, xyz4)["group"] == "C1"


def test_symmetry_search_is_fast_for_large_models():
    w = np.array([[0, 0, 0.12], [0.76, 0, -0.47], [-0.76, 0, -0.47]])
    rng = np.random.default_rng(0)
    pts, Z = [], []
    for _ in range(115):
        s = rng.normal(size=3) * 8
        s[0] = abs(s[0]) + 1
        for sx in (1, -1):
            pts += list(w * [sx, 1, 1] + s * [sx, 1, 1])
            Z += [8, 1, 1]
    t = time.perf_counter()
    pg = sg.point_group(Z, [15.995 if z == 8 else 1.008 for z in Z], np.array(pts))
    assert pg["group"] == "Cs" and time.perf_counter() - t < 10
    t = time.perf_counter()
    assert sg.point_group(Z, [1.0] * len(Z), np.array(pts) + rng.normal(scale=0.2, size=(len(Z), 3)))["group"] == "C1"
    assert time.perf_counter() - t < 10


# ----------------------------------------------------------------------------- isolation rule (spec V9a, V9b)
def test_harmonic_job_gets_screen_and_no_gvpt2(embedded, tmp_path):
    job, coords = embedded
    res = pedro.analyze(job, pedro.Options(subsystem=[0, 1, 2], case="F", coords=None))
    assert res.anharmonic["status"] == "not_applicable"
    assert res.resonance_screen is not None and res.resonance_screen["status"] == "harmonic_screen"


def test_screen_finds_planted_coincidence(water):
    job, *_ = water
    m = harmonic.analyze(job.hessian, job.masses, job.xyz_bohr)
    f = m.freqs
    sc = screen(m, np.full(m.n, 100.0), ["a", "b", "c"], 3, delta=abs(f[2] - 2 * f[0]) + 1, min_ir=0)
    assert any(r["partner"] == "2ν1" and r["fundamental"] == 3 for r in sc["candidates"])
    assert all(abs(r["delta_cm1"]) <= sc["delta_cm1"] for r in sc["candidates"])


def _fake_gvpt2_log(path, e_harm):
    rows = "\n".join(f"      {k + 1}(1)                  {e:10.3f}  {e - 30:9.3f}     1.0     1.0" for k, e in enumerate(e_harm))
    path.write_text(" Normal modes will be PRINTED in DESCENDING order\n PT2 model: Deperturbed VPT2 (DVPT2)\n"
                    "     Anharmonic Infrared Spectroscopy\n\n Fundamental Bands\n -----------------\n"
                    "   Mode(n)                  E(harm)   E(anharm)        I(harm)       I(anharm)\n" + rows + "\n\n"
                    " Normal termination of Gaussian 16\n", encoding="latin-1")


def test_gvpt2_layer_fails_closed_on_imaginary_modes(water, tmp_path):
    job, coords, F = water
    H = job.hessian.copy()
    m0 = harmonic.analyze(H, job.masses, job.xyz_bohr)
    v = m0.Lq[:, 0] * np.repeat(np.sqrt(job.masses), 3)
    H2 = H - 2.5 * m0.eigenvalues[0] * np.outer(v, v)           # turn the softest mode imaginary
    log = tmp_path / "fake.log"
    _fake_gvpt2_log(log, [3000, 2000, 1000])
    job2 = make_job(job.elements, job.xyz_bohr, H2)
    job2.kind, job2.source, job2.anharmonic_present = "gaussian_log", str(log), True
    t = time.perf_counter()
    res = pedro.analyze(job2)
    assert time.perf_counter() - t < 10
    assert res.modes.n_imaginary == 1
    assert res.anharmonic["status"] == "failed" and "imaginary" in res.anharmonic["reason"]
    assert res.ted.shape == (3, 3) and res.checks["ted_closure_max_pp"] < 1e-9      # harmonic results intact


def test_gvpt2_layer_maps_and_labels_states(water, tmp_path):
    job, *_ = water
    m = harmonic.analyze(job.hessian, job.masses, job.xyz_bohr)
    log = tmp_path / "fake.log"
    _fake_gvpt2_log(log, sorted(m.freqs, reverse=True))
    job.kind, job.source, job.anharmonic_present = "gaussian_log", str(log), True
    a = run(job, m, ["low", "mid", "high"])
    assert a["status"] in ("ok", "incomplete")
    assert a["mapping"]["gaussian_to_pedro_1based"] == {"1": 3, "2": 2, "3": 1}
    s = {x["gaussian_quanta"]: x for x in a["states"]}
    assert s["1(1)"]["parentage"] == "[high]" and s["1(1)"]["pedro_quanta"] == "ν₃"


# ----------------------------------------------------------------------------- reference fixtures (spec V8)
@pytest.fixture
def gvpt2_dir():
    p = os.environ.get("PEDRO_GVPT2_DATA")
    if not p or not Path(p).exists():
        pytest.skip("set PEDRO_GVPT2_DATA to the GVPT2 reference fixture folder")
    return Path(p)


def test_c2v_neutral_reference(gvpt2_dir):
    a = pedro.analyze(gvpt2_dir / "c2v_neutral_gvpt2.out").anharmonic
    assert a["status"] == "ok" and a["counts"] == {"fundamental": 48, "overtone": 48, "combination": 1128}
    q = a["quality"]
    assert q["counts"] == {"ANOMALOUS_SHIFT": 14, "QUARTIC_OUTLIER": 16, "RESONANCE_MIXED": 19, "COUPLED_TO_FLAGGED": 13}
    assert q["modes"]["5"]["contamination_cm-1"] == pytest.approx(-5.477, abs=0.01)
    g = a["symmetry"]
    assert g["group"] == "C2v" and g["irrep_counts"] == {"A1": 17, "A2": 8, "B1": 7, "B2": 16}
    assert g["variational_parse_check_max_abs_cm-1"] < 0.05
    assert abs(g["modes"]["5"]["nu_symmetry_filtered"] - 1819.411) < 0.05
    assert abs(g["modes"]["6"]["nu_symmetry_filtered"] - 1821.399) < 0.05
    st = next(s for s in a["states"] if s["gaussian_quanta"] == "5(1)")
    assert st["irrep"] == "B2" and any(f.startswith("FORBIDDEN_COUPLING") for f in st["flags"])
    assert a["mapping"]["max_residual_cm1"] < 0.01


def test_c2v_anion_reference(gvpt2_dir):
    a = pedro.analyze(gvpt2_dir / "c2v_anion_gvpt2.log").anharmonic
    assert a["quality"]["counts"] == {"ANOMALOUS_SHIFT": 9, "QUARTIC_OUTLIER": 9, "RESONANCE_MIXED": 11,
                                      "COUPLED_TO_FLAGGED": 11}
    assert not any(c["class"] == "symmetry-forbidden" for c in a["symmetry"]["variational_couplings"])


def test_water_dimer_symmetry(gvpt2_dir):
    a = pedro.analyze(gvpt2_dir / "water_dimer_gvpt2.out").anharmonic
    assert a["symmetry"]["group"] == "Cs"
