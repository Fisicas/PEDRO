"""VEDA file interchange: .fmu export, .fmv import, .dd2 export in VEDA's layout, and (optionally) VEDA's own
example files, which are read from PEDRO_VEDA_EXAMPLES and never copied into this repository."""

import os
from pathlib import Path

import numpy as np
import pytest

import pedro
from pedro.coordinates import export_dd2, import_dd2
from pedro.io import veda

ORCA = Path(__file__).parent / "data" / "water_hf_def2svp_orca611.hess"


def test_fmu_roundtrip(tmp_path):
    job = pedro.load(ORCA)
    p = veda.write_fmu(job, tmp_path / "w.fmu")
    back = pedro.load(p)
    assert back.elements == job.elements
    np.testing.assert_allclose(back.hessian, job.hessian, atol=1e-7)
    np.testing.assert_allclose(back.xyz_bohr, job.xyz_bohr, atol=1e-7)
    with pytest.raises(FileExistsError):
        veda.write_fmu(job, p)


def test_fmv_reads_veda_layout(tmp_path):
    job = pedro.load(ORCA)
    n = job.natoms
    lines = ["Atomic orientation and masses"]
    for i, (e, x, m) in enumerate(zip(job.elements, job.xyz_bohr * 0.529177210903, job.masses), 1):
        lines.append(f"{i:3d} {['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O'].index(e) + 1:3d} {x[0]:14.8f} {x[1]:14.8f} {x[2]:14.8f} {m:14.8f}")
    lines += ["", "Force constants in Cartesian coordinates"]
    v = (job.hessian * veda.FMV_MDYN_PER_ANG_PER_AU).ravel()
    lines += [" ".join(f"{x:15.8E}" for x in v[a:a + 5]) for a in range(0, len(v), 5)]
    p = tmp_path / "w.fmv"
    p.write_text("\n".join(lines) + "\n")
    back = pedro.load(p)
    assert back.kind == "veda_fmv"
    np.testing.assert_allclose(back.hessian, job.hessian, atol=1e-8)
    r = pedro.analyze(p, pedro.Options(anharmonic="off"))
    np.testing.assert_allclose(r.modes.freqs, job.printed_freqs, atol=0.01)


def test_dd2_export_matches_veda_layout(tmp_path):
    job = pedro.load(ORCA)
    r = pedro.analyze(ORCA, pedro.Options(anharmonic="off"))
    text = export_dd2(r.basis.internal, job.elements, job.xyz_bohr, tmp_path / "w.dd2")
    lines = text.splitlines()
    assert lines[0][:15] == "s 1     1.00   " and lines[0][15:19] in ("STRE", "BEND", "TORS", "OUT ")   # VEDA's columns
    raw = (tmp_path / "w.dd2").read_bytes()
    assert b"\r\n" in raw and b"\n\n" not in raw
    back = import_dd2(tmp_path / "w.dd2", job.elements)

    def direction(cs):     # coefficients up to the scale of each coordinate
        out = []
        for c in cs:
            s = max(abs(w) for _, w in c.components)
            out.append([(p.kind, p.atoms, round(w / s, 2)) for p, w in c.components])
        return out
    assert direction(back) == direction(r.basis.internal)
    r2 = pedro.analyze(ORCA, pedro.Options(coords=str(tmp_path / "w.dd2"), anharmonic="off"))
    np.testing.assert_allclose(r2.ted, r.ted, atol=1e-6)


@pytest.fixture
def veda_examples():
    p = os.environ.get("PEDRO_VEDA_EXAMPLES")
    if not p or not Path(p).exists():
        pytest.skip("set PEDRO_VEDA_EXAMPLES to VEDA 4's folder holding examp_1/ and examp_2/")
    return Path(p)


@pytest.mark.parametrize("folder,inp,stem,ted_within_1pp", [
    ("examp_1/examp_1", "a1-91.fmv", "a1-91", 0.995), ("examp_1/examp_1", "a1-91.fmv", "a1-91a", 0.995),
    ("examp_2/examp_2", "neopent.fmu", "neopent", 0.98), ("examp_2/examp_2", "neopent.fmu", "neopent1", 0.98),
    ("examp_2/examp_2", "neopent.fmu", "neopent2", 0.98), ("examp_2/examp_2", "neopent.fmu", "neopent3", 0.99)])
def test_veda_own_examples(veda_examples, tmp_path, folder, inp, stem, ted_within_1pp):
    src = veda_examples / folder
    p = tmp_path / (stem + Path(inp).suffix)
    p.write_bytes((src / inp).read_bytes())
    for ext in (".dd2", ".ved"):
        (tmp_path / (stem + ext)).write_bytes((src / (stem + ext)).read_bytes())
    v = pedro.analyze(p, pedro.Options(anharmonic="off", screen=False)).validation
    assert v["veda_version"].startswith("4.99")
    assert v["ted_fraction_within_1pp"] >= ted_within_1pp
    assert v["ted_cluster_summed_max_abs_diff_pp"] < 5
