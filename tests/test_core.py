import json
import math

import numpy as np
import pytest

import pedro
from pedro import basis as bs
from pedro import coordinates as co
from pedro import distributions as ds
from pedro import harmonic
from pedro.constants import FREQ_FACTOR
from pedro.io.gaussian_log import parse_archive

from conftest import write_set


# ----------------------------------------------------------------------------- B matrix
@pytest.mark.parametrize("kind,n", [("stretch", 2), ("bend", 3), ("torsion", 4), ("oop", 4)])
def test_analytic_b_rows_match_finite_differences(kind, n):
    rng = np.random.default_rng(42)
    worst = 0.0
    for _ in range(300):
        x = rng.normal(size=(6, 3)) * 1.5
        p = co.Primitive(kind, tuple(rng.permutation(6)[:n]))
        if kind == "bend" and abs(math.cos(co.value(p, x))) > 0.98:
            continue
        a, b = co.b_row(p, x), co.b_row_numeric(p, x, 1e-6)
        worst = max(worst, np.linalg.norm(a - b) / np.linalg.norm(b))
    assert worst < 1e-6


def test_b_rows_are_translation_and_rotation_invariant():
    rng = np.random.default_rng(1)
    x = rng.normal(size=(4, 3))
    for kind in ("stretch", "bend", "torsion", "oop"):
        p = co.Primitive(kind, tuple(range(co.NATOMS[kind])))
        row = co.b_row(p, x).reshape(-1, 3)
        assert np.allclose(row.sum(0), 0, atol=1e-10)                       # translations
        assert np.allclose(np.cross(x, row).sum(0), 0, atol=1e-9)           # infinitesimal rotations


# ----------------------------------------------------------------------------- harmonic analysis
def test_frequencies_equal_gf_eigenvalues(water):
    job, coords, F = water
    m = harmonic.analyze(job.hessian, job.masses, job.xyz_bohr)
    B = co.b_matrix(coords, job.xyz_bohr)
    G = (B / np.repeat(job.masses, 3)[None, :]) @ B.T
    lam = np.sort(np.linalg.eigvals(G @ F).real)
    assert np.allclose(m.eigenvalues, lam, rtol=1e-9)
    assert m.n == 3 and m.n_rigid == 6


def test_ir_intensities_are_mode_phase_independent(water):
    job, *_ = water
    m = harmonic.analyze(job.hessian, job.masses, job.xyz_bohr, job.dipole_derivatives)
    dmu = (-m.Lx).T @ job.dipole_derivatives
    assert np.allclose(m.ir_km_mol, 974.8801 * (dmu**2).sum(1))


# ----------------------------------------------------------------------------- distributions
def test_force_field_is_recovered_and_ted_closes(water):
    job, coords, F = water
    m = harmonic.analyze(job.hessian, job.masses, job.xyz_bohr)
    L = ds.mode_matrix(co.b_matrix(coords, job.xyz_bohr), m.Lx)
    Linv, cond = ds.invert(L)
    assert np.allclose(ds.force_field(Linv, m.eigenvalues), F, atol=1e-10)
    T = ds.ted(L, Linv)
    assert np.allclose(T.sum(1), 100) and np.allclose(T.sum(0), 100)
    # PED (row partition, from F) equals TED (from L alone): docs/theory.md section 4
    Fr = ds.force_field(Linv, m.eigenvalues)
    ped_row = (100 * L * (Fr @ L) / m.eigenvalues[None, :]).T
    assert np.allclose(ped_row, T, atol=1e-9)
    # Morino-Kuchitsu diagonal + pair terms close to 100 for each mode
    mk = ds.morino_kuchitsu(L, Fr, m.eigenvalues)
    for k in range(3):
        assert math.isclose(mk[k].sum() + ds.pair_terms(L, Fr, m.eigenvalues, k).sum(), 100, abs_tol=1e-9)


def test_ted_is_invariant_to_coordinate_scaling_and_mode_sign(water):
    job, coords, _ = water
    m = harmonic.analyze(job.hessian, job.masses, job.xyz_bohr)
    B = co.b_matrix(coords, job.xyz_bohr)
    L = B @ m.Lx
    T = ds.ted(L, np.linalg.inv(L))
    D = np.diag([3.0, -0.2, 7.0])
    S = np.diag([1, -1, 1])
    L2 = D @ L @ S
    assert np.allclose(ds.ted(L2, np.linalg.inv(L2)), T)


def test_group_total_invariant_under_reparametrisation_inside_group(water):
    job, coords, _ = water
    m = harmonic.analyze(job.hessian, job.masses, job.xyz_bohr)
    L = co.b_matrix(coords, job.xyz_bohr) @ m.Lx
    T = ds.ted(L, np.linalg.inv(L))
    s2 = 1 / math.sqrt(2)
    C = np.array([[s2, s2, 0], [s2, -s2, 0], [0, 0, 1]])        # symmetric/antisymmetric stretches
    L2 = C @ L
    T2 = ds.ted(L2, np.linalg.inv(L2))
    assert np.allclose(T[:, :2].sum(1), T2[:, :2].sum(1))
    assert not np.allclose(T[:, 0], T2[:, 0])


def test_veda_balancing_has_unit_marginals():
    rng = np.random.default_rng(5)
    T = rng.normal(size=(6, 6)) * 30
    P, info = ds.veda_compatible_ped(T, np.sign(T))
    assert info["converged"]
    assert np.allclose(np.abs(P).sum(0), 100, atol=1e-5) and np.allclose(np.abs(P).sum(1), 100, atol=1e-5)
    assert np.array_equal(np.sign(P), np.sign(T))


def test_clusters():
    cid = ds.clusters(np.array([100.0, 100.4, 102.0, 300.0, 300.9, 301.7]), 1.0)
    assert list(cid) == [0, 0, 1, 2, 2, 2]     # consecutive gaps < 1 cm-1 chain into one cluster


# ----------------------------------------------------------------------------- embedded bases
@pytest.mark.parametrize("case", ["Fprime", "P"])
def test_block_diagonal_identity(embedded, case, tmp_path):
    job, coords = embedded
    opts = pedro.Options(coords=write_set(tmp_path, coords, job.elements), subsystem=[0, 1, 2], case=case,
                         mobile=list(range(6)) if case == "P" else None)
    res = pedro.analyze(job, opts)
    assert res.checks["identity_max_dev_pp"] < 1e-8
    assert res.checks["ted_closure_max_pp"] < 1e-9


def test_case_F_is_square_and_closes(embedded, tmp_path):
    job, coords = embedded
    res = pedro.analyze(job, pedro.Options(coords=write_set(tmp_path, coords, job.elements), subsystem=[0, 1, 2]))
    assert res.basis.case == "F"
    assert res.ted.shape == (12, 12)
    assert res.checks["ted_closure_max_pp"] < 1e-9


def test_partial_hessian_matches_block_diagonalisation(embedded):
    job, _ = embedded
    m = harmonic.analyze(job.hessian, job.masses, job.xyz_bohr, method="partial", mobile=[0, 1, 2])
    c = np.arange(9)
    ms = np.repeat(np.sqrt(job.masses[:3]), 3)
    lam = np.linalg.eigvalsh(job.hessian[np.ix_(c, c)] / np.outer(ms, ms))
    assert np.allclose(np.sort(m.eigenvalues), lam)


def test_rigid_rows_are_exact_for_rigid_motions(embedded):
    job, _ = embedded
    R, _ = bs.rigid_rows([0, 1, 2], job.xyz_bohr, job.masses, 18)
    d = np.zeros((6, 3))
    d[:3] = [0.1, -0.2, 0.3]                                      # translate the subsystem only
    t = R @ d.reshape(-1)
    assert np.allclose(t[:3], [0.1, -0.2, 0.3]) and np.allclose(t[3:], 0, atol=1e-12)


# ----------------------------------------------------------------------------- I/O
def test_dd2_round_trip(water, tmp_path):
    job, coords, _ = water
    p = tmp_path / "w.dd2"
    co.export_dd2(coords, job.elements, job.xyz_bohr, p)
    back = co.import_dd2(p, job.elements)
    assert [(c.components[0][0], c.components[0][1]) for c in back] == [(c.components[0][0], 1.0) for c in coords]


def test_archive_parser_on_synthetic_log(water):
    job, *_ = water
    n3 = 9
    lt = job.hessian[np.tril_indices(n3)]
    xa = job.xyz_bohr * 0.529177210903
    geo = "\\".join(f"{e},{x:.8f},{y:.8f},{z:.8f}" for e, (x, y, z) in zip(job.elements, xa))
    arc = (f"1\\1\\GINC-TEST\\Freq\\RB3LYP\\6-31G(d)\\H2O1\\USER\\12-Nov-2026\\0\\\\# freq b3lyp/6-31g(d)\\\\water\\\\0,1\\"
           f"{geo}\\\\Version=TEST\\NImag=0\\\\{','.join(f'{v:.8f}' for v in lt)}\\\\0.,0.,0.,0.,0.,0.,0.,0.,0.\\\\\\@")
    text = "\n".join(" " + arc[i:i + 70] for i in range(0, len(arc), 70))
    a = parse_archive("header\n" + text + "\n Normal termination of Gaussian\n")
    assert a["elements"] == job.elements
    assert np.allclose(a["hessian"], job.hessian, atol=1e-7)
    assert np.allclose(a["xyz_angstrom"], xa, atol=1e-7)


def test_npz_end_to_end_never_overwrites(water, tmp_path):
    job, *_ = water
    p = tmp_path / "w.npz"
    np.savez(p, elements=np.array(job.elements), xyz=job.xyz_bohr, hessian=job.hessian,
             dipole_derivatives=job.dipole_derivatives)
    res = pedro.analyze(p)
    out = res.write()
    for f in ("modes.csv", "ted_subsystem.csv", "manifest.json", "summary.txt", "report.html", "coordinates.dd2"):
        assert (out / f).exists()
    out2 = res.write()
    assert out2 != out
    man = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    assert man["checks"]["case"] == "isolated"
    assert abs(res.modes.freqs[-1] - math.sqrt(res.modes.eigenvalues[-1]) * FREQ_FACTOR) < 1e-9
