"""Coordinate transfer by graph isomorphism, coordinate classes, mass options and band metrics."""

import json

import numpy as np
import pytest

import pedro
from pedro import coordinates as co
from pedro import mapping, metrics

from conftest import make_job, valence_hessian

ANG = 1 / 0.529177210903


def ethanal():
    """CH3-CHO: a methyl group gives 6 isomorphisms (H permutations); only one superposes."""
    el = ["C", "C", "O", "H", "H", "H", "H"]
    x = np.array([[0.0, 0.0, 0.0], [1.50, 0.0, 0.0], [2.10, 1.05, 0.0], [2.05, -0.95, 0.0],
                  [-0.38, 1.02, 0.0], [-0.38, -0.51, 0.88], [-0.38, -0.51, -0.88]]) * ANG
    return el, x


def test_isomorphism_recovers_scrambled_order():
    el, x = ethanal()
    rng = np.random.default_rng(4)
    perm = rng.permutation(len(el))                      # target atom perm[i] holds source atom i
    tgt_el = [None] * len(el)
    tgt_x = np.zeros_like(x)
    for i, p in enumerate(perm):
        tgt_el[p], tgt_x[p] = el[i], x[i]
    R = np.linalg.qr(rng.normal(size=(3, 3)))[0]
    tgt_x = tgt_x @ R.T + 3.0
    m = mapping.map_molecule(el, x, tgt_el, tgt_x, list(range(len(el))))
    assert m["n_isomorphisms"] == 6
    assert list(m["map"]) == list(perm) and m["rmsd_angstrom"] < 1e-6


def test_isomorphism_refuses_different_molecules():
    el, x = ethanal()
    el2 = list(el)
    el2[2] = "N"
    with pytest.raises(mapping.MappingError):
        mapping.map_molecule(el, x, el2, x, list(range(len(el))))


def test_transfer_inside_larger_structure(embedded, tmp_path):
    """A water coordinate set defined on an isolated water is placed on the water inside the 6-atom model."""
    job, coords = embedded
    iso = tmp_path / "w.npz"
    np.savez(iso, elements=np.array(job.elements[:3]), xyz=job.xyz_bohr[:3] + 1.0, hessian=np.eye(9))
    dd2 = tmp_path / "w.dd2"
    co.export_dd2(coords, job.elements[:3], job.xyz_bohr[:3], dd2)
    res = pedro.analyze(job, pedro.Options(coords=str(dd2), coords_molecule=str(iso), subsystem=[0, 1, 2], case="Fprime"))
    tr = res.coord_report["transfer"]
    assert tr["source_to_target_1based"] in ([1, 2, 3], [1, 3, 2]) and tr["rmsd_angstrom"] < 1e-6
    assert res.checks["identity_max_dev_pp"] < 1e-8


def test_coordinate_classes_and_metrics(water, tmp_path):
    job, coords, _ = water
    p = tmp_path / "set.json"
    p.write_text(json.dumps(co.export_json(coords, job.elements)), encoding="utf-8")
    res = pedro.analyze(job, pedro.Options(coords=str(p), coord_classes=[("OH1", "stretch", [0, 1]),
                                                                          ("OH2", "stretch", [0, 2]), ("rest", "rest", [])]))
    assert res.basis.groups[:3] == ["OH1", "OH2", "rest"]
    res2 = pedro.analyze(job, pedro.Options(coords=str(p), coord_classes=[("H2", "any", [2]), ("rest", "rest", [])]))
    assert res2.basis.groups[:3] == ["rest", "H2", "H2"]          # the bend touches H2 as well
    j = metrics.stretch_index(res, 0, 1)
    assert j == 0
    bm = metrics.band_metrics(res, j, (0, 1e5))
    k = bm["strongest"]["index"]
    w = res.ted[:, j] * res.modes.ir_km_mol
    assert k == int(np.argmax(w))
    assert bm["centroid_ted_ir"] == pytest.approx((w * res.modes.freqs).sum() / w.sum())
    tot = metrics.class_totals(res, k)
    assert sum(tot.values()) == pytest.approx(100.0)


def test_mass_override_and_following(water, tmp_path):
    job, coords, _ = water
    p = tmp_path / "set.json"
    p.write_text(json.dumps(co.export_json(coords, job.elements)), encoding="utf-8")
    a = pedro.analyze(job, pedro.Options(coords=str(p)))
    b = pedro.analyze(job, pedro.Options(coords=str(p), mass_overrides={1: 2.01410, 2: 2.01410}))
    assert "overridden" in b.job.mass_source
    assert np.all(b.modes.freqs < a.modes.freqs)                       # D2O is lower
    assert np.allclose(job.masses[1], 1.00782503223)                    # caller's job is unchanged
    m = metrics.match(a, b, [0, 1, 2])
    assert sorted(x["b"] for x in m) == [1, 2, 3] and all(x["overlap"] > 0.8 for x in m)
    assert m[0]["b"] == 1                                                # the bend stays the lowest mode
    f = metrics.follow(a, b, 0, [0, 1, 2])[0]
    assert f["index"] == 0 and f["overlap"] > 0.9


def test_printed_masses_option(water):
    job, *_ = water
    with pytest.raises(ValueError):
        pedro.analyze(job, pedro.Options(masses="printed"))
    job.printed_masses = np.round(job.masses, 5)
    r = pedro.analyze(job, pedro.Options(masses="printed"))
    assert np.allclose(r.job.masses, np.round(job.masses, 5))


def test_band_metrics_without_ir(water, tmp_path):
    job, coords, _ = water
    job.dipole_derivatives = None
    p = tmp_path / "set.json"
    p.write_text(json.dumps(co.export_json(coords, job.elements)), encoding="utf-8")
    res = pedro.analyze(job, pedro.Options(coords=str(p)))
    bm = metrics.band_metrics(res, 0, (0, 1e5))
    assert bm["weights"].startswith("TED only") and bm["centroid_ted_ir"] is None
    assert bm["strongest"]["index"] == int(np.argmax(res.ted[:, 0]))
