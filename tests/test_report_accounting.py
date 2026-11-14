"""Regression checks for the scientific meaning of exported subsystem/environment budgets."""

import csv
import json

import numpy as np
import pytest

import pedro
from pedro import report

from conftest import write_set


@pytest.mark.parametrize("case,mobile", [("F", None), ("Fprime", None), ("P", None), ("P", list(range(6)))])
def test_export_keeps_subsystem_rigid_motion_out_of_environment(embedded, tmp_path, monkeypatch, case, mobile):
    job, coords = embedded
    res = pedro.analyze(job, pedro.Options(coords=write_set(tmp_path, coords, job.elements),
                                        subsystem=[0, 1, 2], case=case, mobile=mobile,
                                        anharmonic="off", screen=False))
    # Plotting is independent of the numerical export under test.
    monkeypatch.setattr(report, "figures", lambda _: {})
    out = res.write(tmp_path / "output")
    with (out / "ted_subsystem.csv").open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        names, rows = reader.fieldnames, list(reader)
    env = np.array([float(row["environment_total"]) for row in rows])
    ni = res.basis.n_internal
    summed = res.ted[:, :ni].sum(1) + env
    if case in ("Fprime", "P"):
        rigid = np.array([float(row["subsystem_rigid_body_total"]) for row in rows])
        np.testing.assert_allclose(rigid, res.group_totals["subsystem rigid-body"], atol=0.00051)
        summed += rigid
        # Independent atom-amplitude identity: compare with the environment eigenvector components.
        env_atoms = [3, 4, 5] if case == "Fprime" or mobile else []
        q = res.modes.Lq.reshape(job.natoms, 3, -1)
        expected_env = 100 * (q[env_atoms] ** 2).sum(axis=(0, 1))
        np.testing.assert_allclose(env, expected_env, atol=0.00051)
        assert "S rigid + environment" not in " ".join(res.descriptions)
        if case == "P" and mobile is None:
            np.testing.assert_array_equal(env, 0)
            assert rigid.max() > 50  # a large subsystem rigid motion must never be called environment
            assert all("[environment" not in d for d in res.descriptions)
    else:
        assert "subsystem_rigid_body_total" not in names
        np.testing.assert_allclose(env, res.ted[:, ni:].sum(1), atol=0.00051)
    np.testing.assert_allclose(summed, 100, atol=0.0011)
    expected_label = "case F; approximate" if case == "F" else "exact identity"
    for filename in ("summary.txt", "report.html"):
        assert expected_label in (out / filename).read_text(encoding="utf-8")
    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["ted_subsystem_columns"]["environment_total"]


def test_isolated_export_has_no_environment_column(water, tmp_path, monkeypatch):
    job, *_ = water
    monkeypatch.setattr(report, "figures", lambda _: {})
    out = pedro.analyze(job).write(tmp_path / "output")
    with (out / "ted_subsystem.csv").open(encoding="utf-8-sig", newline="") as f:
        names = next(csv.reader(f))
    assert "environment_total" not in names and "subsystem_rigid_body_total" not in names
