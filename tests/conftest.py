"""Synthetic test systems with known answers. They need no external program or private data."""

import json
import os
from pathlib import Path

import numpy as np
import pytest

from pedro import coordinates as co
from pedro.constants import ISOTOPE_MASS
from pedro.job import Job

ANG = 1 / 0.529177210903


def make_job(elements, xyz_bohr, H, dd=None, layers=None, name="synthetic"):
    return Job(source=name, sha256="0" * 64, kind="synthetic", elements=list(elements),
               masses=np.array([ISOTOPE_MASS[e] for e in elements]), mass_source="test", xyz_bohr=np.asarray(xyz_bohr, float),
               hessian=H, dipole_derivatives=dd, layers=layers)


def valence_hessian(coords, x, F):
    """Cartesian Hessian of V = 1/2 ds^T F ds at its minimum: H = B^T F B (no gradient term)."""
    B = co.b_matrix(coords, x)
    return B.T @ F @ B


@pytest.fixture
def water():
    x = np.array([[0.0, 0.0, 0.0], [0.757, 0.587, 0.0], [-0.757, 0.587, 0.0]]) * ANG
    el = ["O", "H", "H"]
    P = co.Primitive
    coords = [co.Coordinate("r1", "stretch", [(P("stretch", (0, 1)), 1.0)]),
              co.Coordinate("r2", "stretch", [(P("stretch", (0, 2)), 1.0)]),
              co.Coordinate("a", "bend", [(P("bend", (1, 0, 2)), 1.0)])]
    F = np.array([[0.50, 0.01, 0.03], [0.01, 0.50, 0.03], [0.03, 0.03, 0.16]])
    H = valence_hessian(coords, x, F)
    rng = np.random.default_rng(3)
    dd = rng.normal(scale=0.3, size=(9, 3))
    return make_job(el, x, H, dd), coords, F


@pytest.fixture
def embedded():
    """Water (subsystem) coupled by springs to three environment atoms (an Ar-like heavy cage)."""
    rng = np.random.default_rng(11)
    xw = np.array([[0.0, 0.0, 0.0], [0.757, 0.587, 0.0], [-0.757, 0.587, 0.0]]) * ANG
    xe = np.array([[0.0, -2.8, 0.3], [2.6, 1.9, -0.4], [-2.5, 2.0, 0.5]]) * ANG
    x = np.vstack([xw, xe])
    el = ["O", "H", "H", "S", "S", "S"]
    P = co.Primitive
    coords = [co.Coordinate("", "stretch", [(P("stretch", (0, 1)), 1.0)]),
              co.Coordinate("", "stretch", [(P("stretch", (0, 2)), 1.0)]),
              co.Coordinate("", "bend", [(P("bend", (1, 0, 2)), 1.0)])]
    springs = [co.Coordinate("", "stretch", [(P("stretch", (i, j)), 1.0)]) for i in range(6) for j in range(3, 6) if i < j]
    allc = coords + springs
    F = np.diag([0.5, 0.5, 0.16] + list(rng.uniform(0.01, 0.05, len(springs))))
    H = valence_hessian(allc, x, F)
    return make_job(el, x, H), coords


@pytest.fixture
def veda_dir():
    """VEDA 4.99e runs on our own Hessians (validation/veda_reference), or another folder of the same layout given by
    PEDRO_VEDA_DATA."""
    p = os.environ.get("PEDRO_VEDA_DATA") or Path(__file__).resolve().parents[1] / "validation" / "veda_reference"
    if not Path(p).exists():
        pytest.skip(f"no VEDA reference folder at {p}")
    return Path(p)


def write_set(tmp_path, coords, elements):
    p = tmp_path / "set.json"
    p.write_text(json.dumps(co.export_json(coords, elements)), encoding="utf-8")
    return str(p)
