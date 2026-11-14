"""Program-independent Hessian input (.npz or .json), so ORCA, Psi4, PySCF or Q-Chem users can run PEDRO.

Schema 'pedro-hessian/1' (all arrays in atomic units except masses):
    elements    list of element symbols, or Z  atomic numbers
    masses      u (optional; most-abundant isotopes otherwise)
    xyz         bohr, shape (N, 3)
    hessian     Eh/bohr^2, shape (3N, 3N)
    dipole_derivatives  e, shape (3N, 3)   (optional)
    gradient    Eh/bohr, shape (3N,)       (optional)
    layers      list of 'H'/'M'/'L'        (optional)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ..constants import ISOTOPE_MASS, SYMBOLS
from ..job import Job
from .common import ReadError, sha256

SCHEMA = "pedro-hessian/1"


def _get(d, key):
    v = d.get(key) if isinstance(d, dict) else (d[key] if key in d.files else None)
    return None if v is None else np.asarray(v)


def load(path) -> Job:
    path = str(Path(path).resolve())
    if path.lower().endswith(".npz"):
        d = np.load(path, allow_pickle=False)
        kind = "npz"
    else:
        d = json.loads(Path(path).read_text(encoding="utf-8"))
        kind = "json"
        if d.get("format") not in (None, SCHEMA):
            raise ReadError(f"{path}: format '{d.get('format')}' is not {SCHEMA}")
    el = _get(d, "elements")
    if el is None:
        Z = _get(d, "Z")
        if Z is None:
            raise ReadError(f"{path}: needs 'elements' or 'Z'")
        el = [SYMBOLS[int(z) - 1] for z in Z]
    el = [str(e) for e in el]
    n = len(el)
    xyz = _get(d, "xyz")
    H = _get(d, "hessian")
    if xyz is None or H is None:
        raise ReadError(f"{path}: needs 'xyz' (bohr) and 'hessian' (Eh/bohr^2)")
    xyz = xyz.astype(float).reshape(n, 3)
    H = H.astype(float).reshape(3 * n, 3 * n)
    m = _get(d, "masses")
    if m is None:
        if any(e not in ISOTOPE_MASS for e in el):
            raise ReadError(f"{path}: masses missing and no default mass for some elements")
        m, src = np.array([ISOTOPE_MASS[e] for e in el]), "default most-abundant isotopes"
    else:
        m, src = m.astype(float), "file"
    dd = _get(d, "dipole_derivatives")
    g = _get(d, "gradient")
    layers = _get(d, "layers")
    return Job(source=path, sha256=sha256(path), kind=kind, elements=el, masses=m, mass_source=src, xyz_bohr=xyz,
               hessian=(H + H.T) / 2, dipole_derivatives=None if dd is None else dd.astype(float).reshape(3 * n, 3),
               gradient=None if g is None else g.astype(float).reshape(-1),
               layers=None if layers is None else [str(v) for v in layers])
