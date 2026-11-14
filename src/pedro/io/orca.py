"""ORCA Hessian file (.hess) reader.

ORCA writes `<job>.hess` for every frequency job. PEDRO takes from it:
    $hessian              Cartesian Hessian, Eh/bohr^2, printed in column blocks
    $atoms                element, mass (u, as ORCA used them) and coordinates (bohr)
    $dipole_derivatives   d(mu)/dx, e, one row per Cartesian coordinate (optional)
    $vibrational_frequencies, $ir_spectrum   ORCA's printed values, used only to check PEDRO's recomputation
    $multiplicity         (optional)
Masses are kept as ORCA wrote them (standard atomic weights by default), so the recomputed frequencies match ORCA's.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..job import Job
from .common import ReadError, sha256


def _sections(text):
    out, name, buf = {}, None, []
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("$"):
            if name is not None:
                out[name] = buf
            name, buf = s[1:].split()[0].lower() if len(s) > 1 else "", []
        elif name is not None and s and not s.startswith("#"):
            buf.append(s)
    if name is not None:
        out[name] = buf
    return out


def _blocked_matrix(lines, name):
    """A square matrix printed as column blocks: a header line of column indices, then 'row v v v ...'."""
    try:
        n = int(lines[0].split()[0])
    except (IndexError, ValueError):
        raise ReadError(f"${name}: missing dimension line")
    M = np.full((n, n), np.nan)
    cols = None
    for line in lines[1:]:
        t = line.split()
        if all(v.lstrip("-").isdigit() for v in t) and "." not in line:
            cols = [int(v) for v in t]
            continue
        if cols is None:
            raise ReadError(f"${name}: data before a column header")
        r = int(t[0])
        vals = [float(v) for v in t[1:]]
        if len(vals) != len(cols):
            raise ReadError(f"${name}: row {r} has {len(vals)} values for {len(cols)} columns")
        M[r, cols] = vals
    if np.isnan(M).any():
        raise ReadError(f"${name}: matrix incomplete ({int(np.isnan(M).sum())} missing entries)")
    return M


def load(path) -> Job:
    path = str(Path(path).resolve())
    sec = _sections(Path(path).read_text(errors="replace"))
    if "orca_hessian_file" not in sec:
        raise ReadError(f"{path}: not an ORCA Hessian file (no $orca_hessian_file line)")
    for need in ("hessian", "atoms"):
        if need not in sec:
            raise ReadError(f"{path}: no ${need} section")
    at = sec["atoms"]
    n = int(at[0].split()[0])
    rows = [l.split() for l in at[1:1 + n]]
    if len(rows) != n or any(len(r) < 5 for r in rows):
        raise ReadError(f"{path}: $atoms section is incomplete")
    el = [r[0].capitalize() for r in rows]
    mass = np.array([float(r[1]) for r in rows])
    xyz = np.array([[float(v) for v in r[2:5]] for r in rows])
    H = _blocked_matrix(sec["hessian"], "hessian")
    if H.shape != (3 * n, 3 * n):
        raise ReadError(f"{path}: Hessian is {H.shape[0]}x{H.shape[1]}, expected {3 * n} for {n} atoms")
    dd = None
    if "dipole_derivatives" in sec:
        v = sec["dipole_derivatives"]
        dd = np.array([[float(x) for x in l.split()[:3]] for l in v[1:1 + 3 * n]])
        if dd.shape != (3 * n, 3):
            dd = None
    pf = pir = None
    if "vibrational_frequencies" in sec:
        f = np.array([float(l.split()[1]) for l in sec["vibrational_frequencies"][1:]])
        keep = np.abs(f) > 1e-6                            # ORCA prints rigid-body modes as exactly 0
        pf = f[keep]
        if "ir_spectrum" in sec:
            ir = np.array([[float(x) for x in l.split()[:3]] for l in sec["ir_spectrum"][1:]])
            if len(ir) == len(f):
                pir = ir[keep, 2]
    mult = int(sec["multiplicity"][0].split()[0]) if sec.get("multiplicity") else None
    job = Job(source=path, sha256=sha256(path), kind="orca_hess", elements=el, masses=mass,
              mass_source="file (masses as ORCA used them)", xyz_bohr=xyz, hessian=(H + H.T) / 2,
              dipole_derivatives=dd, multiplicity=mult, printed_freqs=pf, printed_ir=pir, program="ORCA",
              termination="normal")
    if "act_atom" in sec and sec["act_atom"] and sec["act_atom"][0].split()[0] not in ("0",):
        job.notes.append("ORCA $act_atom is set; check whether this is a partial-Hessian job")
    return job
