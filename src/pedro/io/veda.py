"""VEDA 4 files: the .fmv input (VEDA's own Hessian file) and the .ved output, which PEDRO compares with its own
results on the same input."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from ..constants import ANG_TO_BOHR, SYMBOLS
from ..job import Job
from .common import ReadError, sha256

# VEDA writes the .fmv Hessian in mdyn/A. Measured against Gaussian's printed Hessian for VEDA's own example a1-91:
# fmv / Gaussian = 15.5691904 with a scatter of 2e-8, so dividing by this factor recovers Gaussian's Eh/bohr^2 values.
FMV_MDYN_PER_ANG_PER_AU = 15.56919


def write_fmu(job, path, atoms=None, title=None):
    """Write a VEDA-readable .fmu (formatted-checkpoint layout) for `job`, or for the Hessian block of `atoms`
    (0-based) when given. Refuses to overwrite. Lets VEDA analyze exactly the Hessian PEDRO read."""
    from ..constants import SYMBOLS as _S
    p = Path(path)
    atoms = list(range(job.natoms)) if atoms is None else list(atoms)
    n = len(atoms)
    rows = (3 * np.asarray(atoms)[:, None] + np.arange(3)).reshape(-1)
    H = job.hessian[np.ix_(rows, rows)]
    Z = [_S.index(job.elements[a]) + 1 for a in atoms]
    mass = (job.printed_masses if job.printed_masses is not None else job.masses)[atoms]
    with open(p, "x", encoding="ascii", newline="\r\n") as out:
        out.write(f"{(title or Path(job.source).name)[:72]}\n")
        out.write(f"{'Freq':<10}{'PEDRO export':<60}\n")
        out.write(f'{"Number of atoms":<43}I{n:17d}\n')
        for label, values, typ in (("Atomic numbers", Z, "I"), ("Current cartesian coordinates", job.xyz_bohr[atoms].ravel(), "R"),
                                   ("Real atomic weights", mass, "R"), ("Cartesian Force Constants", H[np.tril_indices(3 * n)], "R")):
            out.write(f"{label:<43}{typ}   N={len(values):12d}\n")
            w = 6 if typ == "I" else 5
            for a in range(0, len(values), w):
                out.write("".join(f"{int(v):12d}" if typ == "I" else f"{v:16.8E}" for v in values[a:a + w]) + "\n")
    return p


def load_fmv(path) -> Job:
    """VEDA .fmv: 'Atomic orientation and masses' (n, Z, x y z in angstrom, mass in u), 'Force constants in Cartesian
    coordinates' (full 3N x 3N matrix, mdyn/A), optional 'Harmonic frequencies' ('F = value  IR = value')."""
    path = str(Path(path).resolve())
    lines = Path(path).read_text(errors="replace").splitlines()
    try:
        a = next(i for i, l in enumerate(lines) if l.strip().startswith("Atomic orientation and masses"))
        f = next(i for i, l in enumerate(lines) if l.strip().startswith("Force constants in Cartesian coordinates"))
    except StopIteration:
        raise ReadError(f"{path}: not a VEDA .fmv file (needs 'Atomic orientation and masses' and 'Force constants')")
    rows = [l.split() for l in lines[a + 1:f] if l.strip()]
    try:
        Z = [int(r[1]) for r in rows]
        xyz = np.array([[float(v) for v in r[2:5]] for r in rows]) * ANG_TO_BOHR
        mass = np.array([float(r[5]) for r in rows])
    except (IndexError, ValueError):
        raise ReadError(f"{path}: cannot read the atom block")
    n = len(Z)
    vals, freqs, irs = [], [], []
    for l in lines[f + 1:]:
        s = l.strip()
        if s.startswith("Harmonic"):
            continue
        m = re.match(r"F\s*=\s*(-?\d+\.\d+)\s+IR\s*=\s*(-?\d+\.\d+)", s)
        if m:
            freqs.append(float(m.group(1)))
            irs.append(float(m.group(2)))
        elif not freqs and s:
            try:
                vals += [float(v.replace("D", "E")) for v in s.split()]
            except ValueError:
                pass
    if len(vals) != (3 * n) ** 2:
        raise ReadError(f"{path}: {len(vals)} force constants for {n} atoms (expected {(3 * n) ** 2})")
    H = np.array(vals).reshape(3 * n, 3 * n) / FMV_MDYN_PER_ANG_PER_AU
    job = Job(source=path, sha256=sha256(path), kind="veda_fmv", elements=[SYMBOLS[z - 1] for z in Z], masses=mass,
              mass_source="file (masses as VEDA wrote them)", xyz_bohr=xyz, hessian=(H + H.T) / 2, program="VEDA .fmv",
              printed_freqs=np.sort(np.array(freqs)) if len(freqs) == 3 * n - 6 else None, termination="normal")
    return job


def _matrix(lines, start, end, n):
    rows = []
    for line in lines[start:end]:
        t = line.split()
        if len(t) == n + 2 and t[0].isdigit() and t[-1] == t[0]:
            rows.append([float(v) for v in t[1:-1]])
    a = np.array(rows)
    if a.shape != (n, n):
        raise ValueError(f"VEDA table has shape {a.shape}, expected ({n}, {n})")
    return a


def read_ved(path) -> dict:
    """Frequencies (cm^-1, VEDA's descending order), PED and TED (integer %, rows = modes in that order)."""
    lines = Path(path).read_text(errors="replace").splitlines()
    pi = next(i for i, l in enumerate(lines) if l.startswith("PED:"))
    ti = next(i for i, l in enumerate(lines) if l.startswith("TED:"))
    freqs = []
    for line in lines[:pi]:
        if re.fullmatch(r"\s*(?:-?\d+\.\d+\s*)+", line):
            freqs.extend(float(v) for v in line.split())
    n = len(freqs)
    return {"freqs": np.array(freqs), "ped": _matrix(lines, pi + 1, ti, n), "ted": _matrix(lines, ti + 1, len(lines), n),
            "version": (re.search(r"version\s+(\S+)", lines[0]) or [None, "?"])[1]}
