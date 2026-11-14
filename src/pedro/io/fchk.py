"""Gaussian formatted-checkpoint reader (also reads VEDA .fmu files, which use the same layout).

Sections are parsed by label, never by line position. The Hessian comes only from 'Cartesian Force Constants';
'Internal Force Constants' (Gaussian's optimizer data) is never used as a substitute. A missing, truncated or
dimensionally inconsistent Hessian is a hard stop.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from ..constants import SYMBOLS
from ..job import Job
from .common import ReadError, sha256

_HEADER = re.compile(r"^(\S.{0,40}?)\s+([IRCL])\s+(?:N=\s*(\d+)|(\S+))\s*$")


def read_sections(path) -> dict:
    lines = Path(path).read_text(errors="replace").splitlines()
    sections: dict = {}
    i = 0
    while i < len(lines):
        m = _HEADER.match(lines[i])
        if not m:
            i += 1
            continue
        label, kind, count, scalar = m.group(1).strip(), m.group(2), m.group(3), m.group(4)
        if count is not None:
            n = int(count)
            j = i + 1
            vals: list[str] = []
            while j < len(lines) and lines[j][:1] == " ":
                if kind in "IR":
                    vals += lines[j].replace("D", "E").split()
                j += 1
            if kind in "IR":
                if len(vals) != n:
                    raise ReadError(f"{path}: section '{label}' is truncated: expected {n} values, found {len(vals)}")
                sections[label] = np.array(vals, dtype=int if kind == "I" else float)
            i = j
        else:
            if kind == "I":
                sections[label] = int(scalar)
            elif kind == "R":
                sections[label] = float(scalar.replace("D", "E"))
            else:
                sections[label] = scalar
            i += 1
    return sections


def load(path, kind="fchk") -> Job:
    path = str(Path(path).resolve())
    S = read_sections(path)
    for req in ("Atomic numbers", "Current cartesian coordinates"):
        if S.get(req) is None:
            raise ReadError(f"{path}: required section '{req}' missing")
    if S.get("Cartesian Force Constants") is None:
        raise ReadError(f"{path}: no 'Cartesian Force Constants' section. Run the frequency job with the checkpoint "
                        "saved and use 'formchk', or drop the Gaussian .log file instead (its archive block holds "
                        "the Hessian).")
    Z = np.asarray(S["Atomic numbers"], int)
    n = len(Z)
    if n == 0 or np.any(Z < 1) or np.any(Z > len(SYMBOLS)):
        raise ReadError(f"{path}: atomic numbers outside the supported range 1..{len(SYMBOLS)}")
    elements = [SYMBOLS[z - 1] for z in Z]
    xyz = np.asarray(S["Current cartesian coordinates"], float)
    if xyz.size != 3 * n:
        raise ReadError(f"{path}: coordinate count does not match {n} atoms")
    packed = np.asarray(S["Cartesian Force Constants"], float)
    if packed.size != 3 * n * (3 * n + 1) // 2:
        raise ReadError(f"{path}: Cartesian Force Constants has {packed.size} values; expected "
                        f"{3 * n * (3 * n + 1) // 2} for {n} atoms")
    H = np.zeros((3 * n, 3 * n))
    H[np.tril_indices(3 * n)] = packed
    H = H + H.T - np.diag(np.diag(H))
    masses = S.get("Real atomic weights")
    if masses is not None and (np.asarray(masses).size != n or np.any(np.asarray(masses) <= 0)):
        raise ReadError(f"{path}: 'Real atomic weights' inconsistent with {n} atoms")
    dd = S.get("Dipole Derivatives")
    if dd is not None and np.asarray(dd).size != 9 * n:
        raise ReadError(f"{path}: Dipole Derivatives has {np.asarray(dd).size} values, expected {9 * n}")
    grad = S.get("Cartesian Gradient")
    layers = S.get("Atom Layers")
    layer_map = {1: "H", 2: "M", 3: "L"}
    job = Job(
        source=path, sha256=sha256(path), kind=kind, elements=elements,
        masses=np.asarray(masses, float) if masses is not None else None,
        mass_source="file (Real atomic weights)" if masses is not None else "",
        xyz_bohr=xyz.reshape(n, 3), hessian=H,
        dipole_derivatives=None if dd is None else np.asarray(dd, float).reshape(3 * n, 3),
        gradient=None if grad is None else np.asarray(grad, float),
        charge=S.get("Charge"), multiplicity=S.get("Multiplicity"),
        layers=[layer_map.get(int(v), "?") for v in layers] if layers is not None and len(set(map(int, layers))) > 1 else None,
    )
    if S.get("Anharmonic Vib-E2") is not None or S.get("Anharmonic Vib-E2 Info") is not None:
        job.anharmonic_present = True
    return job
