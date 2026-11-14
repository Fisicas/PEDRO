"""Gaussian .log/.out reader.

The Hessian, geometry, dipole derivatives and forces are read from the archive entry at the end of the job
(the block that starts '1\\1\\' and ends '\\\\@'). Gaussian writes it for every frequency job, including ONIOM
jobs, so no checkpoint file is needed. The archive geometry is in the input orientation and the Hessian is in
the same frame; the printed normal-mode vectors (standard orientation) are never mixed with it.

The body of the log is read for validation and context only: termination, route, ONIOM layers and frozen
flags from the input echo, isotope masses, printed frequencies and IR intensities, and whether a GVPT2
(Freq=Anharmonic) section exists.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from ..constants import ANG_TO_BOHR, ISOTOPE_MASS
from ..job import Job
from .common import ReadError, sha256

_ARCHIVE_START = re.compile(r"^ 1([\\|])1\1", re.M)
_NUM = r"-?\d+\.\d+"


# ----------------------------------------------------------------------------- archive block
def _archive_blocks(text):
    """Yield (separator, joined archive string) for every archive entry, last one first."""
    starts = list(_ARCHIVE_START.finditer(text))
    for m in reversed(starts):
        sep = m.group(1)
        stop = text.find("\n\n", m.start())          # the archive is one block of 70-column lines
        raw = text[m.start():len(text) if stop < 0 else stop]
        joined = "".join(line[1:] for line in raw.splitlines())
        end = joined.find(sep + sep + "@")              # the closing '\\@' may wrap onto its own line
        if end >= 0:
            yield sep, joined[:end + 3]


def _floats(s):
    return np.array([float(v.replace("D", "E")) for v in s.split(",") if v.strip()], float)


def parse_archive(text, path="log"):
    """Return a dict with elements, xyz (angstrom), hessian, dipole derivatives, forces, route, version.

    The last archive entry that contains a complete Cartesian Hessian is used."""
    tried = 0
    for sep, arc in _archive_blocks(text):
        tried += 1
        secs = arc.split(sep + sep)
        if len(secs) < 5:
            continue
        head = secs[0].split(sep)
        geo_fields = secs[3].split(sep)
        charge_mult = [int(v) for v in re.findall(r"-?\d+", geo_fields[0])]
        elements, xyz = [], []
        for g in geo_fields[1:]:
            t = [x.strip() for x in g.split(",")]
            if len(t) < 4:
                continue
            sym = re.match(r"([A-Z][a-z]?)", t[0])
            if not sym:
                continue
            nums = t[2:5] if (len(t) >= 5 and re.fullmatch(r"-?\d+", t[1])) else t[1:4]
            elements.append(sym.group(1))
            xyz.append([float(v) for v in nums])
        n = len(elements)
        if n == 0:
            continue
        props = {}
        for kv in secs[4].split(sep):
            if "=" in kv:
                k, v = kv.split("=", 1)
                props[k.strip()] = v
        nh = 3 * n * (3 * n + 1) // 2
        hess = forces = None
        for s in secs[5:]:
            if not s or s.startswith("@"):
                continue
            try:
                vals = _floats(s)
            except ValueError:
                continue
            if vals.size == nh and hess is None:
                hess = vals
            elif vals.size == 3 * n and hess is not None and forces is None:
                forces = vals
        if hess is None:
            continue
        H = np.zeros((3 * n, 3 * n))
        H[np.tril_indices(3 * n)] = hess
        H = H + H.T - np.diag(np.diag(H))
        dd = None
        if "DipoleDeriv" in props:
            d = _floats(props["DipoleDeriv"])
            if d.size == 9 * n:
                dd = d.reshape(3 * n, 3)
        return {
            "elements": elements, "xyz_angstrom": np.array(xyz), "hessian": H, "dipole_derivatives": dd,
            "forces": forces, "route": secs[1].strip(), "title": secs[2].strip(),
            "charge_mult": charge_mult, "version": props.get("Version", ""),
            "job_type": head[3] if len(head) > 3 else "", "nimag": props.get("NImag"),
        }
    raise ReadError(
        f"{path}: no archive entry with a Cartesian Hessian was found ({tried} archive entries checked). "
        "PEDRO needs a frequency job (Freq). If the job finished but the archive was suppressed, convert the "
        "checkpoint with 'formchk' and load the .fchk file instead.")


# ----------------------------------------------------------------------------- body of the log
def _termination(text):
    tail = text[-5000:]
    if "Normal termination of Gaussian" in tail:
        return "normal"
    if "Error termination" in tail:
        return "error"
    return "incomplete"


def _input_atoms(text):
    """Atoms of the first 'Symbolic Z-matrix' echo: element, freeze flag and ONIOM layer (H/M/L) if present."""
    m = re.search(r"Symbolic Z-matrix:\n((?: Charge =[^\n]*\n)+)(.*?)\n\s*\n", text, re.S)
    if not m:
        return []
    atoms = []
    for line in m.group(2).splitlines():
        t = line.split()
        if not t or not re.match(r"[A-Za-z]", t[0]):
            continue
        rest = t[1:]
        freeze = None
        if len(rest) >= 4 and re.fullmatch(r"-?\d+", rest[0]):
            freeze = int(rest[0])
            rest = rest[1:]
        layer = rest[3] if len(rest) >= 4 and rest[3] in ("H", "M", "L") else None
        atoms.append({"freeze": freeze, "layer": layer})
    return atoms


def _printed_masses(text, n):
    vals = re.findall(r"^ Atom\s+\d+ has atomic number\s+\d+ and mass\s+(" + _NUM + ")", text, re.M)
    return [float(v) for v in vals[-n:]] if len(vals) >= n else None


def _printed_freqs(text):
    i = text.rfind("Harmonic frequencies (cm**-1)")
    if i < 0:
        return None, None
    tail = text[i:]
    freqs = [float(v) for l in re.findall(r"^\s*Frequencies --(?!-)(.*)$", tail, re.M) for v in l.split()]
    ir = [float(v) for l in re.findall(r"^\s*IR Inten\s+--(?!-)(.*)$", tail, re.M) for v in l.split()]
    return (np.array(freqs) if freqs else None), (np.array(ir) if len(ir) == len(freqs) and ir else None)


def _route(text):
    m = re.search(r"^ -{5,}\n( #[^\n]*\n(?: [^\n-][^\n]*\n)*?) -{5,}$", text, re.M)
    return "".join(l[1:] for l in m.group(1).splitlines()).strip() if m else ""


def load(path) -> Job:
    path = str(Path(path).resolve())
    text = Path(path).read_text(encoding="latin-1")
    arc = parse_archive(text, path)
    n = len(arc["elements"])
    notes = []

    masses = None
    printed = _printed_masses(text, n)
    default = [ISOTOPE_MASS.get(e) for e in arc["elements"]]
    if printed is not None:
        if all(d is not None and abs(p - d) < 6e-6 for p, d in zip(printed, default)):
            masses, mass_source = np.array(default), "log (printed masses match the most-abundant isotopes; exact values used)"
        else:
            masses, mass_source = np.array(printed), "log (printed isotope masses, 5 decimals)"
    elif all(d is not None for d in default):
        masses, mass_source = np.array(default), "default most-abundant isotopes (no masses printed in the log)"
        notes.append("Isotope masses were not printed in the log; most-abundant isotope masses were assumed.")
    else:
        missing = sorted({e for e, d in zip(arc["elements"], default) if d is None})
        raise ReadError(f"{path}: no printed masses and no default mass for {missing}; supply masses explicitly")

    atoms = _input_atoms(text)
    layers = frozen = None
    if len(atoms) == n:
        lay = [a["layer"] for a in atoms]
        if all(lay) and len(set(lay)) > 1:
            layers = lay
        frozen = [i for i, a in enumerate(atoms) if a["freeze"] is not None and a["freeze"] < 0]
    elif atoms:
        notes.append(f"Input echo lists {len(atoms)} atoms but the archive has {n}; ONIOM layers not read.")

    pf, pir = _printed_freqs(text)
    cm = arc["charge_mult"]
    charge, mult = (cm[0], cm[1]) if len(cm) >= 2 else (None, None)
    # ONIOM: the archive stores the real-system (low-level) charge; the subsystem of interest is the model system,
    # whose charge and multiplicity are printed in the input echo.
    hm = re.search(r"Charge =\s*(-?\d+) Multiplicity =\s*(\d+) for high\s+level calculation on model system", text)
    if hm:
        notes.append(f"ONIOM: real system charge {charge}, multiplicity {mult}; model system (high level) charge "
                     f"{hm.group(1)}, multiplicity {hm.group(2)} (reported).")
        charge, mult = int(hm.group(1)), int(hm.group(2))
    route = _route(text) or arc["route"]
    anh = ("Second-order Perturbative Anharmonic Analysis" in text) or ("Anharmonic Infrared Spectroscopy" in text)
    return Job(
        source=path, sha256=sha256(path), kind="gaussian_log", elements=arc["elements"], masses=masses,
        mass_source=mass_source, xyz_bohr=arc["xyz_angstrom"] * ANG_TO_BOHR, hessian=arc["hessian"],
        dipole_derivatives=arc["dipole_derivatives"], gradient=arc["forces"], charge=charge, multiplicity=mult,
        printed_freqs=pf, printed_ir=pir, printed_masses=None if printed is None else np.array(printed), layers=layers, frozen=frozen or [], route=route,
        program=arc["version"], title=arc["title"], termination=_termination(text), anharmonic_present=anh,
        notes=notes,
    )
