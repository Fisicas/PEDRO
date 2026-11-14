"""Symmetry gauge for GVPT2 output (dev20): separate numerical from physical anharmonic couplings using the molecule's
point group.

A force constant or resonance coupling whose symmetry product is not totally symmetric must vanish exactly. Gaussian
run with NoSymm prints them anyway; their size is therefore a direct measure of the numerical noise, and a resonance that
symmetry forbids is an artifact (in one C2v test molecule, a 13.8 cm-1 forbidden coupling between the A1 and B2 C=O
stretches, 0.24 cm-1 apart, produced a spurious 57 % "mixing").

Scope: abelian point groups (C1, Cs, Ci, C2, C2v, C2h, D2, D2h). dev21: candidate symmetry elements come from the
principal axes and from the atoms themselves, so near-symmetric tops (water dimer) work.
That covers planar model compounds (Cs/C2v). Molecules with degenerate principal moments
(symmetric tops) or no symmetry (C1, e.g. molecules with flexible side chains) get status 'unavailable' or a C1
result: the gauge then says nothing and the audit alone applies.

Provides, from one GVPT2 log step plus the checkpoint's normal modes:
  point group and per-mode irreps (characters from the Cartesian normal modes);
  noise floor: statistics of printed symmetry-forbidden cubic/quartic constants (reduced, cm-1) and per-mode maxima;
  classification of Gaussian's Fermi / Darling-Dennison lists and variational couplings: 'symmetry-forbidden',
    'involves corrupted mode' (audit ANOMALOUS_SHIFT/QUARTIC_OUTLIER) or 'allowed, clean';
  symmetry-filtered GVPT2: Gaussian's printed variational matrix re-diagonalized keeping only allowed, clean couplings
    (a full re-diagonalization must reproduce Gaussian's final fundamentals: parse check);
  contamination per fundamental: 1/2 sum_j x_ij over corrupted partners j (the part of the shift inherited from
    corrupted modes; removing it gives a bound, not a correction).
Mode numbers are Gaussian's anharmonic numbering (as in pedro.gvpt2_quality).
"""
from __future__ import annotations

import re

import numpy as np

from . import quality as gvpt2_quality

BOHR_TO_ANG = 0.529177210903
MAP_TOL_ANG = 0.02          # atom-mapping tolerance for a symmetry operation
SIZE_LIMIT = 200            # atoms; above this the gauge reports unavailable at once (spec D10, overridable)
PARSE_TOL_CM = 0.1          # full re-diagonalization vs Gaussian final fundamentals
SEVERE_KINDS = ("ANOMALOUS_SHIFT", "QUARTIC_OUTLIER")
GROUP_NAMES = {1: "C1"}


# ----------------------------------------------------------------------------- point group
def _equivalence_classes(Z, r, tol):
    """Atoms that a symmetry operation could exchange: same element and same distance from the center of mass."""
    d = np.linalg.norm(r, axis=1)
    classes = {}
    for i in np.argsort(d):
        key = int(Z[i])
        for (z, d0), members in classes.items():
            if z == key and abs(d[i] - d0) < tol:
                members.append(int(i))
                break
        else:
            classes[(key, float(d[i]))] = [int(i)]
    return list(classes.values())


def _candidate_directions(Z, r, axes, tol=MAP_TOL_ANG, limit=20000):
    """Unit vectors that can be a mirror normal or a C2 axis.

    Any off-element mirror swaps pairs of equivalent atoms i, j, so its normal is along r_i - r_j (or it is the plane
    containing all atoms, i.e. a principal axis). A C2 axis either contains an atom (along r_i) or swaps an equivalent
    pair (along r_i + r_j, or perpendicular to r_i and to the other pairs when r_i = -r_j, covered by r_i x r_k).
    Candidates therefore come from principal axes, atom directions and pairs WITHIN equivalence classes (element and
    distance from the center of mass), and are de-duplicated on a hashed direction grid. This replaces dev21's
    all-pairs search (about O(N^4) with the linear de-duplication), which made a 682-atom model effectively hang."""
    seen, out = set(), []

    def add(v):
        nv = np.linalg.norm(v)
        if nv < 1e-6:
            return
        v = v / nv
        k = int(np.argmax(np.abs(v)))
        if v[k] < 0:
            v = -v
        key = tuple(np.round(v * 2000).astype(int))       # ~0.03 degree grid
        if key not in seen:
            seen.add(key)
            out.append(v)
    for a in axes.T:
        add(a)
    for i in range(len(r)):
        add(r[i])
    for cls in _equivalence_classes(Z, r, tol):
        for a in range(len(cls)):
            for b in range(a + 1, len(cls)):
                i, j = cls[a], cls[b]
                add(r[i] - r[j])
                add(r[i] + r[j])
                if np.linalg.norm(r[i] + r[j]) < tol:
                    for k in cls:
                        if k not in (i, j):
                            add(np.cross(r[i], r[k]))
                if len(out) > limit:
                    return out
    return out


class _Mapper:
    """Fast test of whether an orthogonal transformation maps the structure onto itself.

    A few probe atoms (the heaviest-element atoms farthest from the center of mass) are checked first, so most
    candidate operations are rejected in microseconds; survivors are checked atom by atom with per-element k-d trees.
    dev21 checked every candidate against every atom pair (O(N^2) per candidate)."""

    def __init__(self, Z, r, n_probe=6):
        from scipy.spatial import cKDTree
        self.Z, self.r = np.asarray(Z), np.asarray(r)
        self.groups = {z: np.where(self.Z == z)[0] for z in np.unique(self.Z)}
        self.trees = {z: cKDTree(self.r[idx]) for z, idx in self.groups.items()}
        d = np.linalg.norm(self.r, axis=1)
        order = np.lexsort((-d, -self.Z))
        self.probe = order[:n_probe]

    def __call__(self, R, tol):
        rr = self.r @ R.T
        for i in self.probe:
            z = self.Z[i]
            dist, _ = self.trees[z].query(rr[i])
            if dist >= tol:
                return None
        perm = np.empty(len(self.r), int)
        err = 0.0
        for z, idx in self.groups.items():
            dist, j = self.trees[z].query(rr[idx])
            if dist.max() >= tol:
                return None
            perm[idx] = idx[j]
            err = max(err, float(dist.max()))
        if len(np.unique(perm)) != len(perm):
            return None
        return perm.tolist(), err


def point_group(atomic_numbers, masses, xyz_ang, tol=MAP_TOL_ANG) -> dict:
    """Abelian point group (C1, Cs, Ci, C2, C2v, C2h, D2, D2h). Returns {'status', 'group', 'operations': {name: perm},
    'rotations': {name: 3x3}, 'max_mapping_error_angstrom', 'planar'}. Groups with more than 8 operations or with
    C3 or higher axes (non-abelian / degenerate irreps) are reported 'unavailable'."""
    Z, m, x = np.asarray(atomic_numbers), np.asarray(masses, float), np.asarray(xyz_ang, float)
    com = (m[:, None] * x).sum(0) / m.sum()
    r = x - com
    inertia = sum(mi * (ri @ ri * np.eye(3) - np.outer(ri, ri)) for mi, ri in zip(m, r))
    w, axes = np.linalg.eigh(inertia)
    normal = axes[:, 2]
    planar = bool(len(r) >= 3 and np.all(np.abs(r @ normal) < tol))
    found = []                                             # (kind, direction, R, perm, err)
    _maps = _Mapper(Z, r)
    inv = _maps(-np.eye(3), tol)
    if inv:
        found.append(("i", None, -np.eye(3), *inv))
    for v in _candidate_directions(Z, r, axes, tol):
        P = np.outer(v, v)
        for kind, R in (("sigma", np.eye(3) - 2 * P), ("C2", 2 * P - np.eye(3))):
            if any(np.allclose(R, f[2], atol=1e-6) for f in found):
                continue
            res = _maps(R, tol)
            if res:
                found.append((kind, v, R, *res))
    if len(found) + 1 > 8:
        return {"status": "unavailable", "group": None, "reason": f"{len(found) + 1} symmetry operations found: "
                "non-abelian or degenerate point group (not supported)"}
    # name operations deterministically; for planar molecules the molecular plane is 'sigma(mol)'
    ops, rots, worst = {"E": list(range(len(r)))}, {"E": np.eye(3).tolist()}, 0.0
    sig = sorted([f for f in found if f[0] == "sigma"], key=lambda f: tuple(np.round(np.abs(f[1]), 6)))
    c2s = sorted([f for f in found if f[0] == "C2"], key=lambda f: tuple(np.round(np.abs(f[1]), 6)))
    k = 0
    for f in sig:
        if planar and abs(f[1] @ normal) > 1 - 1e-4:
            name = "sigma(mol)"
        else:
            k += 1
            name = f"sigma{k}"
        ops[name], rots[name], worst = f[3], f[2].tolist(), max(worst, f[4])
    for n_, f in enumerate(c2s, 1):
        ops[f"C2_{n_}"], rots[f"C2_{n_}"], worst = f[3], f[2].tolist(), max(worst, f[4])
    if inv:
        ops["i"], rots["i"], worst = inv[0], (-np.eye(3)).tolist(), max(worst, inv[1])
    n = len(ops)
    names = list(ops)
    planes = [s for s in names if s.startswith("sigma")]
    nrot = sum(s.startswith("C2") for s in names)
    group = {1: "C1", 8: "D2h"}.get(n)
    if n == 2:
        group = "Cs" if planes else "Ci" if "i" in names else "C2"
    elif n == 4:
        group = "C2h" if "i" in names else ("D2" if nrot == 3 else "C2v")
    elif n not in (1, 2, 4, 8):
        return {"status": "unavailable", "group": None, "reason": f"{n} operations: not an abelian point group"}
    return {"status": "ok", "group": group, "operations": ops, "rotations": rots, "planar": planar,
            "max_mapping_error_angstrom": round(worst, 6), "principal_moments_amu_ang2": w.tolist()}


def mode_characters(pg: dict, modes_cart: np.ndarray) -> list[dict]:
    """Characters of each Cartesian normal mode (n_modes x N x 3) under the group's operations (identity omitted)."""
    names = [k for k in pg["operations"] if k != "E"]
    out = []
    for d in modes_cart:
        chars = {}
        for k in names:
            R, perm = np.asarray(pg["rotations"][k]), pg["operations"][k]
            dd = np.zeros_like(d)
            for i, j in enumerate(perm):
                dd[j] = R @ d[i]
            chars[k] = float((dd * d).sum() / (d * d).sum())
        out.append(chars)
    return out


def irrep_label(group: str, chars: dict) -> str:
    """Mulliken label for Cs and C2v (planar C2v: B2 = in-plane antisymmetric, i.e. symmetric under sigma(mol); B1 out of
    plane; non-planar C2v: B1 symmetric under 'sigma1'); otherwise the sign pattern over the operations."""
    s = {k: int(np.sign(v)) for k, v in chars.items()}
    if group == "Cs":
        return "A'" if next(v for k, v in s.items() if k.startswith("sigma")) > 0 else "A''"
    if group == "C2v":
        c2 = next(v for k, v in s.items() if k.startswith("C2"))
        sig = sorted(k for k in s if k.startswith("sigma"))
        if c2 > 0:
            return "A1" if s[sig[0]] > 0 else "A2"
        if "sigma(mol)" in s:
            return "B2" if s["sigma(mol)"] > 0 else "B1"
        return "B1" if s["sigma1"] > 0 else "B2"
    return "Γ(" + ",".join(f"{k}{'+' if v > 0 else '-'}" for k, v in sorted(s.items())) + ")"


# ----------------------------------------------------------------------------- log tables
def force_constants(text, title, nidx):
    if title not in text:
        return None
    i = text.index(title)
    e = text.find("Num. of", i)
    out = []
    for line in text[i:e if e > 0 else i + 3_000_000].splitlines():
        f = line.split()
        if len(f) == nidx + 3 and all(x.isdigit() for x in f[:nidx]):
            out.append((tuple(int(x) for x in f[:nidx]), float(f[nidx])))
    return out


def resonance_lists(text):
    if "Resonance Analysis" not in text:
        return [], [], []
    i = text.index("Resonance Analysis")
    e = text.find("Anharmonic X Matrix", i)
    fermi, dd22, dd11, sec = [], [], [], None
    for line in text[i:e if e > 0 else i + 200000].splitlines():
        s = line.strip()
        if s.startswith("Fermi resonances"):
            sec = "f"; continue
        if s.startswith("(2-2)"):
            sec = "22"; continue
        if s.startswith("(1-1)"):
            sec = "11"; continue
        f = s.split()
        if sec == "f" and len(f) == 6 and all(x.isdigit() for x in f[:3]):
            fermi.append(((int(f[0]), int(f[1]), int(f[2])), float(f[3]), float(f[4]), float(f[5])))
        elif sec == "22" and len(f) == 6 and all(x.isdigit() for x in f[:4]):
            dd22.append((tuple(int(x) for x in f[:4]), float(f[4]), float(f[5])))
        elif sec == "11" and len(f) == 4 and all(x.isdigit() for x in f[:2]):
            dd11.append(((int(f[0]), int(f[1])), float(f[2]), float(f[3])))
    return fermi, dd22, dd11


def _state_key(s):
    return tuple(sorted((int(a), int(b)) for a, b in re.findall(r"(\d+)\((\d)\)", s)))


def variational_blocks(text):
    i = text.find("Variational Matrix Definition")
    e = text.find("Projection of DVPT2 states", i)
    if i < 0 or e < 0:
        return []
    lines = text[i:e].splitlines()
    out = []
    for k, line in enumerate(lines):
        m = re.match(r"^\s+(1-1|1-2|2-2|1-3)\s+\|\s+(.+?)\s+\|\s+\|\s+(.+?)\s*$", line)
        if m and k + 1 < len(lines):
            v = re.match(r"^\s+\|\s+(\S+)\s+\|\s+(\S+)\s+\|\s+(\S+)\s*$", lines[k + 1])
            if v:
                f = [float(x.replace("D", "E")) for x in v.groups()]
                out.append({"type": m.group(1), "s1": _state_key(m.group(2)), "s2": _state_key(m.group(3)),
                            "e1": f[0], "w": f[1], "e2": f[2]})
    return out


def _diagonalise(blocks, keep):
    states = sorted({b["s1"] for b in blocks} | {b["s2"] for b in blocks})
    idx = {s: n for n, s in enumerate(states)}
    H = np.zeros((len(states), len(states)))
    for b in blocks:
        H[idx[b["s1"]], idx[b["s1"]]] = b["e1"]
        H[idx[b["s2"]], idx[b["s2"]]] = b["e2"]
        if keep(b):
            H[idx[b["s1"]], idx[b["s2"]]] = H[idx[b["s2"]], idx[b["s1"]]] = b["w"]
    vals, vecs = np.linalg.eigh(H)
    # One-to-one labels (maximum total squared overlap), as Gaussian assigns them. Labeling each eigenvector by its
    # largest component alone gives two eigenvalues the same label when a resonance mixes states strongly (seen for
    # benzaldehyde's Darling-Dennison 1-1 couplings, 51 %/36 %), and the second silently overwrote the first.
    from scipy.optimize import linear_sum_assignment
    rows, cols = linear_sum_assignment(-(vecs ** 2))      # rows: basis states, cols: eigenvectors
    return {states[int(r)]: float(vals[int(c)]) for r, c in zip(rows, cols)}


# ----------------------------------------------------------------------------- main entry
def gauge(text: str, quality: dict, atomic_numbers, masses, xyz_bohr, mode_vectors: dict, size_limit=SIZE_LIMIT) -> dict:
    """Symmetry gauge of one GVPT2 log step.

    quality      = quality.audit(text)
    mode_vectors = {Gaussian anharmonic mode number: Cartesian displacement (N, 3)} from PEDRO's own harmonic
                   analysis of the same Hessian and geometry (no checkpoint needed)."""
    if not quality or quality.get("status") != "audited":
        return {"status": "unavailable", "reason": "no audited GVPT2 step"}
    if len(atomic_numbers) > size_limit:
        return {"status": "unavailable", "reason": f"size guard: {len(atomic_numbers)} atoms > {size_limit}"}
    pg = point_group(atomic_numbers, masses, np.asarray(xyz_bohr) * BOHR_TO_ANG)
    if pg["status"] != "ok":
        return {"status": "unavailable", "reason": pg["reason"]}
    if pg["group"] == "C1":
        return {"status": "no_symmetry", "group": "C1",
                "reason": "C1: no symmetry-forbidden quantities; the symmetry gauge cannot separate noise from physics"}
    qmodes = {int(m): d for m, d in quality["modes"].items()}
    if not set(qmodes) <= set(mode_vectors):
        return {"status": "unavailable", "reason": "not every Gaussian anharmonic mode could be mapped to a harmonic mode"}
    names = sorted(mode_vectors)
    chars_list = mode_characters(pg, np.array([mode_vectors[m] for m in names]))
    chars_h = dict(zip(names, chars_list))
    char = {m: chars_h[m] for m in qmodes}
    max_dev = max(abs(abs(v) - 1) for c in char.values() for v in c.values())
    if max_dev > 0.05:
        return {"status": "unavailable", "group": pg["group"],
                "reason": f"normal modes are not symmetry-adapted (character deviation {max_dev:.3f}); geometry or modes break the group"}
    sign = {m: np.array([np.sign(char[m][k]) for k in sorted(char[m])]) for m in qmodes}
    irrep = {m: irrep_label(pg["group"], char[m]) for m in qmodes}

    def allowed(idx):
        prod = np.ones(len(next(iter(sign.values()))))
        for i in idx:
            prod = prod * sign[i]
        return bool(np.all(prod > 0))
    severe = {m for m, d in qmodes.items() if any(f.startswith(SEVERE_KINDS) for f in d["flags"])}

    def classify(idx):
        if not allowed(idx):
            return "symmetry-forbidden"
        return "involves corrupted mode" if set(idx) & severe else "allowed, clean"
    # noise floor
    floor, per_mode = {}, {m: 0.0 for m in qmodes}
    for kind, title, n in (("cubic", "CUBIC FORCE CONSTANTS IN NORMAL MODES", 3), ("quartic", "QUARTIC FORCE CONSTANTS IN NORMAL MODES", 4)):
        fc = force_constants(text, title, n)
        if fc is None:
            floor[kind] = None
            continue
        fc = [(i, v) for i, v in fc if all(x in qmodes for x in i)]
        al = np.array([abs(v) for i, v in fc if allowed(i)])
        fo = [(i, v) for i, v in fc if not allowed(i)]
        fa = np.array([abs(v) for _, v in fo])
        for i, v in fo:
            for x in set(i):
                per_mode[x] = max(per_mode[x], abs(v))
        floor[kind] = {"n_allowed": int(len(al)), "n_forbidden": len(fo),
                       "allowed_abs_median": float(np.median(al)) if len(al) else None,
                       "forbidden_abs_p95": float(np.percentile(fa, 95)) if len(fa) else 0.0,
                       "forbidden_abs_max": float(fa.max()) if len(fa) else 0.0,
                       "largest_forbidden": [{"indices": list(i), "value_cm-1": v} for i, v in sorted(fo, key=lambda t: -abs(t[1]))[:5]]}
    # resonances
    fermi, dd22, dd11 = resonance_lists(text)
    res = ([{"type": "Fermi 1-2", "modes": list(i), "freq_diff": d, "reduced_cubic": c, "pt2_minus_variational": p, "class": classify(i)} for i, d, c, p in fermi]
           + [{"type": "DD 2-2", "modes": list(i), "freq_diff": d, "coupling": c, "class": classify(i)} for i, d, c in dd22]
           + [{"type": "DD 1-1", "modes": list(i), "freq_diff": d, "coupling": c, "class": classify(i)} for i, d, c in dd11])
    # variational matrix
    blocks = variational_blocks(text)
    for b in blocks:
        b["class"] = classify(tuple(x for s in (b["s1"], b["s2"]) for x, nq in s for _ in range(nq)))
    full = _diagonalise(blocks, lambda b: True) if blocks else {}
    filt = _diagonalise(blocks, lambda b: b["class"] == "allowed, clean") if blocks else {}
    parse = [abs(full[((m, 1),)] - qmodes[m]["nu_final"]) for m in qmodes if ((m, 1),) in full]
    parse_max = max(parse) if parse else None
    filtered_ok = parse_max is None or parse_max <= PARSE_TOL_CM
    # contamination from the printed X matrix
    X = gvpt2_quality._tri(text.splitlines(), "Total Anharmonic X Matrix", max(qmodes))
    mres = {}
    for m, d in qmodes.items():
        forb = [b for b in blocks if b["class"] == "symmetry-forbidden" and ((m, 1),) in (b["s1"], b["s2"])]
        nu_f = filt.get(((m, 1),), d["nu_depert"]) if (blocks and filtered_ok) else None
        mres[str(m)] = {"irrep": irrep[m], "max_forbidden_fc_cm-1": round(per_mode[m], 3),
                        "nu_symmetry_filtered": None if nu_f is None else round(nu_f, 3),
                        "contamination_cm-1": d.get("contamination_cm-1"),     # dev21: computed by the audit
                        "forbidden_couplings": [{"type": b["type"], "partner": [list(x) for x in (b["s2"] if b["s1"] == ((m, 1),) else b["s1"])],
                                                 "coupling_cm-1": b["w"]} for b in forb]}
    counts = {}
    for r in res:
        counts.setdefault(r["type"], {}).setdefault(r["class"], 0)
        counts[r["type"]][r["class"]] += 1
    return {"status": "ok", "group": pg["group"], "max_mapping_error_angstrom": pg["max_mapping_error_angstrom"],
            "mode_numbering": "Gaussian anharmonic numbering", "irrep_counts": {k: list(irrep.values()).count(k) for k in sorted(set(irrep.values()))},
            "noise_floor": floor, "resonances": res, "resonance_classes": counts,
            "variational_couplings": [{"type": b["type"], "state1": [list(x) for x in b["s1"]], "state2": [list(x) for x in b["s2"]],
                                       "coupling_cm-1": b["w"], "class": b["class"]} for b in blocks],
            "variational_parse_check_max_abs_cm-1": None if parse_max is None else round(parse_max, 4),
            "symmetry_filtered_available": bool(blocks) and filtered_ok,
            "severe_modes": sorted(severe), "modes": mres,
            "notes": "symmetry-filtered = Gaussian's printed variational matrix re-diagonalized with only symmetry-allowed couplings "
                     "among states without corrupted modes; contamination = 1/2 sum_j x_ij over corrupted partners (a bound, "
                     "not a correction). Record positions remain Gaussian's values."}


def state_flags(g: dict, quanta_gaussian: dict) -> list[str]:
    """FORBIDDEN_COUPLING(...) for a fundamental that Gaussian coupled to a state symmetry forbids."""
    if not g or g.get("status") != "ok":
        return []
    qs = {int(k): int(v) for k, v in quanta_gaussian.items()}
    if not (len(qs) == 1 and next(iter(qs.values())) == 1):
        return []
    m = g["modes"].get(str(next(iter(qs))), {})
    fc = m.get("forbidden_couplings", [])
    if not fc:
        return []
    return ["FORBIDDEN_COUPLING(" + ",".join(f"{c['type']}:" + "+".join(f"{a}({b})" for a, b in c["partner"]) + f"={c['coupling_cm-1']:+.1f}"
                                            for c in fc) + ")"]
