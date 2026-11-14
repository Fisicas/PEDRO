"""Internal coordinates: primitives, analytic Wilson B matrix, automatic non-redundant sets, VEDA .dd2 I/O.

Primitive conventions (all angles in radians, lengths in bohr; atom indices 0-based):
  stretch(a, b)        r = |x_b - x_a|
  bend(a, b, c)        angle a-b-c with b the apex, in [0, pi]
  torsion(a, b, c, d)  dihedral between planes (a, b, c) and (b, c, d), in (-pi, pi]; IUPAC sign
  oop(a, b, c, d)      VEDA 'OUT A B C D': angle between the bond vector D->A and the plane (B, C, D);
                       D is the central atom; arcsin of the unit normal (B-D)x(C-D) dotted with the unit D->A
B-matrix rows are analytic first derivatives (docs/theory.md, section 2). Central finite differences of
`value` are used only in the tests, as an independent check.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .constants import BOHR_TO_ANG, COVALENT_RADII

GREEK = {"stretch": "ν", "bend": "δ", "torsion": "τ", "oop": "γ"}
ASCII = {"stretch": "nu", "bend": "delta", "torsion": "tau", "oop": "gamma"}
NATOMS = {"stretch": 2, "bend": 3, "torsion": 4, "oop": 4}
DASH = "–"


@dataclass(frozen=True)
class Primitive:
    kind: str
    atoms: tuple[int, ...]

    def label(self, elements, ascii_only=False) -> str:
        names = [f"{elements[a]}{a + 1}" for a in self.atoms]
        if ascii_only:
            return f"{ASCII[self.kind]}({'-'.join(names)})"
        return f"{GREEK[self.kind]}({DASH.join(names)})"


@dataclass
class Coordinate:
    label: str
    kind: str
    components: list[tuple[Primitive, float]]
    rule: str = "primitive"
    id: str = ""

    @property
    def atoms(self) -> set[int]:
        return {a for p, _ in self.components for a in p.atoms}

    @property
    def is_cartesian(self) -> bool:
        return self.kind.startswith("cart")


# ----------------------------------------------------------------------------- values
def _unit(v):
    return v / np.linalg.norm(v)


def value(p: Primitive, x) -> float:
    a = p.atoms
    if p.kind == "stretch":
        return float(np.linalg.norm(x[a[1]] - x[a[0]]))
    if p.kind == "bend":
        u, v = x[a[0]] - x[a[1]], x[a[2]] - x[a[1]]
        c = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
        return float(math.acos(max(-1.0, min(1.0, c))))
    if p.kind == "torsion":
        b1, b2, b3 = x[a[1]] - x[a[0]], x[a[2]] - x[a[1]], x[a[3]] - x[a[2]]
        n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
        return float(math.atan2(np.linalg.norm(b2) * np.dot(b1, n2), np.dot(n1, n2)))
    if p.kind == "oop":
        A, B, C, D = (x[i] for i in a)
        n = np.cross(B - D, C - D)
        return float(math.asin(max(-1.0, min(1.0, np.dot(_unit(n), _unit(A - D))))))
    raise ValueError(p.kind)


# ----------------------------------------------------------------------------- analytic B rows
def b_row(p: Primitive, x) -> np.ndarray:
    """d(value)/d(Cartesian), flattened to length 3N (bohr^-0 for stretches, rad/bohr for angles)."""
    row = np.zeros_like(x, dtype=float)
    a = p.atoms
    if p.kind == "stretch":
        e = _unit(x[a[1]] - x[a[0]])
        row[a[0]], row[a[1]] = -e, e
    elif p.kind == "bend":
        u, v = x[a[0]] - x[a[1]], x[a[2]] - x[a[1]]
        lu, lv = np.linalg.norm(u), np.linalg.norm(v)
        eu, ev = u / lu, v / lv
        c = float(np.dot(eu, ev))
        s = math.sqrt(max(1e-300, 1.0 - c * c))
        da = (c * eu - ev) / (lu * s)
        dc = (c * ev - eu) / (lv * s)
        row[a[0]], row[a[2]], row[a[1]] = da, dc, -(da + dc)
    elif p.kind == "torsion":
        # Blondel & Karplus, J. Comput. Chem. 17, 1132 (1996): singularity-free form
        F, G, H = x[a[0]] - x[a[1]], x[a[1]] - x[a[2]], x[a[3]] - x[a[2]]
        A, B = np.cross(F, G), np.cross(H, G)
        A2, B2, Gn = A @ A, B @ B, np.linalg.norm(G)
        fg, hg = F @ G, H @ G
        row[a[0]] = -Gn / A2 * A
        row[a[1]] = Gn / A2 * A + fg / (A2 * Gn) * A - hg / (B2 * Gn) * B
        row[a[2]] = hg / (B2 * Gn) * B - fg / (A2 * Gn) * A - Gn / B2 * B
        row[a[3]] = Gn / B2 * B
    elif p.kind == "oop":
        A, B, C, D = (x[i] for i in a)
        u, v, w = A - D, B - D, C - D
        n = np.cross(v, w)
        nu, nn = np.linalg.norm(u), np.linalg.norm(n)
        s = float(n @ u) / (nn * nu)
        cos = math.sqrt(max(1e-300, 1.0 - s * s))
        ds_du = n / (nn * nu) - s * u / nu**2
        g = u / (nn * nu) - s * n / nn**2           # ds/dn
        row[a[0]] = ds_du / cos
        row[a[1]] = np.cross(w, g) / cos
        row[a[2]] = np.cross(g, v) / cos
        row[a[3]] = -(row[a[0]] + row[a[1]] + row[a[2]])
    else:
        raise ValueError(p.kind)
    return row.reshape(-1)


def b_row_numeric(p: Primitive, x, h=1e-5) -> np.ndarray:
    """Central finite differences of `value` (test oracle only)."""
    flat = x.reshape(-1).astype(float)
    out = np.zeros(flat.size)
    for i in {3 * k + c for k in p.atoms for c in range(3)}:
        xp, xm = flat.copy(), flat.copy()
        xp[i] += h
        xm[i] -= h
        d = value(p, xp.reshape(x.shape)) - value(p, xm.reshape(x.shape))
        if p.kind == "torsion":
            d = (d + math.pi) % (2 * math.pi) - math.pi
        out[i] = d / (2 * h)
    return out


def coordinate_row(c: Coordinate, x) -> np.ndarray:
    if c.is_cartesian:
        row = np.zeros(x.size)
        atom, axis = c.components[0][0].atoms[0], "xyz".index(c.kind[-1])
        row[3 * atom + axis] = 1.0
        return row
    return sum(w * b_row(p, x) for p, w in c.components)


def b_matrix(coords: list[Coordinate], x) -> np.ndarray:
    return np.vstack([coordinate_row(c, x) for c in coords])


R_REF_BOHR = 1.0 / BOHR_TO_ANG   # angular rows (rad/bohr) are scaled by 1 angstrom for conditioning only


def row_scale(c: Coordinate) -> float:
    return 1.0 if c.is_cartesian or all(p.kind == "stretch" for p, _ in c.components) else R_REF_BOHR


def scaled_b_matrix(coords, x):
    return b_matrix(coords, x) * np.array([row_scale(c) for c in coords])[:, None]


def cartesian_coordinates(atoms, elements) -> list[Coordinate]:
    out = []
    for a in atoms:
        for ax in "xyz":
            out.append(Coordinate(f"{ax}({elements[a]}{a + 1})", f"cart_{ax}", [(Primitive("cart", (a,)), 1.0)],
                                  "environment_cartesian"))
    return out


# ----------------------------------------------------------------------------- connectivity
def infer_bonds(elements, xyz_bohr, atoms=None, scale=1.25) -> list[tuple[int, int]]:
    """Covalent bonds among `atoms` (all atoms if None): d <= scale * (r_i + r_j) with Cordero radii."""
    idx = np.arange(len(elements)) if atoms is None else np.array(sorted(atoms))
    xa = xyz_bohr[idx] * BOHR_TO_ANG
    r = np.array([COVALENT_RADII.get(elements[i], 1.5) for i in idx])
    d = np.linalg.norm(xa[:, None, :] - xa[None, :, :], axis=-1)
    lim = scale * (r[:, None] + r[None, :])
    ii, jj = np.nonzero(np.triu(d <= lim, 1))
    return [(int(idx[i]), int(idx[j])) for i, j in zip(ii, jj)]


def _neighbors(n, bonds):
    nb = [set() for _ in range(n)]
    for i, j in bonds:
        nb[i].add(j)
        nb[j].add(i)
    return nb


def fragments(atoms, bonds) -> list[list[int]]:
    """Connected components of the bond graph restricted to `atoms`."""
    atoms = set(atoms)
    adj = {a: set() for a in atoms}
    for i, j in bonds:
        if i in atoms and j in atoms:
            adj[i].add(j)
            adj[j].add(i)
    seen, out = set(), []
    for s in sorted(atoms):
        if s in seen:
            continue
        comp, stack = [], [s]
        seen.add(s)
        while stack:
            a = stack.pop()
            comp.append(a)
            for b in adj[a] - seen:
                seen.add(b)
                stack.append(b)
        out.append(sorted(comp))
    return sorted(out, key=len, reverse=True)


# ----------------------------------------------------------------------------- automatic coordinate set
def primitive_candidates(elements, x, bonds, linear_deg=175.0) -> list[Primitive]:
    nb = _neighbors(len(elements), bonds)
    out = [Primitive("stretch", (i, j)) for i, j in bonds]
    for b in sorted({a for bond in bonds for a in bond}):
        for a, c in itertools.combinations(sorted(nb[b]), 2):
            if math.degrees(value(Primitive("bend", (a, b, c)), x)) >= linear_deg:
                raise NotImplementedError(f"Near-linear angle {elements[a]}{a+1}-{elements[b]}{b+1}-{elements[c]}{c+1}: "
                                          "linear-bend coordinates are not implemented; supply a .dd2 coordinate set")
            out.append(Primitive("bend", (a, b, c)))
    for d in sorted({a for bond in bonds for a in bond}):
        if len(nb[d]) == 3:
            ring = sorted(nb[d])
            for k in range(3):
                a = ring[k]
                b, c = [r for r in ring if r != a]
                out.append(Primitive("oop", (a, b, c, d)))
    for b, c in bonds:
        for a in sorted(nb[b] - {c}):
            for d in sorted(nb[c] - {b}):
                if a != d:
                    out.append(Primitive("torsion", (a, b, c, d)))
    return out


def natural_xh2(elements, bonds) -> list[Coordinate]:
    """Terminal X(H)2 groups on a two-coordinate center: symmetric/antisymmetric stretch and scissor."""
    nb = _neighbors(len(elements), bonds)
    s2 = 1 / math.sqrt(2)
    out = []
    for X in sorted({a for b in bonds for a in b}):
        Hs = sorted(h for h in nb[X] if elements[h] == "H" and len(nb[h]) == 1)
        if len(Hs) == 2 and len(nb[X]) == 2:
            h1, h2 = Hs
            tag = f"{elements[X]}{X+1}{DASH}{elements[h1]}{h1+1},{elements[h2]}{h2+1}"
            r1, r2 = Primitive("stretch", (X, h1)), Primitive("stretch", (X, h2))
            out += [Coordinate(f"νₛ({tag})", "stretch_sym", [(r1, s2), (r2, s2)], "XH2_terminal"),
                    Coordinate(f"νₐₛ({tag})", "stretch_antisym", [(r1, s2), (r2, -s2)], "XH2_terminal"),
                    Coordinate(f"δ({elements[h1]}{h1+1}{DASH}{elements[X]}{X+1}{DASH}{elements[h2]}{h2+1})",
                               "bend", [(Primitive("bend", (h1, X, h2)), 1.0)], "XH2_terminal")]
    return out


def _candidates(elements, x, bonds):
    nat = natural_xh2(elements, bonds)
    covered = {(p.kind, p.atoms) for c in nat for p, _ in c.components}
    order = {"stretch": 0, "bend": 1, "oop": 2, "torsion": 3}
    prims = sorted((p for p in primitive_candidates(elements, x, bonds) if (p.kind, p.atoms) not in covered),
                   key=lambda p: order[p.kind])
    stretches = [Coordinate(p.label(elements), p.kind, [(p, 1.0)], "bond_stretch") for p in prims if p.kind == "stretch"]
    rest = [Coordinate(p.label(elements), p.kind, [(p, 1.0)], "primitive_fill") for p in prims if p.kind != "stretch"]
    return [c for c in nat if c.kind.startswith("stretch")] + stretches + [c for c in nat if c.kind == "bend"] + rest


def select_nonredundant(coords, x, target, tol):
    """Greedy Gram-Schmidt on dimensionally scaled B rows, in chemical priority order."""
    B = scaled_b_matrix(coords, x)
    basis, kept = [], []
    for c, row in zip(coords, B):
        r = row.copy()
        for q in basis:
            r -= (r @ q) * q
        rel = np.linalg.norm(r) / max(np.linalg.norm(row), 1e-300)
        if rel > tol:
            basis.append(r / np.linalg.norm(r))
            kept.append(c)
            if len(kept) == target:
                break
    return kept


TOL_CASCADE = (0.1, 0.05, 0.02, 0.01, 1e-3)


def build(elements, x, atoms=None, bonds=None, linear=False):
    """Automatic non-redundant set for `atoms` (all atoms if None). Returns (coordinates, report).

    Target size is 3n-6 (3n-5 for linear fragments). The strictest tolerance in TOL_CASCADE that still gives a
    complete set is used (the 0.1 level is the most stable under small geometry perturbations; measured on
    a planar C2v test molecule in dev versions)."""
    atoms = list(range(len(elements))) if atoms is None else sorted(atoms)
    if bonds is None:
        bonds = infer_bonds(elements, x, atoms)
    frs = fragments(atoms, bonds)
    if len(frs) > 1:
        raise ValueError(f"The selected atoms form {len(frs)} disconnected fragments (sizes {[len(f) for f in frs]}); "
                         "an automatic internal-coordinate set needs one bonded molecule")
    target = 3 * len(atoms) - (5 if linear else 6)
    cand = _candidates(elements, x, bonds)
    tried = []
    for tol in TOL_CASCADE:
        kept = select_nonredundant(cand, x, target, tol)
        tried.append({"tolerance": tol, "n_selected": len(kept)})
        if len(kept) == target:
            break
    else:
        raise ValueError(f"Could build only {len(kept)} independent internal coordinates; need {target}")
    for k, c in enumerate(kept, 1):
        c.id = f"S{k}"
    return kept, {"source": "automatic", "n_candidates": len(cand), "n_target": target, "tolerance": tol,
                  "cascade": tried, "n_bonds": len(bonds), "set_id": set_id(kept, elements)}


def set_id(coords, elements) -> str:
    canon = [[c.kind, [[p.kind, list(map(int, p.atoms)), round(float(w), 6)] for p, w in c.components]] for c in coords]
    return "cs-" + hashlib.sha256(json.dumps({"e": list(elements), "c": canon}).encode()).hexdigest()[:16]


# ----------------------------------------------------------------------------- VEDA .dd2 and JSON I/O
_DD2_KIND = {"STRE": "stretch", "BEND": "bend", "TORS": "torsion", "OUT": "oop"}
_KIND_DD2 = {v: k for k, v in _DD2_KIND.items()}


def _check(p, n):
    k = NATOMS[p.kind]
    if len(p.atoms) != k or len(set(p.atoms)) != k or not all(0 <= a < n for a in p.atoms):
        raise ValueError(f"invalid atoms {[a + 1 for a in p.atoms]} for {p.kind} with {n} atoms")


def import_dd2(path, elements, atom_map=None) -> list[Coordinate]:
    """VEDA4 .dd2: 's' and 'k' coordinate lines with optional coefficient continuation lines; everything after
    '****' or 'alternative' is ignored. `atom_map` (list: VEDA atom number - 1 -> target 0-based index) applies the
    set to the same molecule embedded in a larger structure."""
    n = len(elements)
    mp = (lambda a: a) if atom_map is None else (lambda a: int(atom_map[a]))
    coords, cur, last = [], None, None
    for line in Path(path).read_text(errors="replace").splitlines():
        if line.startswith("****") or line.startswith("alternative"):
            break
        m = re.match(r"^[sk]\s+(\d+)\s+(-?\d+\.\d+)\s+(STRE|BEND|TORS|OUT)\s+(.*)$", line)
        if m:
            last = _DD2_KIND[m.group(3)]
            atoms = tuple(mp(int(v) - 1) for v in m.group(4).split()[:NATOMS[last]])
            p = Primitive(last, atoms)
            _check(p, n)
            cur = Coordinate("", last, [(p, float(m.group(2)))], "veda_dd2", f"S{m.group(1)}")
            coords.append(cur)
            continue
        if line[:1].isalpha():
            cur = None
            continue
        m = re.match(r"^\s+(-?\d+\.\d+)\s+(?:(STRE|BEND|TORS|OUT)\s+)?(\d+(?:\s+\d+){1,3})\b", line)
        if m and cur is not None:
            if m.group(2):
                last = _DD2_KIND[m.group(2)]
            atoms = tuple(mp(int(v) - 1) for v in m.group(3).split()[:NATOMS[last]])
            p = Primitive(last, atoms)
            _check(p, n)
            cur.components.append((p, float(m.group(1))))
    if not coords:
        raise ValueError(f"no coordinate definitions found in {path}")
    for c in coords:
        c.label = " + ".join((f"{w:+.2f}·" if len(c.components) > 1 else "") + p.label(elements)
                             for p, w in c.components)
    return coords


def export_dd2(coords, elements, x, path=None) -> str:
    """Write a VEDA-readable .dd2 in VEDA's own column layout (as in the files VEDA 4.99 writes): the coordinate type
    only on the first line of a combination, values in angstrom (6 decimals) or degrees (2 decimals), CRLF line
    endings. Each combination is scaled so its largest coefficient is +-1 (TED and PED do not depend on the scale of a
    coordinate, theory section 5), so that VEDA's two-decimal coefficients keep symmetric combinations exact.
    Cartesian rows are skipped."""
    lines = []
    for k, c in enumerate([c for c in coords if not c.is_cartesian], 1):
        s = max(abs(w) for _, w in c.components) or 1.0
        for j, (p, w) in enumerate(c.components):
            w = w / s
            atoms = "".join(f"{a + 1:5d}" for a in p.atoms)
            v = value(p, x)
            v = f"{v * BOHR_TO_ANG:.6f}" if p.kind == "stretch" else f"{math.degrees(v):.2f}"
            tag = "".join(elements[a] for a in p.atoms)
            head = f"s {k:<4d}" if j == 0 else " " * 6
            kind = f"   {_KIND_DD2[p.kind]:<4s}" if j == 0 else " " * 7
            lines.append(f"{head}{w:6.2f}{kind}{atoms}   {tag:<6s}{v}")
    text = "\r\n".join(lines) + "\r\n****\r\n"
    if path is not None:
        Path(path).write_bytes(text.encode("ascii"))
    return text


FORMAT = "pedro-coordinate-set/1"


def export_json(coords, elements) -> dict:
    return {"format": FORMAT, "set_id": set_id(coords, elements), "elements": list(elements),
            "coordinates": [{"id": c.id, "label": c.label, "kind": c.kind, "rule": c.rule,
                             "components": [{"primitive": p.kind, "atoms_1based": [a + 1 for a in p.atoms],
                                             "coefficient": w} for p, w in c.components]}
                            for c in coords if not c.is_cartesian]}


def import_json(path, elements) -> list[Coordinate]:
    d = json.loads(Path(path).read_text(encoding="utf-8"))
    if d.get("format") != FORMAT:
        raise ValueError(f"{path}: not a {FORMAT} file")
    out = []
    for k, c in enumerate(d["coordinates"], 1):
        comps = []
        for comp in c["components"]:
            p = Primitive(comp["primitive"], tuple(a - 1 for a in comp["atoms_1based"]))
            _check(p, len(elements))
            comps.append((p, float(comp["coefficient"])))
        out.append(Coordinate(c["label"], c["kind"], comps, c.get("rule", "imported"), f"S{k}"))
    return out


def load_set(path, elements, atom_map=None) -> list[Coordinate]:
    if str(path).lower().endswith(".dd2"):
        return import_dd2(path, elements, atom_map)
    return import_json(path, elements)


# ----------------------------------------------------------------------------- classification for summaries
def coordinate_class(c: Coordinate, elements) -> str:
    """Chemical class used for group totals, e.g. 'ν(C–O)', 'ν(C–H)', 'δ', 'τ', 'γ', 'environment'."""
    if c.is_cartesian:
        return "environment"
    kinds = {p.kind for p, _ in c.components}
    if kinds == {"stretch"}:
        pairs = {DASH.join(sorted((elements[p.atoms[0]], elements[p.atoms[1]]))) for p, _ in c.components}
        return f"ν({'/'.join(sorted(pairs))})"
    return {"bend": "δ (bend)", "torsion": "τ (torsion)", "oop": "γ (out-of-plane)"}.get(next(iter(kinds)), "mixed") \
        if len(kinds) == 1 else "mixed"
