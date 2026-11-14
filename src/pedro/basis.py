"""Square coordinate bases for isolated molecules and for a subsystem embedded in a larger Hessian.

Cases (docs/theory.md, section 8):
  isolated  internal coordinates of the whole molecule                               3N-6 coordinates, projected modes
  F         internals of subsystem S + Cartesians of every environment atom E         3N-6, projected modes (default)
  Fprime    internals of S + 6 rigid-body coordinates of S + Cartesians of E         3N, unprojected modes
  P         partial Hessian, mobile set M (S, or S plus part of E): internals of S + 6 rigid-body of S
            + Cartesians of the mobile environment atoms                             3N_M, partial modes
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import coordinates as co


@dataclass
class Basis:
    case: str
    B: np.ndarray                    # (n_coords, 3N)
    labels: list[str]
    groups: list[str]                # group name per coordinate (for group totals)
    internal: list[co.Coordinate]    # internal coordinates of S, in basis order (first len(internal) rows)
    subsystem: np.ndarray            # 0-based atoms of S
    environment: np.ndarray          # 0-based atoms of E
    info: dict = field(default_factory=dict)

    @property
    def n_internal(self):
        return len(self.internal)

    def group_index(self):
        out: dict[str, list[int]] = {}
        for i, g in enumerate(self.groups):
            out.setdefault(g, []).append(i)
        return {g: np.array(v) for g, v in out.items()}


def atom_ranges(atoms) -> str:
    """[0,1,2,5,7,8] -> '1-3,6,8-9' (1-based, as in Gaussian)."""
    a = sorted(int(x) + 1 for x in atoms)
    out, start = [], None
    for k, v in enumerate(a):
        if start is None:
            start = v
        if k + 1 == len(a) or a[k + 1] != v + 1:
            out.append(f"{start}" if start == v else f"{start}-{v}")
            start = None
    return ",".join(out)


def _formula(atoms, elements):
    c = {}
    for a in atoms:
        c[elements[a]] = c.get(elements[a], 0) + 1
    keys = ["C", "H"] + sorted(k for k in c if k not in ("C", "H")) if "C" in c else sorted(c)
    return "".join(f"{k}{c[k] if c[k] > 1 else ''}" for k in keys if k in c)


def _peptide_like(atoms, elements, bonds):
    """True if the fragment contains a peptide backbone motif N-C-C(=O): a protein piece, not a cofactor."""
    s = set(atoms)
    nb = {a: set() for a in s}
    for i, j in bonds:
        if i in s and j in s:
            nb[i].add(j)
            nb[j].add(i)
    for n in (a for a in s if elements[a] == "N"):
        for ca in (x for x in nb[n] if elements[x] == "C"):
            for c in (x for x in nb[ca] if elements[x] == "C" and x != n):
                if any(elements[o] == "O" and len(nb[o]) == 1 for o in nb[c]):
                    return True
    return False


def describe_fragments(job, atoms=None):
    """Bonded fragments (largest first) with formula, Gaussian atom numbers and a protein/water tag."""
    atoms = list(range(job.natoms)) if atoms is None else list(atoms)
    b = co.infer_bonds(job.elements, job.xyz_bohr, atoms)
    out = []
    for f in co.fragments(atoms, b):
        form = _formula(f, job.elements)
        tag = "water" if form == "H2O" else "peptide" if _peptide_like(f, job.elements, b) else ""
        out.append({"atoms": f, "n": len(f), "formula": form, "ranges": atom_ranges(f), "tag": tag})
    return out


def auto_subsystem(job, bonds=None):
    """Subsystem chosen without user input: the largest bonded fragment that is not water and has no peptide
    backbone (N-C-C=O), taken from the ONIOM high layer when layers are known. Returns (atoms, explanation,
    fragment list); the explanation and the list are printed so the choice can be overridden."""
    n = job.natoms
    scope = [i for i, l in enumerate(job.layers) if l == "H"] if job.layers else list(range(n))
    where = "ONIOM high layer" if job.layers else "structure"
    frs = describe_fragments(job, scope)
    if len(frs) == 1 and not job.layers:
        return np.arange(n), "whole structure (one bonded molecule)", frs
    pick = next((f for f in frs if not f["tag"]), frs[0])
    skipped = [f for f in frs[:frs.index(pick)]]
    why = (f"largest bonded fragment of the {where} that is not protein or water: {pick['formula']}, atoms "
           f"{pick['ranges']}" + (f" (skipped {', '.join(f['formula'] + ' [' + f['tag'] + ']' for f in skipped)})"
                                  if skipped else ""))
    return np.array(pick["atoms"]), why, frs


def rigid_rows(atoms, xyz, masses, n3):
    """B rows of the 3 center-of-mass translations and 3 linearized rotations (about the principal axes) of `atoms`.
    Rotation rows are scaled by 1 angstrom so all rows have length units; TED does not depend on row scaling."""
    atoms = np.asarray(atoms)
    m = masses[atoms]
    M = m.sum()
    r = xyz[atoms] - np.average(xyz[atoms], axis=0, weights=m)
    I = np.einsum("a,ab,ac->bc", m, r, r)
    I = np.trace(I) * np.eye(3) - I
    w, axes = np.linalg.eigh(I)
    rows = []
    for e in np.eye(3):
        row = np.zeros((n3 // 3, 3))
        row[atoms] = (m / M)[:, None] * e[None, :]
        rows.append(row.reshape(-1))
    Iinv = np.linalg.pinv(I)
    for k in range(3):
        ax = axes[:, k]
        # theta_k = ax . I^-1 sum_a m_a r_a x dx_a  ->  d theta_k / d dx_a = m_a (I^-1 ax) x r_a ... written as a row:
        v = Iinv @ ax
        row = np.zeros((n3 // 3, 3))
        row[atoms] = m[:, None] * np.cross(r, v[None, :]) * (-1.0)
        rows.append(row.reshape(-1) * co.R_REF_BOHR)
    return np.vstack(rows), ["T_x(S)", "T_y(S)", "T_z(S)", "R_1(S)", "R_2(S)", "R_3(S)"]


def _env_group(a, job):
    if job.layers:
        return {"H": "environment: QM layer", "M": "environment: middle layer", "L": "environment: MM layer"}.get(
            job.layers[a], "environment")
    return "environment"


def build(job, case, internal, subsystem, mobile=None, user_groups=None):
    """Assemble the basis. `internal` are the internal coordinates of S (global atom indices)."""
    n = job.natoms
    n3 = 3 * n
    S = np.array(sorted(subsystem))
    E = np.array(sorted(set(range(n)) - set(S.tolist())))
    el = job.elements
    rows = [co.coordinate_row(c, job.xyz_bohr) for c in internal]
    labels = [c.label for c in internal]
    groups = [co.coordinate_class(c, el) for c in internal]
    info = {}
    if case in ("Fprime", "P"):
        R, names = rigid_rows(S, job.xyz_bohr, job.masses, n3)
        rows += list(R)
        labels += names
        groups += ["subsystem rigid-body"] * 6
    if case == "isolated":
        if len(E):
            raise ValueError("isolated case requires the subsystem to be the whole structure")
    elif case in ("F", "Fprime"):
        env = E
    elif case == "P":
        mob = np.array(sorted(set(range(n)) if mobile is None else set(mobile)))
        if not set(S.tolist()) <= set(mob.tolist()):
            raise ValueError("the mobile set must contain the whole subsystem")
        env = np.array(sorted(set(mob.tolist()) - set(S.tolist())))
        info["mobile"] = mob.tolist()
    else:
        raise ValueError(f"unknown case {case!r}")
    if case != "isolated":
        for a in env:
            for k, ax in enumerate("xyz"):
                r = np.zeros(n3)
                r[3 * a + k] = 1.0
                rows.append(r)
                labels.append(f"{ax}({el[a]}{a + 1})")
                g = user_groups.get(int(a)) if user_groups else None
                groups.append(g or _env_group(a, job))
    return Basis(case=case, B=np.vstack(rows), labels=labels, groups=groups, internal=list(internal),
                 subsystem=S, environment=E, info=info)
