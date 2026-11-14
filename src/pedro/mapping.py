"""Transfer a coordinate set defined on one structure onto the same molecule inside another (spec section 6.3).

Typical use: VEDA generated a .dd2 for an isolated molecule (its .fmu); the same molecule sits inside a large ONIOM
model. PEDRO finds every graph isomorphism between the two bonded molecules (same elements, same bonds), keeps the
one whose atoms superpose with the lowest all-atom RMSD (this fixes the arbitrary order of methyl and CH2 hydrogens),
and renumbers the coordinate set. Nothing is reassigned silently: if the molecules are not isomorphic the run stops
and says why.
"""

from __future__ import annotations

import numpy as np

from .constants import BOHR_TO_ANG
from . import coordinates as co

MAX_ISOMORPHISMS = 200000


class MappingError(Exception):
    pass


def _adjacency(n, bonds):
    adj = [set() for _ in range(n)]
    for i, j in bonds:
        adj[i].add(j)
        adj[j].add(i)
    return adj


def _refine(colors_a, adj_a, colors_b, adj_b):
    """Joint Weisfeiler-Lehman color refinement of two graphs (colors comparable across graphs)."""
    ca, cb = list(colors_a), list(colors_b)
    for _ in range(len(ca) + 1):
        table = {}
        sig_a = [(ca[i], tuple(sorted(ca[j] for j in adj_a[i]))) for i in range(len(ca))]
        sig_b = [(cb[i], tuple(sorted(cb[j] for j in adj_b[i]))) for i in range(len(cb))]
        for s in sorted(set(sig_a) | set(sig_b), key=repr):
            table[s] = len(table)
        na, nb = [table[s] for s in sig_a], [table[s] for s in sig_b]
        if len(set(na)) == len(set(ca)) and len(set(nb)) == len(set(cb)):
            return na, nb
        ca, cb = na, nb
    return ca, cb


def kabsch_rmsd(P, Q):
    """Rotation R and RMSD after optimal superposition of P onto Q (rows are points)."""
    p, q = P - P.mean(0), Q - Q.mean(0)
    U, _, Vt = np.linalg.svd(p.T @ q)
    D = np.diag([1.0, 1.0, np.sign(np.linalg.det(U @ Vt))])
    R = U @ D @ Vt
    return R, float(np.sqrt(((p @ R - q) ** 2).sum(1).mean()))


def isomorphisms(el_a, adj_a, el_b, adj_b, limit=MAX_ISOMORPHISMS):
    """Yield every element-preserving graph isomorphism a -> b as a list (index in a -> index in b)."""
    n = len(el_a)
    ca, cb = _refine(el_a, adj_a, el_b, adj_b)
    if sorted(ca) != sorted(cb):
        return
    # order source atoms: rarest color first, then breadth-first so each new atom has mapped neighbors
    freq = {c: ca.count(c) for c in ca}
    start = min(range(n), key=lambda i: (freq[ca[i]], i))
    order, seen = [], {start}
    queue = [start]
    while queue:
        a = queue.pop(0)
        order.append(a)
        for b in sorted(adj_a[a], key=lambda j: (freq[ca[j]], j)):
            if b not in seen:
                seen.add(b)
                queue.append(b)
    order += [i for i in range(n) if i not in seen]
    cand = {i: [j for j in range(n) if cb[j] == ca[i]] for i in range(n)}
    mp, used, count = [-1] * n, [False] * n, [0]

    def rec(k):
        if count[0] >= limit:
            return
        if k == n:
            count[0] += 1
            yield list(mp)
            return
        a = order[k]
        mapped_nb = [mp[x] for x in adj_a[a] if mp[x] >= 0]
        for b in cand[a]:
            if used[b]:
                continue
            if any(t not in adj_b[b] for t in mapped_nb):
                continue
            if sum(1 for x in adj_b[b] if used[x]) != len(mapped_nb):
                continue
            mp[a], used[b] = b, True
            yield from rec(k + 1)
            mp[a], used[b] = -1, False
    yield from rec(0)


def map_molecule(src_elements, src_xyz_bohr, tgt_elements, tgt_xyz_bohr, tgt_atoms):
    """Map every atom of the source molecule onto `tgt_atoms` of the target structure.

    Returns {'map': array source index -> target atom (0-based, global), 'rmsd_angstrom', 'n_isomorphisms'}."""
    tgt_atoms = list(tgt_atoms)
    n = len(src_elements)
    if n != len(tgt_atoms):
        raise MappingError(f"the coordinate set's molecule has {n} atoms, the subsystem {len(tgt_atoms)}")
    el_b = [tgt_elements[a] for a in tgt_atoms]
    if sorted(src_elements) != sorted(el_b):
        raise MappingError("the two molecules have different element compositions")
    xb = np.asarray(tgt_xyz_bohr)[tgt_atoms]
    bonds_a = co.infer_bonds(list(src_elements), np.asarray(src_xyz_bohr))
    bonds_b = co.infer_bonds(el_b, xb)
    if len(bonds_a) != len(bonds_b):
        raise MappingError(f"different bond counts ({len(bonds_a)} vs {len(bonds_b)}): not the same bonded molecule")
    adj_a, adj_b = _adjacency(n, bonds_a), _adjacency(n, bonds_b)
    xa = np.asarray(src_xyz_bohr) * BOHR_TO_ANG
    xb = xb * BOHR_TO_ANG
    best, count = None, 0
    for iso in isomorphisms(list(src_elements), adj_a, el_b, adj_b):
        count += 1
        _, r = kabsch_rmsd(xa, xb[iso])
        if best is None or r < best[0]:
            best = (r, iso)
    if best is None:
        raise MappingError("the two molecules are not isomorphic (different bonding)")
    return {"map": np.array([tgt_atoms[j] for j in best[1]]), "rmsd_angstrom": best[0], "n_isomorphisms": count,
            "truncated": count >= MAX_ISOMORPHISMS}


def transfer(coords, mapping):
    """Renumber a coordinate set with mapping[source atom] = target atom."""
    return [co.Coordinate(c.label, c.kind, [(co.Primitive(p.kind, tuple(int(mapping[a]) for a in p.atoms)), w)
                                            for p, w in c.components], c.rule, c.id) for c in coords]


def value_differences(src_coords, src_xyz, tgt_coords, tgt_xyz):
    """Largest change of any primitive value between the two structures (angstrom, degrees)."""
    out = {"stretch_angstrom": 0.0, "bend_degrees": 0.0, "torsion_or_oop_degrees": 0.0}
    for cs, ct in zip(src_coords, tgt_coords):
        for (ps, _), (pt, _) in zip(cs.components, ct.components):
            a, b = co.value(ps, src_xyz), co.value(pt, tgt_xyz)
            if not np.isfinite(b):
                raise MappingError(f"coordinate {ct.label} is undefined at the target geometry")
            if ps.kind == "stretch":
                out["stretch_angstrom"] = max(out["stretch_angstrom"], abs(a - b) * BOHR_TO_ANG)
            else:
                d = abs(np.degrees((a - b + np.pi) % (2 * np.pi) - np.pi))
                key = "bend_degrees" if ps.kind == "bend" else "torsion_or_oop_degrees"
                out[key] = max(out[key], float(d))
    return out
