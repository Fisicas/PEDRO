"""Band-level quantities derived from a TED, as used in papers (docs/theory.md, section 11).

    stretch_index(res, a, b)          index of the single stretch coordinate between atoms a and b
    band_metrics(res, j, window)      strongest mode (largest TED x IR), TED- and TED x IR-weighted centroids, TED sum
    class_totals(res, k)              TED of mode k summed over named coordinate classes and environment groups
    follow(parent, child, k, atoms)   where parent mode k goes in another calculation (isotopologue, partial Hessian)
    match(a, b, atoms)                one-to-one assignment of modes between two calculations (Hungarian)
All overlaps use mass-weighted mode vectors restricted to `atoms` and normalized there.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import linear_sum_assignment


def stretch_index(res, a, b):
    for j, c in enumerate(res.basis.internal):
        if len(c.components) == 1 and c.components[0][0].kind == "stretch" and set(c.components[0][0].atoms) == {a, b}:
            return j
    raise KeyError(f"no single stretch coordinate between atoms {a + 1} and {b + 1}")


def band_metrics(res, j, window, scale=None, strict=True):
    """Within a scaled-wavenumber window: the mode with the largest TED_j x IR ('strongest'), the TED-weighted and
    TED x IR-weighted mean positions, and the summed TED. TED enters as a signed fraction, as computed."""
    s = res.options.scale if scale is None else scale
    f = res.modes.freqs * s
    has_ir = res.modes.ir_km_mol is not None
    ir = res.modes.ir_km_mol if has_ir else np.ones(res.modes.n)   # no intensities: rank by TED alone
    t = res.ted[:, j] / 100.0
    lo, hi = window
    m = (f > lo) & (f < hi) if strict else (f >= lo) & (f <= hi)
    idx = np.flatnonzero(m)
    if not idx.size:
        return None
    w = t[idx] * ir[idx]
    i = int(idx[np.argmax(w)])
    return {"window": [lo, hi], "scale": s, "weights": "TED x IR" if has_ir else "TED only (no IR intensities)",
            "strongest": {"mode": i + 1, "index": i, "scaled_cm1": float(f[i]),
                          "ir_km_mol": float(ir[i]) if has_ir else None, "ted_pct": float(100 * t[i])},
            "centroid_ted": float((t[idx] * f[idx]).sum() / t[idx].sum()),
            "centroid_ted_ir": float((w * f[idx]).sum() / w.sum()) if has_ir else None,
            "ted_sum_pct": float(100 * t[idx].sum()),
            "modes_ted_ge5": [(int(k + 1), float(f[k]), float(100 * t[k]), float(ir[k]) if has_ir else None)
                              for k in idx if 100 * t[k] >= 5]}


def class_totals(res, k):
    out = {}
    for g, idx in res.basis.group_index().items():
        out[g] = float(res.ted[k, idx].sum())
    return out


def _vectors(res, atoms):
    rows = (3 * np.asarray(atoms)[:, None] + np.arange(3)[None, :]).reshape(-1)
    v = res.modes.Lq[rows]
    n = np.linalg.norm(v, axis=0)
    n[n == 0] = 1
    return v / n


def overlaps(a, b, atoms):
    """|<a_k|b_l>| for every pair of modes, mass-weighted on `atoms` (each calculation with its own masses)."""
    return np.abs(_vectors(a, atoms).T @ _vectors(b, atoms))


def follow(parent, child, k, atoms, n_best=2):
    o = overlaps(parent, child, atoms)[k]
    best = np.argsort(-o)[:n_best]
    return [{"mode": int(b + 1), "index": int(b), "freq_cm1": float(child.modes.freqs[b]), "overlap": float(o[b])}
            for b in best]


def match(a, b, atoms, threshold=0.8):
    """One-to-one assignment maximizing total overlap; pairs below `threshold` are flagged, not dropped."""
    O = overlaps(a, b, atoms)
    r, c = linear_sum_assignment(-O)
    return [{"a": int(i + 1), "b": int(j + 1), "overlap": float(O[i, j]), "flag": bool(O[i, j] < threshold)}
            for i, j in zip(r, c)]
