"""Energy distributions on a square, non-redundant coordinate basis. Derivations: docs/theory.md, sections 3-6.

With Cartesian normal modes Lx (columns) and the B matrix of the basis,
    L = B Lx                         coordinate amplitude per unit normal coordinate (basis x modes)
    TED[k, i] = 100 L_ik (L^-1)_ki   total energy distribution (row partition; Keresztury-Jalsovszky)
The TED is defined for every mode, including imaginary ones, and uses no force constants. Each mode row and
each coordinate column sums to exactly 100. Entries are signed and are never clipped or renormalized.

Force-constant conventions use F = L^-T diag(lambda) L^-1 and are reported only for modes with lambda > 0:
    morino_kuchitsu[k, i] = 100 F_ii L_ik^2 / lambda_k   (diagonal PED; rows need not sum to 100)
    pair terms            = 200 F_ij L_ik L_jk / lambda_k (i < j); diagonal + pairs = 100
"""

from __future__ import annotations

import numpy as np


class DistributionError(Exception):
    pass


def mode_matrix(B, Lx):
    L = B @ Lx
    if L.shape[0] != L.shape[1]:
        raise DistributionError(f"The coordinate basis is not square: {L.shape[0]} coordinates for {L.shape[1]} modes")
    return L


def invert(L):
    s = np.linalg.svd(L, compute_uv=False)
    cond = float(s[0] / s[-1]) if s[-1] > 0 else float("inf")
    if not np.isfinite(cond) or cond > 1e12:
        raise DistributionError(f"The mode matrix L = B Lx is singular (condition number {cond:.2e}); the "
                                "coordinates are not independent for these modes")
    return np.linalg.solve(L, np.eye(L.shape[0])), cond


def ted(L, Linv):
    """TED[mode, coordinate] in percent."""
    return 100.0 * L.T * Linv


def force_field(Linv, eigenvalues):
    F = Linv.T @ (np.asarray(eigenvalues)[:, None] * Linv)
    return (F + F.T) / 2


def morino_kuchitsu(L, F, eigenvalues):
    """Diagonal PED[mode, coordinate] in percent; NaN rows for modes with lambda <= 0."""
    lam = np.asarray(eigenvalues, float)
    out = np.full((lam.size, L.shape[0]), np.nan)
    ok = lam > 0
    out[ok] = (100.0 * np.diag(F)[:, None] * L**2 / np.where(ok, lam, 1.0)[None, :]).T[ok]
    return out


def pair_terms(L, F, eigenvalues, k):
    """Unordered pair contributions (upper triangle) for mode k; with the diagonal terms they close to 100."""
    lam = float(eigenvalues[k])
    if lam <= 0:
        raise DistributionError("pair terms are defined only for modes with a positive eigenvalue")
    v = L[:, k]
    return np.triu(200.0 * F * np.outer(v, v) / lam, 1)


VEDA_TOL = 1e-8
VEDA_MAX_IT = 200000


def veda_compatible_ped(T, phase, tol=VEDA_TOL, max_iterations=VEDA_MAX_IT):
    """VEDA 4.99e-compatible PED: Sinkhorn balancing of |TED| so every mode row and every coordinate column of
    the magnitudes sums to 100, with signs taken from the coordinate phase sign(L). This rule was inferred from
    VEDA's printed tables for several organic test molecules; VEDA's source code is not available."""
    A = np.abs(np.asarray(T, float))
    if A.shape[0] != A.shape[1] or np.any(A.sum(0) == 0) or np.any(A.sum(1) == 0):
        raise DistributionError("VEDA-compatible balancing needs a square TED without all-zero rows or columns")
    X = A.copy()
    err = np.inf
    for it in range(1, max_iterations + 1):
        X /= X.sum(axis=1, keepdims=True)
        X /= X.sum(axis=0, keepdims=True)
        err = max(np.max(np.abs(X.sum(1) - 1)), np.max(np.abs(X.sum(0) - 1)))
        if err <= tol:
            break
    return 100.0 * np.sign(phase) * X, {"converged": bool(err <= tol), "iterations": it, "max_marginal_error_pp": 100 * err}


# ----------------------------------------------------------------------------- near-degenerate clusters
def clusters(freqs, threshold_cm1=1.0):
    """Cluster id per mode: consecutive modes (ascending) closer than the threshold share an id."""
    f = np.asarray(freqs)
    order = np.argsort(f)
    cid = np.zeros(f.size, int)
    c = 0
    for j in range(1, f.size):
        if f[order[j]] - f[order[j - 1]] >= threshold_cm1:
            c += 1
        cid[order[j]] = c
    return cid


def cluster_sum(T, cid):
    """TED summed over each cluster, broadcast back to the member rows (invariant to rotations in the cluster)."""
    out = np.empty_like(T)
    for c in np.unique(cid):
        m = cid == c
        out[m] = T[m].sum(0)
    return out


def group_totals(T, groups: dict[str, np.ndarray]):
    """Sum TED columns per named group. A group total is invariant to any invertible re-parameterization of the
    coordinates inside the group (docs/theory.md, section 5.4)."""
    return {g: T[:, idx].sum(1) for g, idx in groups.items() if len(idx)}


def top_contributions(row, labels, threshold=10.0, limit=6):
    idx = np.argsort(-np.abs(row))
    out = [(labels[i], float(row[i])) for i in idx[:limit] if abs(row[i]) >= threshold]
    return out
