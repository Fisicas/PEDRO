"""Harmonic normal-mode analysis of a Cartesian Hessian.

full (default)   mass-weight, project out the 6 (5 if linear) rigid-body motions, diagonalize: 3N-6 modes
unprojected      mass-weight and diagonalize without projection: 3N modes (6 near zero)
partial          only the atoms in `mobile` move; all other atoms are infinitely heavy. The H_MM block is
                 mass-weighted and diagonalized without projection: 3N_M modes (partial-Hessian analysis)
See docs/theory.md, sections 1 and 8.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import null_space

from .constants import FREQ_FACTOR, IR_KMMOL_PER_E2_AMU


@dataclass
class Modes:
    freqs: np.ndarray        # cm^-1, ascending; imaginary as negative
    eigenvalues: np.ndarray  # Eh / (bohr^2 u)
    Lq: np.ndarray           # mass-weighted orthonormal eigenvectors (3N, n_modes), zero rows for fixed atoms
    Lx: np.ndarray           # Cartesian displacement per unit normal coordinate, M^-1/2 Lq  (bohr / (bohr u^1/2))
    ir_km_mol: np.ndarray | None
    n_rigid: int             # rigid-body vectors removed (0 for partial/unprojected)
    method: str              # full | unprojected | partial
    mobile: np.ndarray | None = None   # 0-based mobile atoms for partial analyses

    @property
    def n(self) -> int:
        return self.freqs.size

    @property
    def n_imaginary(self) -> int:
        return int(np.sum(self.freqs < 0))


def rigid_basis(xyz, masses):
    """Orthonormal mass-weighted translation/rotation basis (3N x 6, or 3N x 5 for a linear molecule)."""
    com = np.average(xyz, axis=0, weights=masses)
    r = xyz - com
    sm = np.sqrt(masses)
    vecs = [(sm[:, None] * e[None, :]).reshape(-1) for e in np.eye(3)]
    vecs += [(sm[:, None] * np.cross(e, r)).reshape(-1) for e in np.eye(3)]
    u, s, _ = np.linalg.svd(np.column_stack(vecs), full_matrices=False)
    return u[:, s > s[0] * 1e-8]


def _finish(lam, Lq, msq, dipole_derivatives, n_rigid, method, mobile=None):
    order = np.argsort(lam)
    lam, Lq = lam[order], Lq[:, order]
    # deterministic phase: the largest-magnitude component of each mode is positive
    idx = np.argmax(np.abs(Lq), axis=0)
    Lq = Lq * np.sign(Lq[idx, np.arange(Lq.shape[1])])[None, :]
    Lx = Lq / msq[:, None]
    freqs = np.sign(lam) * np.sqrt(np.abs(lam)) * FREQ_FACTOR
    ir = None
    if dipole_derivatives is not None:
        dmu = Lx.T @ dipole_derivatives
        ir = IR_KMMOL_PER_E2_AMU * np.sum(dmu**2, axis=1)
    return Modes(freqs=freqs, eigenvalues=lam, Lq=Lq, Lx=Lx, ir_km_mol=ir, n_rigid=n_rigid, method=method,
                 mobile=mobile)


def analyze(hessian, masses, xyz, dipole_derivatives=None, *, method="full", mobile=None) -> Modes:
    masses = np.asarray(masses, float)
    msq = np.repeat(np.sqrt(masses), 3)
    if method == "partial":
        if mobile is None or len(mobile) == 0:
            raise ValueError("partial analysis needs a non-empty set of mobile atoms")
        mobile = np.array(sorted(set(int(a) for a in mobile)))
        c = (3 * mobile[:, None] + np.arange(3)[None, :]).reshape(-1)
        Hs = hessian[np.ix_(c, c)] / np.outer(msq[c], msq[c])
        lam, W = np.linalg.eigh((Hs + Hs.T) / 2)
        Lq = np.zeros((hessian.shape[0], W.shape[1]))
        Lq[c] = W
        return _finish(lam, Lq, msq, dipole_derivatives, 0, "partial", mobile)
    Hmw = hessian / np.outer(msq, msq)
    Hmw = (Hmw + Hmw.T) / 2
    if method == "unprojected":
        lam, Lq = np.linalg.eigh(Hmw)
        return _finish(lam, Lq, msq, dipole_derivatives, 0, "unprojected")
    if method != "full":
        raise ValueError(f"unknown method {method!r}")
    Q = rigid_basis(xyz, masses)
    V = null_space(Q.T)
    lam, W = np.linalg.eigh(V.T @ Hmw @ V)
    return _finish(lam, V @ W, msq, dipole_derivatives, Q.shape[1], "full")


def compare_printed(freqs, printed, window=None):
    """Max |computed - printed| over modes in ascending order (all modes and optionally a window, cm^-1)."""
    if printed is None or len(printed) == 0:
        return None
    a, b = np.sort(freqs), np.sort(np.asarray(printed, float))
    n = min(a.size, b.size)
    out = {"n_computed": int(a.size), "n_printed": int(b.size), "max_abs_cm1": float(np.max(np.abs(a[:n] - b[:n])))}
    if window is not None:
        sel = (b[:n] >= window[0]) & (b[:n] <= window[1])
        out["window_cm1"] = list(window)
        out["max_abs_window_cm1"] = float(np.max(np.abs(a[:n][sel] - b[:n][sel]))) if sel.any() else None
    return out
