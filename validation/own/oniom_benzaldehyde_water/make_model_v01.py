"""Build a ~600-atom ONIOM test model: benzaldehyde (high layer, B3LYP/6-31G(d)) in a sphere of water (low layer, UFF).

    python make_model_v01.py N_WATERS

Waters are placed with random orientations at random points inside a sphere around the solute, rejecting any
heavy-atom contact closer than 2.6 A (2.2 A to H). Deterministic (seed 20261113). Writes
benzaldehyde_water_oniom.gjf; refuses to overwrite. The structure is a test object for PEDRO's ONIOM handling,
timing and identities, not a physical solvation model.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SOLUTE = [  # benzaldehyde, angstrom (same starting structure as the GVPT2 example)
    ("C", 0.000, 1.395, 0.0), ("C", 1.208, 0.6975, 0.0), ("C", 1.208, -0.6975, 0.0), ("C", 0.000, -1.395, 0.0),
    ("C", -1.208, -0.6975, 0.0), ("C", -1.208, 0.6975, 0.0), ("H", 2.148, 1.240, 0.0), ("H", 2.148, -1.240, 0.0),
    ("H", 0.000, -2.480, 0.0), ("H", -2.148, -1.240, 0.0), ("H", -2.148, 1.240, 0.0), ("C", 0.000, 2.875, 0.0),
    ("O", 1.052, 3.4825, 0.0), ("H", -0.961, 3.430, 0.0)]
WATER = np.array([[0.0, 0.0, 0.0], [0.7572, 0.5865, 0.0], [-0.7572, 0.5865, 0.0]])


def rot(rng):
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    a, b, c, d = q
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)], [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])


def main():
    nw = int(sys.argv[1])
    out = HERE / "benzaldehyde_water_oniom.gjf"
    if out.exists():
        sys.exit("refusing to overwrite " + out.name)
    rng = np.random.default_rng(20261113)
    el = [a[0] for a in SOLUTE]
    xyz = np.array([a[1:] for a in SOLUTE])
    xyz -= xyz.mean(axis=0)
    radius = (3 * (nw / 0.0334 + 150) / (4 * np.pi)) ** (1 / 3)
    placed, tries = 0, 0
    while placed < nw:
        tries += 1
        if tries > 2_000_000:
            sys.exit(f"only {placed} waters placed")
        p = rng.uniform(-radius, radius, 3)
        if np.linalg.norm(p) > radius:
            continue
        w = WATER @ rot(rng).T + p
        d = np.linalg.norm(xyz[:, None, :] - w[None, :, :], axis=2)
        heavy = np.array([e != "H" for e in el])
        if d[heavy][:, 0].min(initial=9) < 2.6 or d.min() < 2.0:
            continue
        el += ["O", "H", "H"]
        xyz = np.vstack([xyz, w])
        placed += 1
    ns = len(SOLUTE)
    lines = ["%nprocshared=8", "%mem=6GB", "%chk=benzaldehyde_water_oniom.chk",
             "#p ONIOM(B3LYP/6-31G(d):UFF) Opt Freq", "",
             f"benzaldehyde in {nw} UFF waters: ONIOM test model for PEDRO ({len(el)} atoms)", "", "0 1 0 1 0 1"]
    for i, (e, x) in enumerate(zip(el, xyz)):
        lines.append(f"{e:<3s} 0 {x[0]:12.6f} {x[1]:12.6f} {x[2]:12.6f} {'H' if i < ns else 'L'}")
    lines += ["", ""]
    out.write_text("\n".join(lines), encoding="ascii")
    print(f"{len(el)} atoms, sphere radius {radius:.1f} A, {tries} trials")


if __name__ == "__main__":
    main()
