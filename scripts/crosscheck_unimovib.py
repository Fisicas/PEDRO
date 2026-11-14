"""Cross-check PEDRO against UniMoVib 1.5.4 on our own 614-atom ONIOM model (benzaldehyde in 200 UFF waters).

    python scripts/crosscheck_unimovib.py make --log model.log --workdir crosscheck_runs
    python scripts/crosscheck_unimovib.py run --log model.log --workdir crosscheck_runs --executable /path/to/unimovib
    python scripts/crosscheck_unimovib.py compare --log model.log --workdir crosscheck_runs

All output paths must be new. The external program is installed separately. No Gaussian calculation is launched.

A. Harmonic frequencies and IR intensities of the whole model: UniMoVib vs PEDRO vs Gaussian's printed values.
B. rev-GSVA (Tao et al. 2021) subsystem Hessian of benzaldehyde (atoms 1-14): UniMoVib's gsva.almode vs PEDRO's
   implementation of the published formula with the same pseudo-inverse cut-off (1e-12, UniMoVib's DiagIv).
UniMoVib writes data files from exactly what PEDRO reads (archive Hessian, geometry, printed masses, dipole derivatives).
"""
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import null_space

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import pedro  # noqa: E402
from pedro import harmonic  # noqa: E402
from pedro.constants import SYMBOLS  # noqa: E402

# Paths are supplied explicitly; no vendor executable or large calculation is bundled.
HERE = LOG = EXE = None
S = []


def block(name, values, per=5):
    v = np.asarray(values, float).reshape(-1)
    return name + "\n" + "\n".join(" ".join(f"{x:20.10E}" for x in v[i:i + per]) for i in range(0, v.size, per)) + "\n"


def write_umv(path, job):
    n = job.natoms
    txt = (f" own ONIOM model, PEDRO export\nNATM\n{n:6d}\n" + block("AMASS", job.printed_masses) +
           block("ZA", [SYMBOLS.index(e) + 1 for e in job.elements]) + block("XYZ", job.xyz_bohr) +
           block("FFXLT", job.hessian[np.tril_indices(3 * n)]))
    txt += block("APT", job.dipole_derivatives.T.reshape(3, -1).T) if job.dipole_derivatives is not None else "NOAPT\n"
    path.write_text(txt + "END\n", encoding="ascii")


def job_inp(path, gsva):
    t = ' UniMoVib job\n\n $contrl\n   qcprog="unimovib"\n   ifsave=.false.\n' + ("   ifgsva=.true.\n" if gsva else "") + " $end\n\n"
    if gsva:
        t += f' $gsva\n   subsystem="1-{len(S)}"\n $end\n\n'
    path.write_text(t + ' $qcdata\n   fchk="model.umv"\n $end\n\n', encoding="ascii")


def umv_freqs(path, after=None):
    t = Path(path).read_text(errors="replace")
    if after:
        t = t[t.index(after):]
    elif "Generalized Subsystem Vibrational Analysis" in t:
        t = t[:t.index("Generalized Subsystem Vibrational Analysis")]
    f = [float(v) for l in re.findall(r"^\s+Frequencies\s+(.*)$", t, re.M) for v in l.split()]
    ir = [float(v) for l in re.findall(r"^\s+IR intensities\s+(.*)$", t, re.M) for v in l.split()]
    return np.array(f), (np.array(ir) if ir else None)


def read_almode(path):
    tok = Path(path).read_text().replace("D", "E").split()
    n = int(tok[tok.index("NATM") + 1])
    get = lambda key, k: np.array([float(x) for x in tok[tok.index(key) + 1: tok.index(key) + 1 + k]])
    return {"n": n, "mass": get("AMASS", n), "xyz": get("XYZ", 3 * n).reshape(n, 3), "H": get("FFX", 9 * n * n).reshape(3 * n, 3 * n)}


def pedro_gsva(job, cut=1e-12):
    w, U = np.linalg.eigh(job.hessian)
    Fp = (U / np.where(np.abs(w) > cut, w, np.inf)) @ U.T
    r = job.xyz_bohr[S] - job.xyz_bohr[S].mean(0)
    e = np.eye(3)
    Q, _ = np.linalg.qr(np.array([np.tile(e[k], len(S)) for k in range(3)] + [np.cross(e[k], r).reshape(-1) for k in range(3)]).T)
    V = null_space(Q.T)
    idx = (3 * np.array(S)[:, None] + np.arange(3)).reshape(-1)
    Vf = np.zeros((job.hessian.shape[0], V.shape[1]))
    Vf[idx] = V
    return V @ np.linalg.inv(Vf.T @ Fp @ Vf) @ V.T


def main():
    import argparse
    from pedro.cli import _atoms
    global HERE, LOG, EXE, S
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("step", choices=["make", "run", "compare"])
    ap.add_argument("--log", type=Path, required=True, help="Gaussian ONIOM frequency log (read only)")
    ap.add_argument("--workdir", type=Path, required=True, help="new work directory for make; same directory for run/compare")
    ap.add_argument("--subsystem", default="1-14", help="1-based atoms, default benzaldehyde benchmark 1-14")
    ap.add_argument("--executable", type=Path, help="UniMoVib executable, required for run only")
    args = ap.parse_args()
    HERE, LOG = args.workdir.resolve(), args.log.resolve()
    EXE = args.executable.resolve() if args.executable else None
    S = _atoms(args.subsystem)
    if not LOG.is_file():
        ap.error("--log must be an existing Gaussian frequency output")
    if not S or min(S) < 0:
        ap.error("--subsystem must be a nonempty 1-based atom list")
    if args.step == "run" and (EXE is None or not EXE.is_file()):
        ap.error("run requires --executable pointing to an existing UniMoVib executable")
    step = args.step
    if step == "make":
        HERE.mkdir(parents=True, exist_ok=False)
        job = pedro.load(LOG)
        if max(S) >= job.natoms:
            ap.error("subsystem atom lies outside the input structure")
        for sub, g in (("full", False), ("gsva", True)):
            (HERE / sub).mkdir()
            write_umv(HERE / sub / "model.umv", job)
            job_inp(HERE / sub / "job.inp", g)
    elif step == "run":
        for sub in ("full", "gsva"):
            with open(HERE / sub / "job.inp") as fi, open(HERE / sub / "job.out", "x") as fo, open(HERE / sub / "job.err", "x") as fe:
                subprocess.run([str(EXE), "-b"], stdin=fi, stdout=fo, stderr=fe, cwd=HERE / sub, check=True)
            print(sub, "done", flush=True)
    elif step == "compare":
        out = HERE / "crosscheck_results_v01.json"
        if out.exists():
            sys.exit("refusing to overwrite")
        job = pedro.load(LOG)
        r = pedro.analyze(job, pedro.Options(masses="printed", case="F", anharmonic="off", screen=False))
        fu, iru = umv_freqs(HERE / "full" / "job.out")
        fp = r.modes.freqs
        if len(fu) == len(fp) + 6:      # UniMoVib also lists the six rigid-body modes; drop the six smallest |f|
            keep = np.sort(np.argsort(np.abs(fu))[6:])
            fu, iru = fu[keep], (None if iru is None else iru[keep])
        res = {"source": LOG.name, "sha256": job.sha256, "n_atoms": job.natoms,
               "pedro_version": pedro.__version__, "masses": "Gaussian printed, 5 decimals",
               "subsystem_1based": [a + 1 for a in S], "pseudoinverse_cutoff_Eh_bohr2": 1e-12,
               "A": {"n_unimovib": len(fu), "n_pedro": len(fp),
                     "max_abs_freq_unimovib_vs_pedro_cm1": float(np.max(np.abs(np.sort(fu) - np.sort(fp)))),
                     "max_abs_freq_pedro_vs_gaussian_cm1": float(np.max(np.abs(np.sort(fp) - np.sort(job.printed_freqs))))}}
        if iru is not None and r.modes.ir_km_mol is not None and len(iru) == len(fp):
            res["A"]["max_abs_ir_unimovib_vs_pedro_km_mol"] = float(np.max(np.abs(iru[np.argsort(fu)] - r.modes.ir_km_mol[np.argsort(fp)])))
        al = read_almode(HERE / "gsva" / "gsva.almode")
        mine = pedro_gsva(job)
        fg, _ = umv_freqs(HERE / "gsva" / "job.out", after="Properties of Intrinsic Fragmental Vibrations")
        mg = harmonic.analyze(al["H"], al["mass"], al["xyz"])
        res["B"] = {"subsystem_atoms": len(S), "hessian_max_abs_diff_Eh_bohr2": float(np.abs(al["H"] - mine).max()),
                    "hessian_max_abs_element_Eh_bohr2": float(np.abs(al["H"]).max()),
                    "n_freq_unimovib": len(fg), "n_imaginary_unimovib": int((fg < 0).sum()),
                    "max_abs_freq_unimovib_vs_pedro_cm1": float(np.max(np.abs(np.sort(fg) - np.sort(mg.freqs))))}
        out.write_text(json.dumps(res, indent=1), encoding="utf-8")
        print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
