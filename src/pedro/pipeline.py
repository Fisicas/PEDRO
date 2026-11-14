"""One call from an input file to a complete, self-describing result: `pedro.analyze(path)`."""

from __future__ import annotations

import csv
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path

import numpy as np

from . import __version__
from . import basis as bs
from . import coordinates as co
from . import distributions as ds
from . import harmonic
from .io import load as load_job
from .job import Job

EMBED_AUTO_ATOMS = 150       # above this many atoms (or with ONIOM layers) the embedded analysis is the default


@dataclass
class Options:
    coords: str | None = None            # .dd2 or pedro JSON; None = sibling .dd2 if present, else automatic
    coords_molecule: str | None = None   # structure the coordinate set was defined on (default: .fmu/.log next to a .dd2);
                                         # when it differs from the input, the set is transferred by graph isomorphism
    atom_map: str | None = None          # CSV (source_atom,target_atom; 1-based) instead of the automatic transfer
    masses: str = "auto"                 # auto (exact isotope masses when the printed ones agree) | printed (as printed)
    mass_overrides: dict | None = None   # {atom (0-based): mass in u}, e.g. 18O or D substitution
    coord_classes: str | list | None = None  # CSV or list of (name, kind, atoms0): named classes for group totals
    subsystem: list[int] | None = None   # 0-based atoms; None = automatic
    case: str = "auto"                   # auto | isolated | F | Fprime | P
    mobile: list[int] | str | None = None  # for case P: None = subsystem only, "layer:H", or 0-based atoms
    scale: float = 1.0                   # display scale factor for frequencies (never used in the analysis)
    cluster_cm1: float = 1.0
    ted_threshold: float = 10.0          # % shown in the VEDA-style mode descriptions
    groups_csv: str | None = None        # atom,group (1-based atoms) to name environment groups
    window: tuple[float, float] | None = None   # scaled cm^-1 window highlighted in reports
    anharmonic: str = "auto"             # auto (run when a GVPT2 step exists) | off
    symmetry: str = "auto"               # GVPT2 symmetry gauge: auto (size-guarded) | off | force
    screen: bool = True                  # harmonic overtone/Fermi-candidate screen when no usable GVPT2 data
    screen_delta_cm1: float = 10.0

    def to_dict(self):
        return asdict(self)


@dataclass
class Result:
    job: Job
    options: Options
    modes: harmonic.Modes
    basis: bs.Basis
    ted: np.ndarray                      # (modes, coordinates) %
    ted_cluster: np.ndarray
    cluster_id: np.ndarray
    group_totals: dict
    subsystem_ted: np.ndarray            # % of each mode on internal coordinates of S
    subsystem_amplitude: np.ndarray      # % of mass-weighted amplitude on atoms of S
    ped_diag: np.ndarray | None          # Morino-Kuchitsu, NaN for lambda <= 0
    veda_ped: np.ndarray | None
    checks: dict
    descriptions: list[str]
    coord_report: dict
    timings: dict = field(default_factory=dict)
    validation: dict | None = None
    anharmonic: dict | None = None        # GVPT2 layer (status, states, quality, symmetry gauge)
    resonance_screen: dict | None = None  # harmonic overtone / Fermi-candidate screen

    # ------------------------------------------------------------------ convenience
    @property
    def scaled(self):
        return self.modes.freqs * self.options.scale

    def mode_table(self):
        rows = []
        ir = self.modes.ir_km_mol
        for k in range(self.modes.n):
            rows.append({
                "mode": k + 1, "freq_cm1": round(float(self.modes.freqs[k]), 4),
                "scaled_cm1": round(float(self.scaled[k]), 4), "ir_km_mol": None if ir is None else round(float(ir[k]), 4),
                "imaginary": bool(self.modes.freqs[k] < 0), "cluster": int(self.cluster_id[k]),
                "subsystem_ted_pct": round(float(self.subsystem_ted[k]), 2),
                "subsystem_amplitude_pct": round(float(self.subsystem_amplitude[k]), 2),
                "description": self.descriptions[k],
            })
        return rows

    def write(self, out_dir=None):
        from .report import write_all
        return write_all(self, out_dir)


def _read_groups(path):
    out = {}
    with open(path, encoding="utf-8-sig") as f:
        for r in csv.reader(f):
            if r and r[0].strip().isdigit():
                out[int(r[0]) - 1] = r[1].strip()
    return out


def _read_map(path):
    pairs = []
    with open(path, encoding="utf-8-sig") as f:
        for r in csv.reader(f):
            if len(r) >= 2 and r[0].strip().isdigit():
                pairs.append((int(r[0]) - 1, int(r[1]) - 1))
    m = np.full(max(a for a, _ in pairs) + 1, -1)
    for a, b in pairs:
        m[a] = b
    return m


def _apply_masses(job, opts, say):
    import dataclasses
    if opts.masses == "auto" and not opts.mass_overrides:
        return job
    m = np.array(job.masses, float)
    src = job.mass_source
    if opts.masses == "printed":
        if job.printed_masses is None or len(job.printed_masses) != job.natoms:
            raise ValueError("masses='printed' needs the isotope masses printed in a Gaussian log")
        m, src = np.array(job.printed_masses, float), "log (as printed, 5 decimals)"
    elif opts.masses != "auto":
        raise ValueError(f"unknown masses option {opts.masses!r}")
    if opts.mass_overrides:
        for a, v in opts.mass_overrides.items():
            m[int(a)] = float(v)
        src += f"; {len(opts.mass_overrides)} atom(s) overridden: " + ", ".join(
            f"{job.elements[int(a)]}{int(a) + 1}={float(v):.5f}" for a, v in sorted(opts.mass_overrides.items()))
        say(f"Masses: {len(opts.mass_overrides)} override(s)")
    return dataclasses.replace(job, masses=m, mass_source=src)


def _source_structure(coord_path, opts):
    if opts.coords_molecule:
        return Path(opts.coords_molecule)
    p = Path(coord_path)
    if p.suffix.lower() != ".dd2":
        return None
    for ext in (".fmu", ".fchk", ".fch", ".log", ".out"):
        if p.with_suffix(ext).exists():
            return p.with_suffix(ext)
    return None


def _auto_transfer(job, coord_path, opts, S, say):
    """If the coordinate set belongs to another structure, map that molecule onto the subsystem by graph isomorphism."""
    from . import mapping
    src = _source_structure(coord_path, opts)
    if src is None or Path(src).resolve() == Path(job.source).resolve():
        return None, None
    sj = load_job(src)
    if sj.natoms == job.natoms and sj.elements == job.elements:
        return None, None                     # same molecule and numbering: no transfer needed
    m = mapping.map_molecule(sj.elements, sj.xyz_bohr, job.elements, job.xyz_bohr, S.tolist())
    src_set = co.load_set(coord_path, sj.elements)
    coords = mapping.transfer(src_set, m["map"])
    dv = mapping.value_differences(src_set, sj.xyz_bohr, coords, job.xyz_bohr)
    say(f"Coordinate set transferred from {Path(src).name} by graph isomorphism "
        f"({m['n_isomorphisms']} isomorphisms; all-atom RMSD {m['rmsd_angstrom']:.3f} Å)")
    return m["map"], {"coords": coords, "from": str(src), "method": "graph isomorphism, lowest all-atom RMSD",
                      "n_isomorphisms": m["n_isomorphisms"], "rmsd_angstrom": m["rmsd_angstrom"],
                      "source_to_target_1based": [int(a) + 1 for a in m["map"]], "max_value_differences": dv}


def coordinate_classes(internal, spec):
    """Name each internal coordinate by the first matching rule (name, kind, atoms), atoms 0-based:
    kind 'stretch' matches a single stretch between exactly those two atoms; 'any' matches a coordinate touching any
    of the atoms; 'rest' matches everything left. A CSV has rows name,kind,atoms with 1-based atom ranges ('1-5,9')."""
    if isinstance(spec, (str, Path)):
        from .cli import _atoms
        rules = []
        with open(spec, encoding="utf-8-sig") as f:
            for r in csv.reader(f):
                if len(r) >= 2 and r[1].strip() in ("stretch", "any", "rest"):
                    rules.append((r[0].strip(), r[1].strip(), _atoms(r[2]) if len(r) > 2 and r[2].strip() else []))
    else:
        rules = [(n, k, list(a)) for n, k, a in spec]
    out = []
    for c in internal:
        name = "other"
        for n, kind, atoms in rules:
            s = set(atoms)
            if kind == "stretch" and len(c.components) == 1 and c.components[0][0].kind == "stretch" \
                    and set(c.components[0][0].atoms) == s:
                name = n
                break
            if kind == "any" and c.atoms & s:
                name = n
                break
            if kind == "rest":
                name = n
                break
        out.append(name)
    return out


def _choose_case(job, opts, S):
    if opts.case != "auto":
        return opts.case
    if len(S) == job.natoms:
        return "isolated"
    return "F"


def analyze(source, options: Options | None = None, progress=None) -> Result:
    """Run the whole analysis. `source` is a path or a Job. `progress(msg)` receives status lines."""
    opts = options or Options()
    say = progress or (lambda m: None)
    t0 = time.perf_counter()
    tm = {}
    job = source if isinstance(source, Job) else load_job(source)
    job = _apply_masses(job, opts, say)
    tm["read_s"] = time.perf_counter() - t0
    say(f"Read {Path(job.source).name}: {job.natoms} atoms ({job.formula}), {job.kind}")

    # ---- subsystem and coordinates
    frags = None
    coord_path = opts.coords or (str(job.attachments["dd2"]) if "dd2" in job.attachments else None)
    if opts.subsystem is not None:
        S, why = np.array(sorted(opts.subsystem)), "user selection"
    elif coord_path and opts.atom_map:
        amap = _read_map(opts.atom_map)
        S, why = np.array(sorted(a for a in amap if a >= 0)), "atoms named in the atom map"
    elif coord_path and job.natoms <= EMBED_AUTO_ATOMS and not job.layers:
        S, why = np.arange(job.natoms), "whole structure (coordinate set supplied)"
    else:
        S, why, frags = bs.auto_subsystem(job)
    case = _choose_case(job, opts, S)
    say(f"Subsystem: {len(S)} atoms ({why}); basis case {case}")

    t1 = time.perf_counter()
    transfer_info = None
    if coord_path:
        amap = _read_map(opts.atom_map) if opts.atom_map else None
        if amap is None:
            amap, transfer_info = _auto_transfer(job, coord_path, opts, S, say)
        internal = co.load_set(coord_path, job.elements, amap) if amap is None or transfer_info is None else \
            transfer_info.pop("coords")
        bad = [c.label for c in internal if not c.atoms <= set(S.tolist())]
        if bad:
            raise ValueError(f"{len(bad)} coordinates use atoms outside the subsystem (first: {bad[0]})")
        coord_report = {"source": str(coord_path), "n": len(internal), "set_id": co.set_id(internal, job.elements),
                        "transfer": transfer_info}
        need = 3 * len(S) - 6
        if len(internal) != need:
            raise ValueError(f"The coordinate set has {len(internal)} coordinates; the subsystem needs {need}")
    else:
        internal, coord_report = co.build(job.elements, job.xyz_bohr, S.tolist())
    tm["coordinates_s"] = time.perf_counter() - t1
    say(f"Coordinates: {len(internal)} internal ({coord_report['source'] if 'source' in coord_report else 'automatic'})")

    # ---- normal modes
    t2 = time.perf_counter()
    mobile = None
    if case == "P":
        if opts.mobile is None:
            mobile = S.tolist()
        elif isinstance(opts.mobile, str) and opts.mobile.startswith("layer:"):
            lay = opts.mobile.split(":", 1)[1]
            if not job.layers:
                raise ValueError("mobile=layer:… needs ONIOM layers in the input")
            mobile = sorted(set(i for i, l in enumerate(job.layers) if l == lay) | set(S.tolist()))
        else:
            mobile = sorted(set(opts.mobile) | set(S.tolist()))
    method = {"isolated": "full", "F": "full", "Fprime": "unprojected", "P": "partial"}[case]
    modes = harmonic.analyze(job.hessian, job.masses, job.xyz_bohr, job.dipole_derivatives, method=method, mobile=mobile)
    tm["modes_s"] = time.perf_counter() - t2
    say(f"Normal modes: {modes.n} ({modes.n_imaginary} imaginary) in {tm['modes_s']:.1f} s")

    # ---- basis and distributions
    t3 = time.perf_counter()
    ug = _read_groups(opts.groups_csv) if opts.groups_csv else None
    B = bs.build(job, case, internal, S, mobile=mobile, user_groups=ug)
    if opts.coord_classes is not None:
        cls = coordinate_classes(internal, opts.coord_classes)
        B.groups[:len(internal)] = cls
    L = ds.mode_matrix(B.B, modes.Lx)
    Linv, cond = ds.invert(L)
    T = ds.ted(L, Linv)
    tm["ted_s"] = time.perf_counter() - t3
    say(f"TED computed for {T.shape[0]} modes x {T.shape[1]} coordinates in {tm['ted_s']:.1f} s")

    gidx = B.group_index()
    gt = ds.group_totals(T, gidx)
    ni = B.n_internal
    s_ted = T[:, :ni].sum(1)
    q2 = (modes.Lq**2).reshape(job.natoms, 3, -1).sum(1)          # (atoms, modes)
    s_amp = 100.0 * q2[S].sum(0)
    cid = ds.clusters(modes.freqs, opts.cluster_cm1)
    Tc = ds.cluster_sum(T, cid)

    ped_diag = veda = None
    veda_info = None
    F = ds.force_field(Linv, modes.eigenvalues)
    ped_diag = ds.morino_kuchitsu(L, F, modes.eigenvalues)
    if case == "isolated":
        veda, veda_info = ds.veda_compatible_ped(T, np.sign(L).T)

    # ---- per-mode description (VEDA style, internal coordinates of S plus environment total)
    has_rigid = case in ("Fprime", "P")
    env_start = ni + 6 if has_rigid else ni
    env_total = T[:, env_start:].sum(1)
    desc = []
    for k in range(modes.n):
        parts = [f"{lab} {val:.0f}" for lab, val in ds.top_contributions(T[k, :ni], B.labels[:ni], opts.ted_threshold)]
        if has_rigid:
            rigid_total = T[k, ni:ni + 6].sum()
            if abs(rigid_total) >= opts.ted_threshold:
                parts.append(f"[subsystem rigid-body {rigid_total:.0f}]")
        if abs(env_total[k]) >= opts.ted_threshold and case != "isolated":
            parts.append(f"[environment {env_total[k]:.0f}]")
        desc.append(" + ".join(parts) if parts else "(no coordinate ≥ threshold)")

    # ---- checks
    checks = {"reconstruction": harmonic.compare_printed(modes.freqs, job.printed_freqs) if method == "full" else None,
              "n_imaginary": modes.n_imaginary, "condition_L": cond,
              "ted_closure_max_pp": float(np.max(np.abs(T.sum(1) - 100))),
              "ted_coordinate_closure_max_pp": float(np.max(np.abs(T.sum(0) - 100))),
              "subsystem_choice": why, "subsystem_atoms": bs.atom_ranges(S), "case": case,
              "fragments": None if frags is None else [{k: f[k] for k in ("n", "formula", "ranges", "tag")} for f in frags[:12]], "n_modes": modes.n, "n_coordinates": T.shape[1],
              "ped_force_field": "available for modes with positive eigenvalues only",
              "veda_ped": None if veda_info is None else {**veda_info,
                          "status": "diagnostic (imaginary modes present)" if modes.n_imaginary else "ok"}}
    if job.gradient is not None:
        g = np.linalg.norm(job.gradient.reshape(-1, 3), axis=1)
        checks["stationarity"] = {"max_atom_force_Eh_bohr": float(g.max()), "rms_atom_force_Eh_bohr": float(np.sqrt((g**2).mean())),
                                  "max_atom_force_subsystem_Eh_bohr": float(g[S].max()),
                                  "note": "Gaussian's default convergence threshold on the maximum force is 4.5e-4 Eh/bohr"}
    if case in ("Fprime", "P"):
        # exact identity for block-diagonal B: environment TED total = environment amplitude fraction
        cart_total = T[:, ni + 6:].sum(1)
        env_atoms = B.environment if case == "Fprime" else [a for a in (mobile or []) if a not in set(S.tolist())]
        env_amp = 100.0 * q2[env_atoms].sum(0) if len(env_atoms) else np.zeros(modes.n)
        checks["identity_max_dev_pp"] = float(np.max(np.abs(cart_total - env_amp)))
        checks["identity_note"] = "exact identity: environment TED total = environment share of mass-weighted amplitude"
    if case == "F":
        env_amp = 100.0 * q2[B.environment].sum(0)
        checks["identity_F_max_dev_pp"] = float(np.max(np.abs(env_total - env_amp)))
        checks["identity_F_note"] = "approximate in case F (projected modes, rigid motion of S carried by E)"
    # ---- resonances and overtones (isolated from the harmonic results: failures become a status)
    t4 = time.perf_counter()
    anh = scr = None
    if opts.anharmonic != "off" and job.anharmonic_present:
        from .anharmonic import run as run_gvpt2
        say("GVPT2 step found: reading resonances, overtones and combination bands")
        anh = run_gvpt2(job, modes, desc, symmetry=opts.symmetry)
        say(f"GVPT2 layer: {anh['status']}" + (f" ({anh.get('reason')})" if anh.get("reason") else
            f", {anh['n_states']} states, verdicts {anh['verdicts']}"))
    else:
        anh = {"status": "not_applicable", "reason": "switched off" if opts.anharmonic == "off" else
               "no GVPT2 (Freq=Anharmonic) step in the input"}
    if opts.screen and anh.get("status") not in ("ok", "incomplete"):
        from .anharmonic.screen import screen
        scr = screen(modes, s_amp, desc, job.natoms, delta=opts.screen_delta_cm1, window=opts.window, scale=opts.scale)
        say(f"Harmonic resonance screen: {scr['n_candidates']} Fermi-resonance candidates among "
            f"{scr['n_fundamentals_screened']} IR-active subsystem fundamentals")
    checks["anharmonic"] = {k: anh.get(k) for k in ("status", "reason", "pt2_model", "counts", "verdicts") if k in anh}
    tm["anharmonic_s"] = time.perf_counter() - t4
    tm["total_s"] = time.perf_counter() - t0
    res = Result(job=job, options=opts, modes=modes, basis=B, ted=T, ted_cluster=Tc, cluster_id=cid, group_totals=gt,
                 subsystem_ted=s_ted, subsystem_amplitude=s_amp, ped_diag=ped_diag, veda_ped=veda, checks=checks,
                 descriptions=desc, coord_report=coord_report, timings=tm, anharmonic=anh, resonance_screen=scr)
    if "ved" in job.attachments and case == "isolated":
        from .validate import compare_with_veda
        res.validation = compare_with_veda(res, job.attachments["ved"])
        say(f"VEDA comparison: {100 * res.validation['ted_fraction_within_1pp']:.2f}% of TED entries within 1 pp")
    say(f"Done in {tm['total_s']:.1f} s")
    return res
