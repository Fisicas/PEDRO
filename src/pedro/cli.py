"""Command line: `pedro analyze FILE`, `pedro gui`, `pedro veda-compare FOLDER`."""

from __future__ import annotations

import argparse
import json
import sys


def _atoms(spec: str | None):
    """'1-31,40,42-44' (1-based, as in Gaussian) -> sorted 0-based list."""
    if not spec:
        return None
    out = set()
    for part in spec.replace(" ", "").split(","):
        if "-" in part:
            a, b = part.split("-")
            out.update(range(int(a) - 1, int(b)))
        elif part:
            out.add(int(part) - 1)
    return sorted(out)


def main(argv=None):
    from . import __version__
    p = argparse.ArgumentParser(prog="pedro", description="Potential/total energy distributions from Cartesian Hessians")
    p.add_argument("--version", action="version", version=f"PEDRO {__version__}")
    sub = p.add_subparsers(dest="cmd")
    a = sub.add_parser("analyze", help="analyze a Gaussian .log/.fchk, VEDA .fmu, or .npz/.json Hessian")
    a.add_argument("file")
    a.add_argument("--coords", help=".dd2 or PEDRO JSON coordinate set (default: sibling .dd2, else automatic)")
    a.add_argument("--coords-molecule", help="structure the coordinate set belongs to (default: .fmu/.log next to the "
                   ".dd2); the set is transferred onto the subsystem by graph isomorphism")
    a.add_argument("--atom-map", help="CSV source_atom,target_atom (1-based) instead of the automatic transfer")
    a.add_argument("--masses", default="auto", choices=["auto", "printed"],
                   help="printed: Gaussian's 5-decimal masses exactly as printed")
    a.add_argument("--mass-override", help="CSV atom,mass (1-based atoms, u), e.g. 18O or D substitution")
    a.add_argument("--coord-classes", help="CSV name,kind,atoms (kind stretch|any|rest) naming coordinate classes")
    a.add_argument("--subsystem", help="atoms of the molecule to describe, 1-based, e.g. '1-31,40'")
    a.add_argument("--case", default="auto", choices=["auto", "isolated", "F", "Fprime", "P"])
    a.add_argument("--mobile", help="partial Hessian (case P): 'layer:H' or 1-based atoms; default subsystem only")
    a.add_argument("--scale", type=float, default=1.0, help="display scale factor for frequencies")
    a.add_argument("--groups", help="CSV atom,group (1-based) naming environment groups")
    a.add_argument("--out", help="output folder (must not exist or be empty)")
    a.add_argument("--anharmonic", default="auto", choices=["auto", "off"], help="GVPT2 resonance/overtone layer")
    a.add_argument("--symmetry", default="auto", choices=["auto", "off", "force"], help="GVPT2 symmetry gauge")
    a.add_argument("--screen-delta", type=float, default=10.0, help="harmonic Fermi-candidate window (cm-1)")
    a.add_argument("--window", type=float, nargs=2, metavar=("LOW", "HIGH"), help="scaled wavenumber window (cm-1)")
    a.add_argument("--quiet", action="store_true")
    g = sub.add_parser("gui", help="open the drag-and-drop window")
    g.add_argument("files", nargs="*", help="optional file to analyze on start")
    v = sub.add_parser("veda-compare", help="run PEDRO next to VEDA on folders holding .fmu + .dd2 + .ved")
    v.add_argument("folders", nargs="+")
    args = p.parse_args(argv)
    if args.cmd in (None, "gui"):
        from .gui.app import main as gui_main
        return gui_main()
    if args.cmd == "analyze":
        from .pipeline import Options, analyze
        mob = args.mobile if (args.mobile or "").startswith("layer:") else _atoms(args.mobile)
        mo = None
        if args.mass_override:
            import csv as _csv
            with open(args.mass_override, encoding="utf-8-sig") as fh:
                mo = {int(r[0]) - 1: float(r[1]) for r in _csv.reader(fh) if len(r) >= 2 and r[0].strip().isdigit()}
        opts = Options(coords=args.coords, atom_map=args.atom_map, coords_molecule=args.coords_molecule,
                       masses=args.masses, mass_overrides=mo, coord_classes=args.coord_classes, subsystem=_atoms(args.subsystem), case=args.case,
                       mobile=mob, scale=args.scale, groups_csv=args.groups, anharmonic=args.anharmonic,
                       symmetry=args.symmetry, screen_delta_cm1=args.screen_delta,
                       window=tuple(args.window) if args.window else None)
        res = analyze(args.file, opts, progress=None if args.quiet else (lambda m: print(m, flush=True)))
        out = res.write(args.out)
        print(f"Results: {out}")
        return 0
    if args.cmd == "veda-compare":
        from pathlib import Path
        from .pipeline import analyze
        rows = []
        for folder in args.folders:
            fmu = next(Path(folder).glob("*.fmu"))
            res = analyze(fmu)
            v = res.validation or {}
            rows.append({"case": Path(folder).name, **{k: v.get(k) for k in (
                "ted_fraction_within_1pp", "ted_max_abs_diff_pp", "ted_cluster_summed_max_abs_diff_pp",
                "ped_fraction_within_1pp", "freq_max_abs_diff_cm1")}, "n_imaginary": res.checks["n_imaginary"]})
        print(json.dumps(rows, indent=1))
        return 0
    p.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
