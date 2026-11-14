"""Reproducible validation run on public data. Writes validation/validation_<date>_vNN.json and docs/validation.md.

    python scripts/run_validation.py [--veda-examples DIR] [--oniom LOG ...]

Everything it reads by default is in this repository and was computed for PEDRO (no third-party or unpublished data):
  validation/veda_reference/<case>/   VEDA 4.99e run on our own Hessians: <case>.fmu, .dd2, .ved
  validation/published/<name>/        our recalculation of a published example + the published values (published.json)
  examples/, tests/data/              our Gaussian and ORCA outputs (frequency reproduction; GVPT2 logs)
Optional:
  --veda-examples DIR   VEDA 4's own example folder (examp_1/, examp_2/ of the VEDA download). These files belong to
                        VEDA's author and are not redistributed; only the agreement statistics are recorded.
  --oniom LOG ...       large ONIOM frequency logs (timing, memory, identities in cases F, Fprime, P)
  --extra-gvpt2 LOG ... more Gaussian Freq=Anharmonic logs
Input files are only read; their SHA-256 hashes are recorded.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import platform
import re
import shutil
import sys
import tempfile
import time
import tracemalloc
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import pedro  # noqa: E402
from pedro.io.veda import read_ved  # noqa: E402

VEDA_EXAMPLES = [("examp_1/examp_1", "a1-91.fmv", "a1-91"), ("examp_1/examp_1", "a1-91.fmv", "a1-91a"),
                 ("examp_2/examp_2", "neopent.fmu", "neopent"), ("examp_2/examp_2", "neopent.fmu", "neopent1"),
                 ("examp_2/examp_2", "neopent.fmu", "neopent2"), ("examp_2/examp_2", "neopent.fmu", "neopent3")]
VKEYS = ("veda_version", "ted_fraction_within_1pp", "ted_fraction_within_0p5pp", "ted_max_abs_diff_pp",
         "ted_cluster_summed_max_abs_diff_pp", "ted_rms_diff_pp", "ped_fraction_within_1pp",
         "ped_fraction_magnitude_rounds_to_veda", "ped_max_abs_diff_pp", "freq_max_abs_diff_cm1")


def rel(p):
    try:
        return Path(p).resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return Path(p).name


def veda_row(fmu, label, source):
    r = pedro.analyze(fmu, pedro.Options(anharmonic="off", screen=False))
    v = r.validation
    return {"case": label, "source": source, "formula": r.job.formula, "n_atoms": r.job.natoms, "n_modes": r.modes.n,
            "n_imaginary": r.modes.n_imaginary, "sha256_input": r.job.sha256, **{k: v[k] for k in VKEYS},
            "largest_differences": v["largest_differences"][:3]}


def published_rows(folder):
    spec = json.loads((folder / "published.json").read_text(encoding="utf-8"))
    out = {"name": folder.name, "title": spec["title"], "reference": spec["reference"],
           "published_level": spec["published_level"], "our_calculation": spec["our_calculation"], "tables": []}
    for tname, t in spec["tables"].items():
        r = pedro.analyze(folder / spec["input"], pedro.Options(coords=str(folder / t["dd2"]), anharmonic="off", screen=False))
        vdir = ROOT / "validation" / "veda_reference" / t.get("veda_case", "")
        vp = None
        if t.get("veda_case") and (vdir / f"{t['veda_case']}.ved").exists():
            vp = np.abs(read_ved(vdir / f"{t['veda_case']}.ved")["ped"][::-1])
        f = r.modes.freqs
        rows = []
        for s, entries in t["ped"].items():
            j = int(s) - 1
            for fr, p in entries:
                k = int(np.argmin(np.abs(f - fr)))
                rows.append({"coord": int(s), "label": r.basis.internal[j].label, "published_cm1": fr,
                             "ours_cm1": round(float(f[k]), 1), "published": p, "pedro": round(abs(float(r.veda_ped[k, j])), 1),
                             "veda_on_ours": None if vp is None else float(vp[k, j])})
        d = [abs(x["pedro"] - x["published"]) for x in rows]
        out["tables"].append({"table": tname, "rows": rows, "n": len(rows), "within_1pp": int(sum(x <= 1 for x in d)),
                              "max_pp": max(d), "freq_max_abs_diff_cm1": max(abs(x["ours_cm1"] - x["published_cm1"]) for x in rows),
                              "pedro_vs_veda_on_ours_max_pp": None if vp is None else
                              max(abs(x["pedro"] - x["veda_on_ours"]) for x in rows)})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--veda-examples")
    ap.add_argument("--oniom", nargs="*", default=[])
    ap.add_argument("--extra-gvpt2", nargs="*", default=[])
    ap.add_argument("--output-dir", type=Path, help="write versioned JSON/Markdown here instead of refreshing docs/validation.md")
    a = ap.parse_args()
    rec = {"pedro": pedro.__version__, "date": dt.datetime.now().isoformat(timespec="seconds"),
           "python": platform.python_version(), "numpy": np.__version__, "platform": platform.platform(),
           "veda_own": [], "veda_examples": [], "published": [], "logs": [], "oniom": [], "gvpt2": []}
    for case in sorted(p for p in (ROOT / "validation" / "veda_reference").iterdir() if p.is_dir()):
        rec["veda_own"].append(veda_row(next(case.glob("*.fmu")), case.name, rel(case)))
        print("veda_own", case.name, flush=True)
    if a.veda_examples:
        base = Path(a.veda_examples)
        with tempfile.TemporaryDirectory() as tmp:
            for folder, inp, stem in VEDA_EXAMPLES:
                d = Path(tmp) / stem
                d.mkdir()
                shutil.copyfile(base / folder / inp, d / (stem + Path(inp).suffix))
                for ext in (".dd2", ".ved"):
                    shutil.copyfile(base / folder / (stem + ext), d / (stem + ext))
                rec["veda_examples"].append(veda_row(d / (stem + Path(inp).suffix), stem, "VEDA 4 distribution, " + folder.split("/")[0]))
                print("veda_examples", stem, flush=True)
    for folder in sorted(p for p in (ROOT / "validation" / "published").iterdir() if (p / "published.json").exists()):
        rec["published"].append(published_rows(folder))
        print("published", folder.name, flush=True)
    logs = sorted((ROOT / "examples").glob("*.log")) + sorted((ROOT / "examples").glob("*.hess")) + \
        sorted((ROOT / "tests" / "data").glob("*.hess")) + sorted((ROOT / "validation" / "published").glob("*/*.log"))
    for log in logs:
        r = pedro.analyze(log, pedro.Options(anharmonic="off", screen=False))
        rec["logs"].append({"file": rel(log), "program": r.job.program or r.job.kind, "sha256": r.job.sha256,
                            "n_atoms": r.job.natoms, "case": r.basis.case,
                            "freq_max_abs_diff_cm1": (r.checks["reconstruction"] or {}).get("max_abs_cm1"),
                            "n_imaginary": r.modes.n_imaginary, "masses": r.job.mass_source})
        print("log", log.name, flush=True)
    for log in a.oniom:
        job = pedro.load(log)
        row = {"file": rel(log), "sha256": job.sha256, "n_atoms": job.natoms, "runs": {}}
        for case, mobile in (("F", None), ("Fprime", None), ("P", "layer:H")):
            tracemalloc.start()
            t = time.perf_counter()
            r = pedro.analyze(job, pedro.Options(case=case, mobile=mobile, anharmonic="off", screen=False))
            dt_s = time.perf_counter() - t
            peak = tracemalloc.get_traced_memory()[1] / 1e6
            tracemalloc.stop()
            c = r.checks
            row["runs"][case] = {"seconds": round(dt_s, 2), "peak_MB": round(peak), "n_modes": r.modes.n,
                                 "n_imaginary": r.modes.n_imaginary, "n_subsystem_atoms": len(r.basis.subsystem),
                                 "condition_L": c["condition_L"], "ted_closure_max_pp": c["ted_closure_max_pp"],
                                 "freq_max_abs_diff_cm1": (c["reconstruction"] or {}).get("max_abs_cm1"),
                                 "identity_max_dev_pp": c.get("identity_max_dev_pp", c.get("identity_F_max_dev_pp"))}
            print("oniom", Path(log).name, case, f"{dt_s:.1f}s", flush=True)
        row["max_atom_force_Eh_bohr"] = c.get("stationarity", {}).get("max_atom_force_Eh_bohr")
        rec["oniom"].append(row)
    gv = [p for p in sorted((ROOT / "examples").glob("*.log")) if "gvpt2" in p.name.lower()] + [Path(p) for p in a.extra_gvpt2]
    for log in gv:
        t = time.perf_counter()
        r = pedro.analyze(log)
        sec = time.perf_counter() - t
        an = r.anharmonic
        row = {"file": rel(log), "sha256": r.job.sha256, "n_atoms": r.job.natoms, "seconds": round(sec, 2),
               "status": an["status"], "reason": an.get("reason")}
        if an.get("states"):
            from pedro.anharmonic import gvpt2_step_text, parse
            from pedro.anharmonic.screen import screen
            mp = {int(k): v for k, v in an["mapping"]["gaussian_to_pedro_1based"].items()}
            true = {(mp[i], tuple(sorted((mp[j], mp[k])))) for (i, j, k), *_ in parse(gvpt2_step_text(log))["fermi"]}
            sc = screen(r.modes, np.full(r.modes.n, 100.0), r.descriptions, r.job.natoms)
            keys = set()
            for c in sc["candidates"]:
                n = [int(v) for v in re.findall(r"\d+", c["partner"])]
                keys.add((c["fundamental"], tuple(sorted(n * 2 if len(n) == 1 else n))))
            row.update({"pt2_model": an["pt2_model"], "counts": an["counts"], "verdicts": an["verdicts"],
                        "mapping_residual_cm1": an["mapping"]["max_residual_cm1"],
                        "quality_counts": an["quality"].get("counts"), "symmetry_group": an["symmetry"].get("group"),
                        "symmetry_status": an["symmetry"].get("status"),
                        "variational_parse_check_cm1": an["symmetry"].get("variational_parse_check_max_abs_cm-1"),
                        "resonances": an["resonance_counts"],
                        "screen": {"delta_cm1": sc["delta_cm1"], "n_candidates": sc["n_candidates"],
                                   "gaussian_fermi": len(true), "recovered": len(true & keys)}})
        rec["gvpt2"].append(row)
        print("gvpt2", Path(log).name, flush=True)
    write(rec, a.output_dir)


def write(rec, output_dir=None):
    stamp = dt.date.today().strftime("%Y%m%d")
    folder = Path(output_dir) if output_dir is not None else ROOT / "validation"
    out = folder / f"validation_{stamp}_v01.json"
    k = 1
    while out.exists():
        k += 1
        out = folder / f"validation_{stamp}_v{k:02d}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=1, ensure_ascii=False), encoding="utf-8")
    report = out.with_suffix(".md") if output_dir is not None else ROOT / "docs" / "validation.md"
    report.write_text(markdown(rec, out.name, external_output=output_dir is not None), encoding="utf-8")
    print("written", out)


def _veda_table(rows):
    L = ["| Case | Source | Atoms | Modes (imag.) | TED entries within 1 pp | TED within 0.5 pp | Max TED diff (pp) | "
         "Max after summing near-degenerate modes (pp) | PED entries within 1 pp | Max PED diff (pp) |",
         "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for v in rows:
        L.append(f"| {v['case']} | {v['source']} | {v['n_atoms']} | {v['n_modes']} ({v['n_imaginary']}) | "
                 f"{100 * v['ted_fraction_within_1pp']:.2f}% | {100 * v['ted_fraction_within_0p5pp']:.2f}% | "
                 f"{v['ted_max_abs_diff_pp']:.2f} | {v['ted_cluster_summed_max_abs_diff_pp']:.2f} | "
                 f"{100 * v['ped_fraction_within_1pp']:.2f}% | {v['ped_max_abs_diff_pp']:.2f} |")
    return L


def markdown(rec, src, external_output=False):
    record_path = src if external_output else f"validation/{src}"
    L = ["# Validation", "",
         f"Generated by `scripts/run_validation.py` on {rec['date'][:10]} with PEDRO {rec['pedro']} "
         f"(Python {rec['python']}, NumPy {rec['numpy']}). Raw record: `{record_path}`. Default fixtures are owned "
         "calculations included in the repository. Optional VEDA examples are read from a local installation and "
         "not redistributed; optional large ONIOM logs are external inputs intended for the release data archive "
         "(its public link is pending). Input files are identified by SHA-256 in the record.", ""]
    L += ["## 1. Agreement with VEDA on identical input", "",
          "Same Hessian, same masses and the same `.dd2` coordinate set, run through VEDA and through PEDRO. VEDA prints "
          "whole-number percentages, so ±0.5 percentage points (pp) is the best agreement any entry can show.", ""]
    if rec["veda_own"]:
        L += ["### 1a. VEDA 4.99e on our own calculations (`validation/veda_reference/`)", "",
              "Methanol with the coordinate sets published by Jamróz (2013, Tables 1A/1B); ethylene and water with the "
              "sets VEDA builds itself (*Create DD2*). VEDA was run by hand: input, *Create DD2* where stated, *Calculate*; "
              "no MIX or Optimize.", ""] + _veda_table(rec["veda_own"]) + [""]
    if rec["veda_examples"]:
        L += ["### 1b. VEDA's own example files (VEDA 4.99f, from the VEDA download)", "",
              "Dimethylamine (a1-91) and neopentane (four coordinate sets, with combinations of up to 12 terms), as "
              "distributed with VEDA by its author: his Hessians, his coordinate sets and his VEDA output.", ""] + \
             _veda_table(rec["veda_examples"]) + \
             ["", "Neopentane has many exactly degenerate modes, whose split between partners is arbitrary; the larger "
              "differences sit in those sets and shrink when degenerate partners are summed.", ""]
    if rec["published"]:
        L += ["## 2. Agreement with published values", ""]
        for p in rec["published"]:
            L += [f"### {p['title']}", "", f"*Published*: {p['reference']} ({p['published_level']}). *Ours*: "
                  f"{p['our_calculation']}, with the published coordinate sets transcribed (`validation/published/{p['name']}/`).", ""]
            for t in p["tables"]:
                L += [f"**{t['table']}.** {t['within_1pp']} of {t['n']} published entries reproduced within 1 pp "
                      f"(largest difference {t['max_pp']:.1f} pp); frequencies within {t['freq_max_abs_diff_cm1']:.1f} cm⁻¹."
                      + ("" if t["pedro_vs_veda_on_ours_max_pp"] is None else
                         f" VEDA run on our Hessian with the same set agrees with PEDRO within {t['pedro_vs_veda_on_ours_max_pp']:.1f} pp "
                         "in every entry, so the remaining differences from the paper come from the underlying calculation."), "",
                      "| s | Coordinate | Published (cm⁻¹) | Ours (cm⁻¹) | PED published | PED VEDA on ours | PED PEDRO |",
                      "|---|---|---:|---:|---:|---:|---:|"]
                for x in t["rows"]:
                    vo = "" if x["veda_on_ours"] is None else f"{x['veda_on_ours']:.0f}"
                    L.append(f"| {x['coord']} | {x['label']} | {x['published_cm1']} | {x['ours_cm1']} | {x['published']} | "
                             f"{vo} | {x['pedro']} |")
                L.append("")
    if rec["logs"]:
        L += ["## 3. Frequencies recomputed from the Hessian vs the program's printed values", "",
              "| File | Program | Atoms | Basis | Max abs. difference (cm⁻¹) | Imaginary modes | Masses |", "|---|---|---:|---|---:|---:|---|"]
        for v in rec["logs"]:
            fq = "" if v["freq_max_abs_diff_cm1"] is None else f"{v['freq_max_abs_diff_cm1']:.4f}"
            L.append(f"| `{v['file']}` | {v['program']} | {v['n_atoms']} | {v['case']} | {fq} | {v['n_imaginary']} | {v['masses']} |")
        L.append("")
    if rec["oniom"]:
        L += ["## 4. Large ONIOM models", "",
              "| File | Atoms | Case | Modes (imag.) | Subsystem atoms | Freq. vs printed (cm⁻¹) | TED closure (pp) | "
              "Identity check (pp) | Time (s) | Peak traced memory (MB) |", "|---|---:|---|---:|---:|---:|---:|---:|---:|---:|"]
        for v in rec["oniom"]:
            for case, r in v["runs"].items():
                fq = "" if r["freq_max_abs_diff_cm1"] is None else f"{r['freq_max_abs_diff_cm1']:.4f}"
                L.append(f"| `{v['file']}` | {v['n_atoms']} | {case} | {r['n_modes']} ({r['n_imaginary']}) | "
                         f"{r['n_subsystem_atoms']} | {fq} | {r['ted_closure_max_pp']:.0e} | {r['identity_max_dev_pp']:.1e} | "
                         f"{r['seconds']} | {r['peak_MB']} |")
        L += ["", "Case F identity values are approximate by construction (docs/theory.md §7); Fprime and P values "
              "test an exact identity. Case P uses the ONIOM high layer as the mobile set.", ""]
    if rec["gvpt2"]:
        L += ["## 5. Resonances and overtones (GVPT2 logs)", "",
              "| File | Atoms | States (fund./overt./comb.) | Mapping residual (cm⁻¹) | Verdicts clean/warning/failed | "
              "Audit flags (anomalous/quartic/mixed/coupled) | Point group | Variational parse check (cm⁻¹) | Time (s) |",
              "|---|---:|---|---:|---|---|---|---:|---:|"]
        for v in rec["gvpt2"]:
            if v["status"] not in ("ok", "incomplete"):
                L.append(f"| `{v['file']}` | {v['n_atoms']} | {v['status']}: {v['reason']} | | | | | | {v['seconds']} |")
                continue
            c, vd, q = v["counts"], v["verdicts"], v["quality_counts"] or {}
            pc = v["variational_parse_check_cm1"]
            L.append(f"| `{v['file']}` | {v['n_atoms']} | {c['fundamental']}/{c['overtone']}/{c['combination']} | "
                     f"{v['mapping_residual_cm1']:.4f} | {vd['clean']}/{vd['warning']}/{vd['failed']} | "
                     f"{q.get('ANOMALOUS_SHIFT')}/{q.get('QUARTIC_OUTLIER')}/{q.get('RESONANCE_MIXED')}/{q.get('COUPLED_TO_FLAGGED')} | "
                     f"{v['symmetry_group'] or v['symmetry_status']} | {'' if pc is None else f'{pc:.3f}'} | {v['seconds']} |")
        sc = [v for v in rec["gvpt2"] if v.get("screen")]
        if sc:
            L += ["", "**How well does the harmonic screen predict Fermi resonances?** The same logs were screened as if only "
                  "the harmonic data existed, and the candidates were compared with the Fermi resonances Gaussian actually "
                  "treated:", "", "| File | Δ (cm⁻¹) | Candidates | Fermi resonances treated by Gaussian | Contained in the screen |",
                  "|---|---:|---:|---:|---:|"]
            for v in sc:
                s_ = v["screen"]
                L.append(f"| `{v['file']}` | {s_['delta_cm1']:g} | {s_['n_candidates']} | {s_['gaussian_fermi']} | "
                         f"{s_['recovered']} ({100 * s_['recovered'] / max(s_['gaussian_fermi'], 1):.0f}%) |")
            L += ["", "The screen is a pointer to where a GVPT2 calculation is needed, not a substitute for one."]
        L.append("")
    L += ["## What this does and does not show", "",
          "- It shows that PEDRO's algebra reproduces VEDA's tables on identical input, that it reproduces published PED "
          "values to the extent the underlying calculation does, and that it recomputes each program's frequencies "
          "from its Hessian.",
          "- It does not validate any structure, coordinate choice or band assignment. Those remain scientific "
          "judgments; PEDRO reports the evidence (forces, imaginary modes, coordinate set) needed to make them.", ""]
    return "\n".join(L)


if __name__ == "__main__":
    main()
