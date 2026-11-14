"""Output folder: CSV tables, NumPy arrays, coordinate files, a manifest and a self-contained HTML report.

A new folder is created for every run (`<input stem>_PEDRO_<date>_vNN`); existing folders are never overwritten.
"""

from __future__ import annotations

import base64
import csv
import datetime as dt
import html
import io
import json
import platform
from pathlib import Path

import numpy as np

from . import __version__
from . import coordinates as co


def _new_dir(res, out_dir):
    if out_dir is not None:
        p = Path(out_dir)
        if p.exists() and any(p.iterdir()):
            raise FileExistsError(f"{p} exists and is not empty; PEDRO never overwrites results")
        p.mkdir(parents=True, exist_ok=True)
        return p
    src = Path(res.job.source)
    stamp = dt.date.today().strftime("%Y%m%d")
    for v in range(1, 1000):
        p = src.parent / f"{src.stem}_PEDRO_{stamp}_v{v:02d}"
        if not p.exists():
            p.mkdir()
            return p
    raise RuntimeError("too many result folders")


def _csv(path, header, rows):
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def _jsonable(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))


def environment_identity_check(checks):
    """Display the case-specific check without treating projected F as an exact identity."""
    if "identity_F_max_dev_pp" in checks:
        return ("Environment TED vs amplitude (case F; approximate)",
                f"max deviation {checks['identity_F_max_dev_pp']:.4g} pp; "
                "projected accounting includes relative rigid motion of the subsystem in the environment coordinates")
    if "identity_max_dev_pp" in checks:
        return ("Environment TED vs amplitude (exact identity)",
                f"max deviation {checks['identity_max_dev_pp']:.4g} pp; "
                "subsystem rigid-body TED is accounted for separately")
    return None


def _partition_columns(res):
    """Non-internal TED totals; rigid motion of S belongs to S, not to E."""
    ni = res.basis.n_internal
    if res.basis.case in ("Fprime", "P"):
        return {"subsystem_rigid_body_total": res.ted[:, ni:ni + 6].sum(1),
                "environment_total": res.ted[:, ni + 6:].sum(1)}
    if res.basis.case == "F":
        return {"environment_total": res.ted[:, ni:].sum(1)}
    return {}


def manifest(res):
    # pip records the exact Git commit for managed installs from a VCS URL, including beta tags.
    from importlib import metadata
    installation = None
    try:
        direct_url = metadata.distribution("pedro-ped").read_text("direct_url.json")
        installation = json.loads(direct_url) if direct_url else None
    except metadata.PackageNotFoundError:
        pass
    return {
        "software": "PEDRO", "version": __version__, "created": dt.datetime.now().isoformat(timespec="seconds"),
        "input": res.job.summary(), "options": res.options.to_dict(), "checks": res.checks,
        "coordinates": res.coord_report, "basis": {"case": res.basis.case, "n_coordinates": len(res.basis.labels),
                                                   "n_internal": res.basis.n_internal,
                                                   "subsystem_atoms_1based": (res.basis.subsystem + 1).tolist()},
        "timings_s": res.timings, "validation": res.validation,
        "installation": installation,
        "anharmonic": res.checks.get("anharmonic"),
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "scipy": __import__("scipy").__version__, "platform": platform.platform()},
        "units": {"frequency": "cm^-1 (unscaled; 'scaled' columns use the display factor only)",
                  "ir": "km/mol", "distributions": "percent, signed, never clipped or renormalized"},
        "ted_subsystem_columns": {
            "internal_coordinates": "signed TED of the subsystem internal coordinates",
            "environment_total": ("projected environment-coordinate TED; only approximately equal to environment amplitude"
                                  if res.basis.case == "F" else
                                  "TED of mobile environment Cartesians, excluding subsystem rigid motion")
                                  if res.basis.case != "isolated" else None,
            "subsystem_rigid_body_total": ("TED of the subsystem's three translations and three rotations"
                                           if res.basis.case in ("Fprime", "P") else None),
        },
    }


def write_all(res, out_dir=None):
    out = _new_dir(res, out_dir)
    m = res.modes
    ni = res.basis.n_internal
    labels = res.basis.labels
    rows = res.mode_table()
    _csv(out / "modes.csv", list(rows[0].keys()), [list(r.values()) for r in rows])
    head = ["mode", "freq_cm1"] + [f"{res.basis.internal[i].id or 'S%d' % (i + 1)} {labels[i]}" for i in range(ni)]
    partitions = _partition_columns(res)
    _csv(out / "ted_subsystem.csv", head + list(partitions),
         [[k + 1, f"{m.freqs[k]:.4f}"] + [f"{v:.3f}" for v in res.ted[k, :ni]]
          + [f"{values[k]:.3f}" for values in partitions.values()] for k in range(m.n)])
    groups = list(res.group_totals)
    _csv(out / "group_totals.csv", ["mode", "freq_cm1"] + groups,
         [[k + 1, f"{m.freqs[k]:.4f}"] + [f"{res.group_totals[g][k]:.3f}" for g in groups] for k in range(m.n)])
    if res.veda_ped is not None:
        _csv(out / "ped_veda_compatible.csv", head,
             [[k + 1, f"{m.freqs[k]:.4f}"] + [f"{v:.3f}" for v in res.veda_ped[k, :ni]] for k in range(m.n)])
    if res.ped_diag is not None:
        _csv(out / "ped_morino_kuchitsu.csv", head,
             [[k + 1, f"{m.freqs[k]:.4f}"] + ["" if np.isnan(v) else f"{v:.3f}" for v in res.ped_diag[k, :ni]]
              for k in range(m.n)])
    np.savez_compressed(out / "ted.npz", ted=res.ted, ted_cluster_summed=res.ted_cluster, freqs_cm1=m.freqs,
                        ir_km_mol=m.ir_km_mol if m.ir_km_mol is not None else np.array([]),
                        labels=np.array(labels), groups=np.array(res.basis.groups), cluster_id=res.cluster_id)
    (out / "coordinates.json").write_text(json.dumps(co.export_json(res.basis.internal, res.job.elements), indent=1,
                                                     ensure_ascii=False), encoding="utf-8")
    co.export_dd2(res.basis.internal, res.job.elements, res.job.xyz_bohr, out / "coordinates.dd2")
    write_anharmonic(res, out)
    (out / "manifest.json").write_text(json.dumps(manifest(res), indent=1, default=_jsonable, ensure_ascii=False),
                                       encoding="utf-8")
    (out / "summary.txt").write_text(summary_text(res), encoding="utf-8")
    (out / "report.html").write_text(html_report(res), encoding="utf-8")
    return out


STATE_COLS = ["kind", "pedro_quanta", "gaussian_quanta", "e_harm_cm1", "e_anharm_cm1", "shift_cm1", "i_harm_km_mol",
              "i_anharm_km_mol", "verdict", "flags", "resonances", "projection_pct", "overlap_indicator", "irrep",
              "nu_symmetry_filtered_cm1", "contamination_cm1", "parentage"]


def write_anharmonic(res, out):
    a = res.anharmonic
    if a and a.get("states"):
        rows = sorted(a["states"], key=lambda s: -s["e_anharm_cm1"])
        _csv(out / "anharmonic_states.csv", STATE_COLS,
             [["; ".join(r[c]) if isinstance(r[c], list) else ("" if r[c] is None else r[c]) for c in STATE_COLS]
              for r in rows])
        audit = {k: v for k, v in a.items() if k != "states"}
        (out / "gvpt2_audit.json").write_text(json.dumps(audit, indent=1, default=_jsonable, ensure_ascii=False),
                                              encoding="utf-8")
    sc = res.resonance_screen
    if sc:
        if sc["candidates"]:
            cols = list(sc["candidates"][0])
            _csv(out / "resonance_screen.csv", cols, [[r[c] for c in cols] for r in sc["candidates"]])
        if sc["overtones"]:
            cols = list(sc["overtones"][0])
            _csv(out / "overtones_harmonic.csv", cols, [[r[c] for c in cols] for r in sc["overtones"]])


def shown_states(a):
    """Fundamentals plus overtones/combinations with I(anharm) >= 1 km/mol, highest energy first."""
    return sorted([s for s in a["states"] if s["kind"] == "fundamental" or (s["i_anharm_km_mol"] or 0) >= 1.0],
                  key=lambda s: -s["e_anharm_cm1"])


def anharmonic_lines(res):
    a, sc = res.anharmonic or {}, res.resonance_screen
    L = ["", "RESONANCES AND OVERTONES"]
    if a.get("states"):
        L.append(f"GVPT2 ({a.get('pt2_model')}): {a['counts']['fundamental']} fundamentals, {a['counts']['overtone']} "
                 f"overtones, {a['counts']['combination']} combinations; verdicts {a['verdicts']}; resonances "
                 f"{a['resonance_counts']}; symmetry gauge: {a['symmetry'].get('group') or a['symmetry'].get('status')}")
        L += [f"Note: {n}" for n in a.get("notes", [])]
        L.append(a["caveat"])
        L.append(f"{'state':>16} {'E(harm)':>9} {'E(anh)':>9} {'I(anh)':>8} {'verdict':>8}  parentage | resonances | flags")
        for s in shown_states(a):
            L.append(f"{s['pedro_quanta']:>16} {s['e_harm_cm1']:9.1f} {s['e_anharm_cm1']:9.1f} {s['i_anharm_km_mol']:8.2f} "
                     f"{s['verdict']:>8}  {s['parentage']} | {'; '.join(s['resonances']) or '-'} | "
                     f"{'; '.join(s['flags']) or '-'}")
        L.append("(fundamentals and states with I(anh) >= 1 km/mol; all states: anharmonic_states.csv; "
                 "audit and symmetry gauge: gvpt2_audit.json)")
    else:
        L.append(f"GVPT2: {a.get('status')} ({a.get('reason')})")
    if sc:
        L.append(f"Harmonic screen: {sc['n_candidates']} Fermi-resonance candidates (|w_k - w_i - w_j| <= {sc['delta_cm1']:g} "
                 f"cm-1) for {sc['n_fundamentals_screened']} IR-active fundamentals among {sc['n_modes_screened']} modes "
                 f"with >= {sc['min_subsystem_amplitude_pct']:g}% amplitude on the subsystem.")
        L.append(sc["caveat"])
        for r in sc["candidates"][:30]:
            L.append(f"  nu{r['fundamental']} {r['omega_k_cm1']:8.1f} ~ {r['partner']:>12} {r['partner_harmonic_cm1']:8.1f}  "
                     f"delta {r['delta_cm1']:+6.1f}  [{r['fundamental_character']}]")
        L.append("All candidates: resonance_screen.csv; harmonic overtones: overtones_harmonic.csv")
    return L


def anharmonic_html(res):
    esc = html.escape
    a, sc = res.anharmonic or {}, res.resonance_screen
    parts = ["<h2>Resonances and overtones</h2>"]
    col = {"clean": "#1d7a3a", "warning": "#a86b00", "failed": "#a33"}
    if a.get("states"):
        parts.append(f"<p>GVPT2 ({esc(str(a.get('pt2_model')))}): {a['counts']} &middot; verdicts {a['verdicts']} &middot; "
                     f"symmetry gauge {esc(str(a['symmetry'].get('group') or a['symmetry'].get('status')))}. "
                     f"{esc(a['caveat'])}</p>")
        parts.append("<table><tr><th>State</th><th>E(harm)</th><th>E(anh)</th><th>I(anh)</th><th>Verdict</th>"
                     "<th>Parentage (TED of parents)</th><th>Resonances</th><th>Flags</th></tr>")
        for s in shown_states(a):
            parts.append(f"<tr><td>{esc(s['pedro_quanta'])}</td><td>{s['e_harm_cm1']:.1f}</td><td>{s['e_anharm_cm1']:.1f}</td>"
                         f"<td>{s['i_anharm_km_mol']:.2f}</td><td style=\"color:{col[s['verdict']]}\">{s['verdict']}</td>"
                         f"<td>{esc(s['parentage'])}</td><td>{esc('; '.join(s['resonances']))}</td>"
                         f"<td>{esc('; '.join(s['flags']))}</td></tr>")
        parts.append("</table><p class=\"note\">Fundamentals and states with I(anh) &ge; 1 km/mol; all states in "
                     "anharmonic_states.csv.</p>")
    else:
        parts.append(f"<p>GVPT2: {esc(str(a.get('status')))} ({esc(str(a.get('reason')))})</p>")
    if sc:
        parts.append(f"<p>Harmonic screen: {sc['n_candidates']} Fermi-resonance candidates. {esc(sc['caveat'])}</p>")
        if sc["candidates"]:
            parts.append("<table><tr><th>Fundamental</th><th>cm⁻¹</th><th>Partner</th><th>cm⁻¹</th><th>Δ</th>"
                         "<th>Fundamental character</th><th>Partner character</th></tr>")
            for r in sc["candidates"][:40]:
                parts.append(f"<tr><td>ν{r['fundamental']}</td><td>{r['omega_k_cm1']:.1f}</td><td>{esc(r['partner'])}</td>"
                             f"<td>{r['partner_harmonic_cm1']:.1f}</td><td>{r['delta_cm1']:+.1f}</td>"
                             f"<td>{esc(r['fundamental_character'])}</td>"
                             f"<td>{esc(r['partner_character'])}</td></tr>")
            parts.append("</table>")
    return "".join(parts)


# ----------------------------------------------------------------------------- text summary
def summary_text(res):
    j, c, m = res.job, res.checks, res.modes
    L = [f"PEDRO {__version__}   {dt.datetime.now():%Y-%m-%d %H:%M}",
         f"Input       {j.source}",
         f"            {j.natoms} atoms, {j.formula}, charge {j.charge}, multiplicity {j.multiplicity}",
         f"Route       {j.route}",
         f"Masses      {j.mass_source}",
         f"Basis       case {c['case']}: {res.basis.n_internal} internal coordinates of {len(res.basis.subsystem)} "
         f"subsystem atoms ({c['subsystem_choice']}); {len(res.basis.labels)} coordinates in total",
         f"Modes       {m.n} ({c['n_imaginary']} imaginary)"]
    for f in c.get("fragments") or []:
        L.append(f"Fragment    {f['n']:4d} atoms  {f['formula']:<14s} {f['tag']:<8s} atoms {f['ranges']}")
    if c.get("reconstruction"):
        L.append(f"Check       max |PEDRO - printed| frequency = {c['reconstruction']['max_abs_cm1']:.4f} cm-1")
    L.append(f"Check       TED closure (max |row sum - 100|) = {c['ted_closure_max_pp']:.1e} pp; cond(L) = {c['condition_L']:.3g}")
    identity = environment_identity_check(c)
    if identity:
        L.append(f"Check       {identity[0]}: {identity[1]}")
    if "stationarity" in c:
        L.append(f"Check       max atomic force = {c['stationarity']['max_atom_force_Eh_bohr']:.2e} Eh/bohr")
    if res.validation and res.validation.get("status") == "compared":
        v = res.validation
        L.append(f"VEDA        {100 * v['ted_fraction_within_1pp']:.2f}% of TED entries within 1 pp of VEDA {v['veda_version']}; "
                 f"max difference {v['ted_max_abs_diff_pp']:.1f} pp ({v['ted_cluster_summed_max_abs_diff_pp']:.1f} pp after "
                 "summing near-degenerate modes)")
    L += ["", f"Scale factor for display: {res.options.scale}", "",
          f"{'mode':>5} {'freq':>10} {'scaled':>10} {'IR':>9} {'S-TED%':>7}  description (TED %, |value| >= "
          f"{res.options.ted_threshold:g})"]
    ir = m.ir_km_mol
    for k in range(m.n - 1, -1, -1):
        L.append(f"{k + 1:5d} {m.freqs[k]:10.2f} {res.scaled[k]:10.2f} {'' if ir is None else f'{ir[k]:9.2f}':>9} "
                 f"{res.subsystem_ted[k]:7.1f}  {res.descriptions[k]}")
    L += anharmonic_lines(res)
    return "\n".join(L) + "\n"


# ----------------------------------------------------------------------------- HTML report
def _fig_png(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    import matplotlib.pyplot as plt
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def figures(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    m = res.modes
    out = {}
    f = res.scaled
    ir = m.ir_km_mol if m.ir_km_mol is not None else np.ones(m.n)
    fig, ax = plt.subplots(figsize=(9, 3.2))
    grid = np.linspace(max(0, f.min() - 50), f.max() + 50, 4000)
    spec = sum(i * (10 / np.pi) / ((grid - x) ** 2 + 100) for x, i in zip(f[f > 0], ir[f > 0]))
    ax.plot(grid, spec / max(spec.max(), 1e-12) * ir.max(), color="#888", lw=0.8)
    sc = ax.scatter(f, ir, c=np.clip(res.subsystem_ted, 0, 100), cmap="viridis", s=14, vmin=0, vmax=100, zorder=3)
    ax.vlines(f, 0, ir, color="#bbb", lw=0.6)
    ax.set_xlabel(f"wavenumber / cm$^{{-1}}$ (×{res.options.scale:g})")
    ax.set_ylabel("IR / km mol$^{-1}$")
    ax.invert_xaxis()
    fig.colorbar(sc, ax=ax, label="subsystem TED %")
    out["spectrum"] = _fig_png(fig)
    # group totals for the strongest-IR modes of the subsystem
    sel = np.argsort(-(ir * np.clip(res.subsystem_amplitude, 0, 100)))[:25]
    sel = np.sort(sel)
    groups = list(res.group_totals)
    fig, ax = plt.subplots(figsize=(9, 3.6))
    bottom_p = np.zeros(sel.size)
    bottom_n = np.zeros(sel.size)
    cmap = plt.get_cmap("tab20")
    for gi, g in enumerate(groups):
        v = res.group_totals[g][sel]
        pos, neg = np.clip(v, 0, None), np.clip(v, None, 0)
        ax.bar(range(sel.size), pos, bottom=bottom_p, color=cmap(gi % 20), label=g)
        ax.bar(range(sel.size), neg, bottom=bottom_n, color=cmap(gi % 20))
        bottom_p += pos
        bottom_n += neg
    ax.set_xticks(range(sel.size), [f"{res.scaled[k]:.0f}" for k in sel], rotation=90, fontsize=7)
    ax.set_ylabel("TED %")
    ax.axhline(0, color="k", lw=0.5)
    ax.legend(fontsize=7, ncol=3, loc="upper left", bbox_to_anchor=(1.01, 1))
    ax.set_title("Group totals for the 25 most IR-intense subsystem modes", fontsize=9)
    out["groups"] = _fig_png(fig)
    return out


def html_report(res):
    j, c = res.job, res.checks
    try:
        figs = figures(res)
    except Exception as e:  # plotting must never block the numerical results
        figs = {"error": str(e)}
    esc = html.escape
    rows = []
    ir = res.modes.ir_km_mol
    order = np.argsort(-res.modes.freqs)
    for k in order:
        if res.basis.case != "isolated" and res.subsystem_amplitude[k] < 20:
            continue
        cls = ' class="imag"' if res.modes.freqs[k] < 0 else ""
        rows.append(f"<tr{cls}><td>{k + 1}</td><td>{res.modes.freqs[k]:.2f}</td><td>{res.scaled[k]:.2f}</td>"
                    f"<td>{'' if ir is None else f'{ir[k]:.1f}'}</td><td>{res.subsystem_ted[k]:.0f}</td>"
                    f"<td>{res.subsystem_amplitude[k]:.0f}</td><td>{esc(res.descriptions[k])}</td></tr>")
    rec = c.get("reconstruction") or {}
    val = res.validation or {}
    checks = [("Frequencies reproduce the printed values", f"max difference {rec['max_abs_cm1']:.4f} cm⁻¹" if rec else "no printed frequencies"),
              ("Imaginary modes (kept, flagged)", str(c["n_imaginary"])),
              ("TED closure", f"{c['ted_closure_max_pp']:.1e} percentage points"),
              ("Condition number of L", f"{c['condition_L']:.3g}")]
    identity = environment_identity_check(c)
    if identity:
        checks.append(identity)
    if "stationarity" in c:
        checks.append(("Largest atomic force", f"{c['stationarity']['max_atom_force_Eh_bohr']:.2e} Eh/bohr"))
    if val.get("status") == "compared":
        checks.append((f"Agreement with VEDA {val['veda_version']}",
                       f"{100 * val['ted_fraction_within_1pp']:.2f}% of TED entries within 1 pp; "
                       f"max {val['ted_max_abs_diff_pp']:.1f} pp ({val['ted_cluster_summed_max_abs_diff_pp']:.1f} pp cluster-summed)"))
    chk = "".join(f"<tr><th>{esc(a)}</th><td>{esc(b)}</td></tr>" for a, b in checks)
    note = ("" if res.basis.case == "isolated" else
            "<p class='note'>Only modes with at least 20% of their mass-weighted amplitude on the subsystem are listed; "
            "all modes are in modes.csv.</p>")
    imgs = "".join(f"<img src='{v}' alt='{k}'>" for k, v in figs.items() if k != "error")
    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><title>PEDRO report: {esc(Path(j.source).name)}</title>
<style>
:root{{--ink:#1b2430;--mute:#5b6675;--line:#d9dee5;--accent:#0f6e74;--bad:#a33}}
body{{font:14px/1.5 system-ui,Segoe UI,sans-serif;color:var(--ink);max-width:1100px;margin:24px auto;padding:0 16px;background:#fff}}
h1{{font-size:22px;margin:0 0 4px}} h2{{font-size:16px;margin:28px 0 8px;color:var(--accent)}}
.meta{{color:var(--mute)}} table{{border-collapse:collapse;width:100%;font-size:13px}}
td,th{{border-bottom:1px solid var(--line);padding:4px 8px;text-align:left;vertical-align:top}}
th{{font-weight:600}} tr.imag td{{color:var(--bad)}} img{{max-width:100%;margin:8px 0}} .note{{color:var(--mute)}}
</style></head><body>
<h1>PEDRO {__version__}: {esc(Path(j.source).name)}</h1>
<div class="meta">{j.natoms} atoms · {esc(j.formula)} · charge {j.charge} · multiplicity {j.multiplicity} · {esc(j.route)}<br>
Basis case <b>{esc(c['case'])}</b>: {res.basis.n_internal} internal coordinates on {len(res.basis.subsystem)} subsystem atoms
({esc(c['subsystem_choice'])}) · display scale {res.options.scale:g}</div>
<h2>Checks</h2><table>{chk}</table>
<h2>Spectrum</h2>{imgs}
<h2>Mode characters (TED %)</h2>{note}
<table><tr><th>#</th><th>cm⁻¹</th><th>scaled</th><th>IR</th><th>S-TED %</th><th>S-amp %</th><th>Description</th></tr>
{''.join(rows)}</table>
{anharmonic_html(res)}
<p class="note">TED = total energy distribution, T<sub>ik</sub> = 100·L<sub>ik</sub>(L⁻¹)<sub>ki</sub>. Values are signed and
sum to 100 for every mode. See docs/theory.md for definitions.</p>
</body></html>"""
