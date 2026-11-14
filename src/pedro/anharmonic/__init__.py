"""Resonances and overtones: the GVPT2 layer of PEDRO (docs/theory.md, section 11).

When a Gaussian Freq=Anharmonic (GVPT2) step is present, PEDRO reports every fundamental, two-quantum overtone and
1+1 combination band with Gaussian's anharmonic position and intensity, the Fermi and Darling-Dennison resonances
each state takes part in, Gaussian's variational projection, a quality verdict, and a harmonic-parentage label built
from PEDRO's TED of the parent modes. Positions and intensities are Gaussian's; PEDRO never recomputes, clips or
rescales them. A parentage label is not an anharmonic energy distribution.

Large-model isolation rule (spec section 9.2): this layer runs only when GVPT2 data are present, every failure is
caught and reported as the layer's status, and nothing here can change or slow the harmonic results.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

from . import quality as gq
from . import symmetry as sg

SHIFT_RATIO_LIMIT = 0.25         # |E_anh - E_harm| / E_harm above this marks a fundamental as a VPT2 breakdown
MIN_FRACTION_OF_HARMONIC = 0.5   # an overtone/combination below this fraction of E_harm is a VPT2 breakdown
MAP_TOL_CM = 0.5                 # Gaussian E(harm) vs PEDRO harmonic frequency for a fundamental
AMBIGUOUS_CM = 0.1               # two PEDRO modes closer than this to one Gaussian mode: mapping is ambiguous
STATIONARY_FORCE = 4.5e-4        # Eh/bohr, Gaussian's default maximum-force threshold
_SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


class AnharmonicError(Exception):
    pass


# ----------------------------------------------------------------------------- locating and parsing the GVPT2 step
def gvpt2_step_text(path) -> str | None:
    """Text of the (last) job step that contains GVPT2 output, or None."""
    text = Path(path).read_text(encoding="latin-1")
    k = text.rfind("Anharmonic Infrared Spectroscopy")
    if k < 0:
        k = text.rfind("Second-order Perturbative Anharmonic Analysis")
        if k < 0:
            return None
    starts = [m.end() for m in re.finditer(r"^ (?:Normal termination of Gaussian|Error termination)[^\n]*$", text[:k], re.M)]
    start = starts[-1] if starts else 0
    end = text.find("Normal termination of Gaussian", k)
    return text[start:(len(text) if end < 0 else end)]


def _table(block, header):
    m = re.search(rf"^ {re.escape(header)}\n -+\n(.*?)(?:\n\s*\n|\Z)", block, re.M | re.S)
    if not m:
        return None
    rows = []
    for line in m.group(1).splitlines():
        if "Mode(n)" in line:
            continue
        quanta = [(int(a), int(b)) for a, b in re.findall(r"(\d+)\((\d+)\)", line)]
        nums = re.findall(r"(?<![\d(])-?\d+\.\d+(?:[EeDd][-+]?\d+)?", line)
        if quanta:
            rows.append({"quanta": quanta, "values": [float(v.replace("D", "E")) for v in nums]})
    return rows


def parse(text) -> dict:
    """Gaussian GVPT2 tables of one step (Gaussian's anharmonic mode numbering, descending harmonic frequency)."""
    out = {}
    m = re.search(r"PT2 model:\s*(.+)", text)
    out["pt2_model"] = m.group(1).strip() if m else None
    m = re.search(r"Normal modes will be PRINTED in (\w+) order", text)
    out["printed_order"] = m.group(1) if m else None
    i = text.find("Anharmonic Infrared Spectroscopy")
    if i < 0:
        raise AnharmonicError("no 'Anharmonic Infrared Spectroscopy' section (GVPT2 incomplete or not requested)")
    ir = text[i:]
    fund, ov, comb = (_table(ir, h) for h in ("Fundamental Bands", "Overtones", "Combination Bands"))
    if not fund:
        raise AnharmonicError("anharmonic 'Fundamental Bands' table not parseable")
    states = []
    for kind, rows in (("fundamental", fund), ("overtone", ov or []), ("combination", comb or [])):
        for r in rows:
            v = r["values"]
            if len(v) < 3:
                continue
            states.append({"kind": kind, "quanta": {a: n for a, n in r["quanta"]},
                           "e_harm": v[0], "e_anharm": v[1],
                           "i_harm": v[2] if kind == "fundamental" else None,
                           "i_anharm": v[3] if kind == "fundamental" and len(v) > 3 else v[-1]})
    # high/low overlap indicators ("H"/"L") from 'Vibrational Energies at Anharmonic Level'
    hl = {}
    j = text.find("Vibrational Energies at Anharmonic Level")
    if j >= 0:
        for line in text[j:i].splitlines():
            mm = re.match(r"^ ([HL])\s+((?:\d+\(\d\)\s+)+)", line)
            if mm:
                hl[_key(mm.group(2))] = mm.group(1)
    # projection of deperturbed states on variational states (< 80 % only)
    proj = {}
    k = text.find("Projection of DVPT2 states on New Variational States")
    if k >= 0:
        for st, pct in re.findall(r"State \|([^>]+)>\s+has overlap of\s+([\d.]+)%", text[k:k + 200000]):
            proj[_key(st)] = float(pct)
    fermi, dd22, dd11 = sg.resonance_lists(text)
    out.update({"states": states, "high_low": hl, "projection": proj, "fermi": fermi, "dd22": dd22, "dd11": dd11,
                "resonance_section": "Resonance Analysis" in text})
    return out


def _key(s):
    return tuple(sorted((int(a), int(b)) for a, b in re.findall(r"(\d+)\((\d+)\)", s)))


def quanta_text(q: dict) -> str:
    """{12: 1, 40: 1} -> 'ν₁₂ + ν₄₀'; {33: 2} -> '2ν₃₃' (PEDRO harmonic mode numbers, ascending, from 1)."""
    return " + ".join(("2" if n == 2 else "") + "ν" + str(h).translate(_SUB) for h, n in sorted(q.items()))


# ----------------------------------------------------------------------------- mapping onto PEDRO's harmonic modes
def map_modes(fund_e_harm: dict, freqs) -> tuple[dict, dict, list]:
    """Gaussian anharmonic mode -> PEDRO mode index (0-based), by matching the printed harmonic energies of the
    fundamentals to PEDRO's frequencies (Hungarian assignment). Returns (map, residual_cm1, ambiguous_modes)."""
    g = sorted(fund_e_harm)
    a = np.array([fund_e_harm[m] for m in g])
    f = np.asarray(freqs)
    C = np.abs(a[:, None] - f[None, :])
    r, c = linear_sum_assignment(C)
    mp = {g[i]: int(c[k]) for k, i in enumerate(r)}
    res = {g[i]: float(C[i, c[k]]) for k, i in enumerate(r)}
    ambiguous = []
    for m, h in mp.items():
        near = np.where(np.abs(f - fund_e_harm[m]) <= res[m] + AMBIGUOUS_CM)[0]
        if len(near) > 1:
            ambiguous.append(m)
    return mp, res, sorted(ambiguous)


def renumber(flag: str, mp: dict) -> str:
    """Rewrite mode numbers inside audit/gauge flags from Gaussian's anharmonic numbering to PEDRO's (ν + 1-based)."""
    nu = lambda m: f"ν{mp[int(m)] + 1}" if int(m) in mp else f"g{m}"
    head, _, body = flag.partition("(")
    if head == "COUPLED_TO_FLAGGED":
        body = re.sub(r"\d+", lambda m: nu(m.group(0)), body)
    elif head == "PARENT_FLAGGED":
        body = re.sub(r"(\d+):", lambda m: nu(m.group(1)) + ":", body)
    elif head == "QUARTIC_OUTLIER":
        body = re.sub(r"k\(([\d,]+)\)", lambda m: "k(" + ",".join(nu(x) for x in m.group(1).split(",")) + ")", body)
    elif head == "FORBIDDEN_COUPLING":
        body = re.sub(r"(\d+)\((\d)\)", lambda m: ("2" if m.group(2) == "2" else "") + nu(m.group(1)), body)
    else:
        return flag
    return f"{head}({body}"


def breakdown(kind, eh, ea, label):
    if kind == "fundamental" and eh > 0 and abs(ea - eh) / eh > SHIFT_RATIO_LIMIT:
        return [f"VPT2_BREAKDOWN(shift {100 * (ea - eh) / eh:+.0f}% of {label})"]
    if kind != "fundamental" and ea < MIN_FRACTION_OF_HARMONIC * eh:
        return [f"VPT2_BREAKDOWN({label} at {ea:.1f} < {100 * MIN_FRACTION_OF_HARMONIC:.0f}% of harmonic {eh:.1f})"]
    if not (math.isfinite(ea) and ea > 0):
        return [f"VPT2_BREAKDOWN(non-positive energy {ea})"]
    return []


# ----------------------------------------------------------------------------- the layer
def run(job, modes, descriptions, *, symmetry="auto", size_limit=sg.SIZE_LIMIT) -> dict:
    """Run the GVPT2 layer for `job` (a Gaussian log) on PEDRO's harmonic `modes`. Never raises: failures are
    returned as {'status': 'failed' | 'incomplete', 'reason': ...}."""
    try:
        return _run(job, modes, descriptions, symmetry, size_limit)
    except AnharmonicError as e:
        return {"status": "failed", "reason": str(e)}
    except Exception as e:  # isolation rule: never let this layer break a harmonic run
        return {"status": "failed", "reason": f"unexpected {type(e).__name__}: {e}"}


def _run(job, modes, descriptions, symmetry, size_limit):
    if job.kind != "gaussian_log":
        return {"status": "not_applicable", "reason": "GVPT2 tables are read from the Gaussian log; load the .log file"}
    text = gvpt2_step_text(job.source)
    if text is None:
        return {"status": "not_applicable", "reason": "no GVPT2 (Freq=Anharmonic) step in this log"}
    if modes.method != "full":
        return {"status": "not_applicable", "reason": "GVPT2 states refer to the full projected harmonic modes; "
                                                       "use the isolated or F basis"}
    if modes.n_imaginary:
        raise AnharmonicError(f"{modes.n_imaginary} imaginary harmonic mode(s): VPT2 is defined only at a minimum")
    p = parse(text)
    notes = []
    if p["printed_order"] not in (None, "DESCENDING"):
        notes.append(f"Gaussian printed modes in {p['printed_order']} order; mapping uses harmonic energies, not order")
    fund = {next(iter(s["quanta"])): s["e_harm"] for s in p["states"] if s["kind"] == "fundamental"}
    if len(fund) != modes.n:
        raise AnharmonicError(f"Gaussian lists {len(fund)} fundamentals for {modes.n} harmonic modes: reduced-dimensional "
                              "GVPT2 is not yet supported (mode numbering unconfirmed)")
    mp, res, ambiguous = map_modes(fund, modes.freqs)
    worst = max(res.values())
    if worst > MAP_TOL_CM:
        raise AnharmonicError(f"Gaussian harmonic energies differ from PEDRO's by up to {worst:.3f} cm⁻¹ (> {MAP_TOL_CM}); "
                              "the GVPT2 step does not belong to this Hessian")
    if job.gradient is not None:
        fmax = float(np.max(np.linalg.norm(job.gradient.reshape(-1, 3), axis=1)))
        if fmax > STATIONARY_FORCE:
            notes.append(f"largest atomic force {fmax:.1e} Eh/bohr exceeds {STATIONARY_FORCE:.1e}: VPT2 assumes a stationary point")

    audit = gq.audit(text)
    gauge = {"status": "skipped", "reason": "symmetry gauge switched off"}
    if symmetry != "off" and audit.get("status") == "audited":
        vecs = {m: modes.Lx[:, h].reshape(-1, 3) for m, h in mp.items()}
        from ..constants import SYMBOLS
        Z = [SYMBOLS.index(e) + 1 for e in job.elements]
        try:
            gauge = sg.gauge(text, audit, Z, job.masses, job.xyz_bohr, vecs,
                             size_limit=10**9 if symmetry == "force" else size_limit)
        except Exception as e:
            gauge = {"status": "error", "reason": f"{type(e).__name__}: {e}"}

    # resonance partners per state (Gaussian numbering keys)
    partners: dict[tuple, list[str]] = {}

    def gq_text(q):
        return quanta_text({mp[a] + 1: n for a, n in q.items()})
    for (i, j, k), d, c, _pv in p["fermi"]:
        two = {j: 2} if j == k else {j: 1, k: 1}
        a, b = ((i, 1),), tuple(sorted(two.items()))
        partners.setdefault(a, []).append(f"Fermi with {gq_text(two)} (Δ={d:+.1f}, k={c:+.1f} cm⁻¹)")
        partners.setdefault(b, []).append(f"Fermi with {gq_text({i: 1})} (Δ={d:+.1f}, k={c:+.1f} cm⁻¹)")
    for (i, j, k, l), d, c in p["dd22"]:
        s1 = tuple(sorted(({i: 2} if i == j else {i: 1, j: 1}).items()))
        s2 = tuple(sorted(({k: 2} if k == l else {k: 1, l: 1}).items()))
        partners.setdefault(s1, []).append(f"DD 2-2 with {gq_text(dict(s2))} (K={c:+.1f} cm⁻¹)")
        partners.setdefault(s2, []).append(f"DD 2-2 with {gq_text(dict(s1))} (K={c:+.1f} cm⁻¹)")
    for (i, j), d, c in p["dd11"]:
        partners.setdefault(((i, 1),), []).append(f"DD 1-1 with {gq_text({j: 1})} (K={c:+.1f} cm⁻¹)")
        partners.setdefault(((j, 1),), []).append(f"DD 1-1 with {gq_text({i: 1})} (K={c:+.1f} cm⁻¹)")

    states = []
    for s in p["states"]:
        key = tuple(sorted(s["quanta"].items()))
        hq = {mp[a] + 1: n for a, n in s["quanta"].items()}            # PEDRO modes, 1-based ascending
        label = " + ".join(("2 × " if n == 2 else "") + f"[{descriptions[h - 1]}]" for h, n in sorted(hq.items()))
        flags = gq.state_flags(audit, {str(a): n for a, n in s["quanta"].items()})
        flags += sg.state_flags(gauge, {str(a): n for a, n in s["quanta"].items()})
        flags = [renumber(f, mp) for f in flags]
        flags += breakdown(s["kind"], s["e_harm"], s["e_anharm"], quanta_text(hq))
        if s["kind"] != "fundamental" and key in p["projection"] and not any(f.startswith("RESONANCE_MIXED") for f in flags):
            flags.append(f"RESONANCE_MIXED({p['projection'][key]:.0f}%)")
        amb = sorted(a for a in s["quanta"] if a in ambiguous)
        if amb:
            flags.append("MAPPING_AMBIGUOUS(near-degenerate parent; label may belong to its partner)")
        mode_q = audit.get("modes", {}).get(str(next(iter(s["quanta"])))) if s["kind"] == "fundamental" else None
        states.append({
            "kind": s["kind"], "pedro_quanta": quanta_text(hq), "pedro_modes": {str(h): n for h, n in sorted(hq.items())},
            "gaussian_quanta": " + ".join(f"{a}({n})" for a, n in sorted(s["quanta"].items())),
            "e_harm_cm1": s["e_harm"], "e_anharm_cm1": s["e_anharm"], "shift_cm1": round(s["e_anharm"] - s["e_harm"], 3),
            "i_harm_km_mol": s["i_harm"], "i_anharm_km_mol": s["i_anharm"],
            "projection_pct": p["projection"].get(key), "overlap_indicator": p["high_low"].get(key),
            "resonances": partners.get(key, []), "parentage": label,
            "irrep": (gauge.get("modes", {}).get(str(next(iter(s["quanta"]))), {}).get("irrep")
                      if s["kind"] == "fundamental" and gauge.get("status") == "ok" else None),
            "nu_symmetry_filtered_cm1": (gauge.get("modes", {}).get(str(next(iter(s["quanta"]))), {}).get("nu_symmetry_filtered")
                                         if s["kind"] == "fundamental" and gauge.get("status") == "ok" else None),
            "contamination_cm1": None if mode_q is None else mode_q.get("contamination_cm-1"),
            "flags": flags, "verdict": gq.verdict(flags)})
    counts = {k: sum(s["kind"] == k for s in states) for k in ("fundamental", "overtone", "combination")}
    verdicts = {v: sum(s["verdict"] == v for s in states) for v in ("clean", "warning", "failed")}
    status = "ok"
    if audit.get("status") != "audited" or audit.get("missing_tables"):
        status = "incomplete"
        notes.append("quality audit incomplete: " + (audit.get("reason") or ", ".join(audit.get("missing_tables", []))))
    return {"status": status, "pt2_model": p["pt2_model"], "n_states": len(states), "counts": counts,
            "verdicts": verdicts, "mapping": {"method": "Hungarian match of Gaussian E(harm) to PEDRO frequencies",
                                              "max_residual_cm1": worst, "ambiguous_gaussian_modes": ambiguous,
                                              "gaussian_to_pedro_1based": {str(a): h + 1 for a, h in sorted(mp.items())}},
            "resonance_counts": {"fermi": len(p["fermi"]), "dd22": len(p["dd22"]), "dd11": len(p["dd11"])},
            "states": states, "quality": audit, "symmetry": gauge, "notes": notes,
            "caveat": "Positions and intensities are Gaussian's GVPT2 values. Parentage labels use PEDRO's harmonic TED "
                      "of the parent modes; they are not an anharmonic energy distribution. Mode numbers (ν) are PEDRO's "
                      "harmonic numbering (ascending frequency); the audit and gauge in gvpt2_audit.json keep Gaussian's "
                      "anharmonic numbering, translated by mapping.gaussian_to_pedro_1based."}
