"""GVPT2 quality audit (T6, dev18): which Gaussian GVPT2 states can be trusted?

Gaussian's deperturbed/variational VPT2 output contains everything needed to judge each state, but prints it in
separate tables. This module reads them from the log text of one Freq=Anharmonic step and returns per-mode flags.
Mode numbers are Gaussian's ANHARMONIC numbering (descending harmonic frequency).

Tables used: 'Fundamental Bands' (omega, final nu), 'Vibrational Energies' (deperturbed nu_DVPT2 and value after
the variational treatment), 'Projection of DVPT2 states on New Variational States' (overlaps < 80 %), 'Variational
Matrix Definition' (which states are coupled), the 'Total / 3rd Deriv. / 4th Deriv. contributions to X Matrix', and
'QUARTIC FORCE CONSTANTS IN NORMAL MODES' (reduced, cm-1).

Flags (thresholds are module constants; every flag carries its numbers):
  ANOMALOUS_SHIFT     400 < omega < 2500 cm-1 and a_DVPT2 = nu_DVPT2/omega outside [A_LOW, A_HIGH]; or
                      omega <= 400 and |a_DVPT2 - 1| > SOFT_TOL. Position unusable.
  QUARTIC_OUTLIER     the mode shares a semi-diagonal quartic constant |k_iijj| > K_OUTLIER cm-1 with a partner,
                      neither an X-H stretch (omega < XH_OMEGA). Typical |k_iijj| are a few cm-1; values this large among
                      heavy-atom modes indicate numerical-differentiation artifacts (or, for open shells, vibronic
                      softening that VPT2 cannot treat).
  RESONANCE_MIXED     Gaussian's projection overlap < 80 %: the final state has mixed character.
  CONTAMINATED        (dev21) |1/2 sum_j x_ij| over corrupted partners j exceeds max(5 cm-1, 20 % of the mode's
                      deperturbed shift): a sizeable part of the shift is inherited from modes VPT2 cannot treat.
  COUPLED_TO_FLAGGED  the variational matrix couples the mode to a state containing an ANOMALOUS_SHIFT or
                      QUARTIC_OUTLIER mode: position and (borrowed) intensity may inherit the artifact.
State flags: a fundamental takes its mode's flags; an overtone/combination gets PARENT_FLAGGED[...] if any parent mode
carries ANOMALOUS_SHIFT, QUARTIC_OUTLIER or COUPLED_TO_FLAGGED. dev19: pedro.analysis adds VPT2_BREAKDOWN(...) to a
state that fails PEDRO's own breakdown check (anharmonic.breakdown_problems); with the audit available, those checks
are per-state verdicts and lower the report-level status only to 'warning'.
Validated on: water B2PLYP-D3/aVTZ (no flags) and a C2v test molecule, neutral and anion, wB97X-D/def2-SVP.
"""
from __future__ import annotations

import re

import numpy as np

A_LOW, A_HIGH, SOFT_TOL = 0.93, 1.01, 0.15
K_OUTLIER, XH_OMEGA = 100.0, 2500.0
OVERLAP_PRINT_LIMIT = 80.0
SEVERE = ("ANOMALOUS_SHIFT", "QUARTIC_OUTLIER", "COUPLED_TO_FLAGGED")
CONTAM_ABS, CONTAM_REL = 5.0, 0.20   # dev21: CONTAMINATED if |1/2 sum x_ij (corrupted j)| > max(5 cm-1, 20 % of |nu_DVPT2 - omega|)


def _tri(lines, title, n):
    hits = [k for k, l in enumerate(lines) if title in l]
    if not hits:
        return None
    i, X, cols = hits[0] + 2, np.zeros((n, n)), None
    while i < len(lines):
        l = lines[i]
        if re.match(r"^\s+(\d+\s+)*\d+\s*$", l):
            cols = [int(c) for c in l.split()]
        elif cols and re.match(r"^\s+\d+\s+-?\d", l):
            f = l.split()
            r = int(f[0])
            for c, v in zip(cols, f[1:]):
                if r <= n and c <= n:
                    X[r - 1, c - 1] = X[c - 1, r - 1] = float(v.replace("D", "E"))
        elif cols and (l.strip() == "" or "Matrix" in l or "contributions" in l):
            break
        i += 1
    return X


def audit(text: str) -> dict:
    """Audit one GVPT2 log step. Returns {'status': 'audited'|'unavailable', 'modes': {mode: {...}}, ...}."""
    if "Fundamental Bands" not in text:
        return {"status": "unavailable", "reason": "no GVPT2 'Fundamental Bands' table"}
    i = text.index(" Fundamental Bands")
    fund = {}
    for line in text[i:i + 20000].splitlines()[3:]:
        m = re.match(r"^\s*(\S*)\s+(\d+)\(1\)\s+active\s+([\d.]+)\s+([\d.-]+)", line)
        if not m:
            if fund:
                break
            continue
        fund[int(m.group(2))] = {"omega": float(m.group(3)), "nu_final": float(m.group(4)), "gaussian_row_flag": m.group(1)}
    if not fund:
        return {"status": "unavailable", "reason": "'Fundamental Bands' table not parseable"}
    n = max(fund)
    missing = []
    j = text.find("Vibrational Energies (cm^-1)")
    if j >= 0:
        for line in text[j:j + 40000].splitlines()[3:]:
            m = re.match(r"^\s+(\d+)\(1\)\s+([\d.-]+)\s+([\d.-]+)\s*$", line)
            if m and int(m.group(1)) in fund:
                fund[int(m.group(1))]["nu_depert"] = float(m.group(2))
            elif line.strip() and not re.match(r"^\s+\d+\(\d\)", line) and "(" not in line:
                break
    else:
        missing.append("deperturbed energies")
    k = text.find("Projection of DVPT2 states on New Variational States")
    overlaps = {}
    if k >= 0:
        overlaps = {int(a): float(b) for a, b in re.findall(r"State \|(\d+)\(1\)>\s+has overlap of\s+([\d.]+)%", text[k:k + 40000])}
    couplings = []
    v0, v1 = text.find("Variational Matrix Definition"), k
    if v0 >= 0 and v1 > v0:
        for typ, s1, s2 in re.findall(r"^\s+(1-1|1-2|2-2|1-3)\s+\|\s+([\d() ]+?)\s+\|\s+\|\s+([\d() ]+?)\s*$", text[v0:v1], re.M):
            couplings.append((typ, [int(a) for a in re.findall(r"(\d+)\(\d\)", s1)], [int(a) for a in re.findall(r"(\d+)\(\d\)", s2)]))
    if "deperturbed energies" in missing and not couplings:
        missing.remove("deperturbed energies")       # Gaussian prints that table only when resonances are treated
    lines = text.splitlines()
    Xt, X3, X4 =(_tri(lines, t, n) for t in ("Total Anharmonic X Matrix", "3rd Deriv. contributions to X Matrix",
                                              "4th Deriv. contributions to X Matrix"))
    if X4 is None:
        missing.append("X-matrix contributions")
    semi = {}
    q = text.find("QUARTIC FORCE CONSTANTS IN NORMAL MODES")
    if q >= 0:
        qb = text[q:]
        e = qb.find("Num. of 4th derivatives")
        for line in (qb[:e] if e > 0 else qb[:800000]).splitlines():
            f = line.split()
            if len(f) >= 7 and all(x.isdigit() for x in f[:4]) and f[0] == f[1] and f[2] == f[3] and f[0] != f[2]:
                semi[(int(f[0]), int(f[2]))] = float(f[4])
    else:
        missing.append("quartic force constants")
    heavy = {m for m, d in fund.items() if d["omega"] < XH_OMEGA}
    outlier = {}
    for (a, b), kv in semi.items():
        if abs(kv) > K_OUTLIER and a in heavy and b in heavy:
            outlier.setdefault(a, []).append((b, kv))
            outlier.setdefault(b, []).append((a, kv))
    modes = {}
    for m, d in fund.items():
        dep = d.get("nu_depert", d["nu_final"])
        a_dep = dep / d["omega"]
        flags = []
        if (400 < d["omega"] < XH_OMEGA and not A_LOW <= a_dep <= A_HIGH) or (d["omega"] <= 400 and abs(a_dep - 1) > SOFT_TOL):
            flags.append(f"ANOMALOUS_SHIFT(a_DVPT2={a_dep:.3f})")
        if m in outlier:
            top = sorted(outlier[m], key=lambda t: -abs(t[1]))[:3]
            flags.append("QUARTIC_OUTLIER(" + ",".join(f"k({m},{m},{b},{b})={kv:+.0f}" for b, kv in top) + ")")
        if m in overlaps:
            flags.append(f"RESONANCE_MIXED({overlaps[m]:.0f}%)")
        modes[m] = {**d, "nu_depert": dep, "a_DVPT2": round(a_dep, 5), "a_final": round(d["nu_final"] / d["omega"], 5),
                    "overlap_pct": overlaps.get(m), "flags": flags,
                    "x_rowsum_total": None if Xt is None else round(float(Xt[m - 1].sum()), 2),
                    "x_rowsum_3rd": None if X3 is None else round(float(X3[m - 1].sum()), 2),
                    "x_rowsum_4th": None if X4 is None else round(float(X4[m - 1].sum()), 2)}
    severe = {m for m, d in modes.items() if any(f.startswith(("ANOMALOUS_SHIFT", "QUARTIC_OUTLIER")) for f in d["flags"])}
    for typ, s1, s2 in couplings:
        for a, b in ((s1, s2), (s2, s1)):
            if len(a) == 1 and a[0] in modes:
                bad = sorted(x for x in b if x in severe and x != a[0])
                modes[a[0]].setdefault("coupled_states", []).append({"type": typ, "partner_modes": b})
                if bad and a[0] not in severe:
                    modes[a[0]]["_bad_partners"] = sorted(set(modes[a[0]].get("_bad_partners", [])) | set(bad))
    for m, d in modes.items():                    # dev21: contamination of every mode by corrupted partners
        if Xt is not None:
            c = 0.5 * sum(float(Xt[m - 1, j - 1]) for j in severe if j != m)
            d["contamination_cm-1"] = round(c, 3)
            shift = d["nu_depert"] - d["omega"]
            if m not in severe and abs(c) > max(CONTAM_ABS, CONTAM_REL * abs(shift)):
                d["flags"].append(f"CONTAMINATED(½Σx={c:+.1f} of shift {shift:+.1f} cm⁻¹)")
    for m, d in modes.items():
        if d.get("_bad_partners"):
            d["flags"].append("COUPLED_TO_FLAGGED(" + ",".join(str(x) for x in d.pop("_bad_partners")) + ")")
    med = float(np.median(np.abs(list(semi.values())))) if semi else None
    return {"status": "audited", "mode_numbering": "Gaussian anharmonic numbering (descending harmonic frequency)",
            "missing_tables": missing, "n_modes": len(modes), "median_abs_semidiagonal_quartic_cm-1": med,
            "thresholds": {"A_LOW": A_LOW, "A_HIGH": A_HIGH, "SOFT_TOL": SOFT_TOL, "K_OUTLIER_cm-1": K_OUTLIER,
                           "XH_OMEGA_cm-1": XH_OMEGA, "overlap_print_limit_pct": OVERLAP_PRINT_LIMIT,
                       "CONTAM_ABS_cm-1": CONTAM_ABS, "CONTAM_REL": CONTAM_REL},
            "n_couplings": len(couplings), "modes": {str(m): d for m, d in sorted(modes.items())},
            "counts": {f: sum(any(x.startswith(f) for x in d["flags"]) for d in modes.values())
                       for f in ("ANOMALOUS_SHIFT", "QUARTIC_OUTLIER", "RESONANCE_MIXED", "COUPLED_TO_FLAGGED")}}


def state_flags(quality: dict, quanta_gaussian: dict) -> list[str]:
    """Flags for a state given its quanta in Gaussian anharmonic numbering ({mode: n})."""
    if not quality or quality.get("status") != "audited":
        return []
    modes = quality["modes"]
    qs = {int(k): int(v) for k, v in quanta_gaussian.items()}
    if len(qs) == 1 and next(iter(qs.values())) == 1:
        return list(modes.get(str(next(iter(qs))), {}).get("flags", []))
    bad = []
    for m in sorted(qs):
        kinds = [f.split("(")[0] for f in modes.get(str(m), {}).get("flags", []) if f.startswith(SEVERE)]
        if kinds:
            bad.append(f"{m}:{'+'.join(kinds)}")
    return [f"PARENT_FLAGGED({','.join(bad)})"] if bad else []


def state_status(flags: list[str]) -> str | None:
    """Record-level status implied by flags: 'failed' when the position is unusable (the state or a parent has an
    anomalous deperturbed shift, or the state itself fails PEDRO's VPT2-breakdown check), 'warning' for any other
    flag, None if clean."""
    if not flags:
        return None
    if any("ANOMALOUS_SHIFT" in f or f.startswith("VPT2_BREAKDOWN") for f in flags):
        return "failed"
    return "warning"


def verdict(flags: list[str]) -> str:
    """Per-state verdict shown to users: 'clean', 'warning' or 'failed' (dev19)."""
    return state_status(flags) or "clean"
