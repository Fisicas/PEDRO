"""Harmonic resonance screen: overtone and combination positions and Fermi-resonance candidates without GVPT2.

For a harmonic frequency job (for example a large ONIOM model, where GVPT2 is not feasible) PEDRO lists where the
two-quantum overtones 2ω_i and 1+1 combinations ω_i + ω_j of the subsystem's modes fall, and which of them nearly
coincide with an IR-active fundamental of the subsystem (|ω_k − ω_i − ω_j| ≤ Δ). Near-coincidence is necessary
for a Fermi resonance but far from sufficient: whether the states mix depends on a cubic force constant and on the
anharmonic shifts, which only an anharmonic calculation provides. Every row is therefore a CANDIDATE.

Measured on the public GVPT2 examples (docs/validation.md): with Δ = 10 cm⁻¹ the screen contains 7 of benzaldehyde's 12
Fermi resonances treated by Gaussian (11 % of its 61 candidates are real) and misses ethylene's single one.
A spatial-overlap ranking was tested and did not separate real resonances from coincidences, so none is offered.
"""

from __future__ import annotations

import numpy as np

DELTA_CM = 10.0            # harmonic near-coincidence window (cm^-1)
MIN_SUBSYSTEM_AMP = 50.0   # % of amplitude on the subsystem for a mode to enter the screen
MIN_IR = 1.0               # km/mol; fundamentals weaker than this are not screened (when intensities are known)


def screen(modes, subsystem_amplitude, descriptions, natoms, delta=DELTA_CM, min_amp=MIN_SUBSYSTEM_AMP,
           window=None, scale=1.0, min_ir=MIN_IR):
    f = modes.freqs
    ir = modes.ir_km_mol
    sel = np.where((f > 0) & (subsystem_amplitude >= min_amp))[0]
    overtones = [{"mode": int(i + 1), "omega_cm1": float(f[i]), "overtone_2omega_cm1": float(2 * f[i]),
                  "scaled_2omega_cm1": float(2 * f[i] * scale), "parent": descriptions[i]} for i in sel]
    fund = sel
    if ir is not None:
        fund = fund[ir[fund] >= min_ir]
    if window is not None:
        fund = fund[(f[fund] * scale >= window[0]) & (f[fund] * scale <= window[1])]
    pair_i, pair_j = np.triu_indices(sel.size)
    I, J = sel[pair_i], sel[pair_j]
    S = f[I] + f[J]
    rows = []
    for k in fund:
        d = f[k] - S
        for n in np.where((np.abs(d) <= delta) & (I != k) & (J != k))[0]:
            i, j = int(I[n]), int(J[n])
            rows.append({"fundamental": int(k + 1), "omega_k_cm1": float(f[k]),
                         "ir_k_km_mol": None if ir is None else float(ir[k]),
                         "partner": f"2ν{i + 1}" if i == j else f"ν{i + 1} + ν{j + 1}",
                         "partner_kind": "overtone" if i == j else "combination",
                         "partner_harmonic_cm1": float(S[n]), "delta_cm1": float(d[n]),
                         "fundamental_character": descriptions[k],
                         "partner_character": (f"2 × [{descriptions[i]}]" if i == j else
                                               f"[{descriptions[i]}] + [{descriptions[j]}]")})
    rows.sort(key=lambda r: (-r["omega_k_cm1"], abs(r["delta_cm1"])))
    return {"status": "harmonic_screen", "delta_cm1": delta, "min_subsystem_amplitude_pct": min_amp,
            "min_ir_km_mol": min_ir if ir is not None else None,
            "n_modes_screened": int(sel.size), "n_fundamentals_screened": int(fund.size), "n_candidates": len(rows),
            "candidates": rows, "overtones": overtones,
            "caveat": "Harmonic near-coincidences only (candidates). On four GVPT2 test cases this screen contained "
                      "60-78% of the Fermi resonances Gaussian treated, and 7-9% of its candidates were among "
                      "them. Coupling strengths and anharmonic shifts need a GVPT2 calculation."}
