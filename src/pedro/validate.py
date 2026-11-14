"""Comparison with VEDA 4 output on identical input (Hessian, masses and .dd2 coordinate set)."""

from __future__ import annotations

import numpy as np

from .distributions import clusters, cluster_sum
from .io.veda import read_ved


def compare_with_veda(res, ved_path, cluster_cm1=1.0):
    v = read_ved(ved_path)
    n = res.ted.shape[0]
    if v["ted"].shape != res.ted.shape:
        return {"status": "skipped", "reason": f"VEDA table {v['ted'].shape} vs PEDRO {res.ted.shape}"}
    # VEDA lists modes in descending frequency; PEDRO ascending. Matching here is by frequency order.
    vt, vp, vf = v["ted"][::-1], v["ped"][::-1], v["freqs"][::-1]
    T = res.ted
    d = np.abs(T - vt)
    cid = clusters(res.modes.freqs, cluster_cm1)
    dc = np.abs(cluster_sum(T, cid) - cluster_sum(vt, cid))
    out = {
        "status": "compared", "ved": str(ved_path), "veda_version": v["version"], "n_modes": n,
        "mode_matching": "ascending frequency order (VEDA rows reversed)",
        "freq_max_abs_diff_cm1": float(np.max(np.abs(np.sort(res.modes.freqs) - vf))),
        "ted_fraction_within_1pp": float(np.mean(d < 1.0)),
        "ted_fraction_within_0p5pp": float(np.mean(d <= 0.5 + 1e-9)),
        "ted_max_abs_diff_pp": float(d.max()),
        "ted_cluster_summed_max_abs_diff_pp": float(dc.max()),
        "ted_rms_diff_pp": float(np.sqrt(np.mean(d**2))),
        "note": "VEDA prints integer percentages, so 0.5 pp is the best agreement possible for any entry",
    }
    if res.veda_ped is not None:
        # The overall sign of a normal mode is arbitrary, so PED rows are compared after choosing, per mode, the
        # sign that best matches VEDA. Relative signs between coordinates within a mode are compared as printed.
        s = np.sign(np.sum(res.veda_ped * vp, axis=1, keepdims=True))
        s[s == 0] = 1
        dp = np.abs(s * res.veda_ped - vp)
        out["ped_fraction_within_1pp"] = float(np.mean(dp < 1.0))
        out["ped_fraction_magnitude_rounds_to_veda"] = float(np.mean(np.round(np.abs(res.veda_ped)) == np.abs(vp)))
        out["ped_max_abs_diff_pp"] = float(dp.max())
        out["ped_sign_note"] = "per-mode overall sign aligned to VEDA (a mode's sign is arbitrary)"
    big = np.argwhere(d >= 1.0)
    worst = []
    for k, i in big[np.argsort(-d[big[:, 0], big[:, 1]])][:10]:
        worst.append({"mode": int(k + 1), "freq_cm1": float(res.modes.freqs[k]), "coordinate": int(i + 1),
                      "label": res.basis.labels[i], "pedro": float(T[k, i]), "veda": float(vt[k, i]),
                      "cluster_size": int(np.sum(cid == cid[k]))})
    out["largest_differences"] = worst
    return out
