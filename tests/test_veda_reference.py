"""Reproduction of VEDA 4 TED/PED tables on identical input.

By default the cases in validation/veda_reference: VEDA 4.99e run on our own Hessians (methanol with Jamroz's
published coordinate sets, ethylene and water with VEDA's automatic sets). PEDRO_VEDA_DATA points to another folder
of the same layout: one sub-folder per case with <stem>.fmu, <stem>.dd2 and <stem>.ved.
"""

import pytest

import pedro


def test_veda_reproduction(veda_dir):
    cases = sorted(p for p in veda_dir.iterdir() if p.is_dir() and list(p.glob("*.fmu")))
    assert len(cases) >= 4
    assert cases
    for case in cases:
        res = pedro.analyze(next(case.glob("*.fmu")))
        v = res.validation
        assert v and v["status"] == "compared", case
        assert v["ted_fraction_within_1pp"] >= 0.99, (case.name, v["ted_fraction_within_1pp"])
        assert v["ped_fraction_within_1pp"] >= 0.99, (case.name, v["ped_fraction_within_1pp"])
        assert res.checks["ted_closure_max_pp"] < 1e-9
