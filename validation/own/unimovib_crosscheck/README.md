# Independent UniMoVib cross-check

`crosscheck_results_v01.json` preserves the original comparison on the owned 614-atom benzaldehyde/200-water
ONIOM benchmark. Its input SHA-256 is the same as
`../oniom_benzaldehyde_water/benzaldehyde_water_oniom.log.sha256`. The record is historical, not a new beta run.

`crosscheck_results_beta_v01.json` records a fresh run of the portable recipe with PEDRO 0.2.0b1 and the same input
hash. All 1,836 full-model frequencies and IR intensities agree with UniMoVib to its print precision (5e-5 cm^-1
and 5e-5 km/mol); the rev-GSVA Hessians differ by at most 6.14e-11 Eh/bohr².

The portable recipe is `scripts/crosscheck_unimovib.py`, adapted from the original local script without changing
the comparison formula. It uses Gaussian's printed five-decimal masses, projects global rigid motion for the
full-model comparison, and uses an absolute Hessian pseudo-inverse cutoff of 1e-12 Eh/bohr² for rev-GSVA. If
UniMoVib prints six additional rigid-motion modes, the six with smallest absolute frequency are omitted from
the full-model comparison; the remaining frequencies are compared in ascending order. No smoothing or scaling
is applied. IR intensities are in km/mol and Hessians in Eh/bohr².

Install UniMoVib separately from [its source repository](https://github.com/zorkzou/UniMoVib); the recorded run
used version 1.5.4. No UniMoVib executable, source or third-party example is bundled here. The 42 MB Gaussian
benchmark log is also external; see its [availability status](../oniom_benzaldehyde_water/README.md).

From the PEDRO repository root, with the log and executable at your own paths:

```powershell
python scripts/crosscheck_unimovib.py make --log C:/data/benzaldehyde_water_oniom.log --workdir C:/work/pedro_crosscheck_v01
python scripts/crosscheck_unimovib.py run --log C:/data/benzaldehyde_water_oniom.log --workdir C:/work/pedro_crosscheck_v01 --executable C:/tools/UniMoVib/unimovib.exe
python scripts/crosscheck_unimovib.py compare --log C:/data/benzaldehyde_water_oniom.log --workdir C:/work/pedro_crosscheck_v01
```

`make` requires a new work directory; `run` refuses to overwrite program output, and `compare` refuses to replace
its JSON record. The subsystem defaults to benchmark atoms 1–14 (1-based); supply `--subsystem` explicitly for
another nonlinear molecule. Full/GSVA inputs and outputs stay in that work directory. The rev-GSVA calculation
is a validation cross-check implemented in the script, not a new PEDRO analysis case or an internal-coordinate
decomposition of the full-model modes.
