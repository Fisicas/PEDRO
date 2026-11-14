# Own ONIOM test model: benzaldehyde in 200 UFF waters (614 atoms)

`make_model_v01.py 200` writes `benzaldehyde_water_oniom.gjf` (deterministic, seed 20261113): benzaldehyde in the high
layer (B3LYP/6-31G(d)), 200 randomly oriented waters in a sphere of 11.4 Å in the low layer (UFF), mechanical
embedding. Gaussian 16 A.03: `ONIOM(B3LYP/6-31G(d):UFF) Opt Freq` (the file here runs on 2 cores; 15 min).
The optimization converged and every frequency is real (lowest 1.4 cm⁻¹).

The 42 MB output log is not stored in Git; its SHA-256 is in `benzaldehyde_water_oniom.log.sha256`.
**Availability: the public data archive is pending.** Add the actual persistent URL/DOI here before advertising
downloadable reproduction of this benchmark; no placeholder DOI is a published record. The supplied `.gjf` and
deterministic model generator are public inputs for regenerating the calculation with Gaussian.

With an existing copy of the log, from the repository root:

```powershell
python scripts/run_validation.py --oniom C:/data/benzaldehyde_water_oniom.log --output-dir C:/work/pedro_validation_v01
```

Use a new output directory. This writes versioned JSON and Markdown records there. Input files are read only;
their hashes are recorded. Check the log hash against the supplied `.sha256` before comparing results.

It is a test object for PEDRO's ONIOM handling (layers, automatic choice of the described molecule, embedded
bases F/Fprime/P, timing, memory), not a physical solvation model.

Collaborator-owned research ONIOM logs are not redistributed with PEDRO. This 614-atom benchmark was generated
specifically for the software and is separate from those research datasets.

## Full-model PED/TED extraction check

The [13 November 2026 check record](ped_ted_check_20261113_v01.json) documents a fresh run of PEDRO 0.2.0b1 on the
42,513,569-byte log with the SHA-256 above. Case F retains all 1,836 vibrational modes, with 36 automatic internal
coordinates for benzaldehyde (atoms 1–14) and 1,800 environment Cartesians. The maximum frequency difference from
Gaussian is 0.0301 cm^-1. Both CSV exports were read back and checked, and unrounded arrays were saved.

TED closes within 5.69e-13 percentage points per mode. Row-partition PED derived directly from the Cartesian
quadratic energy agrees with TED within 4.54e-8 percentage points. The separately exported Morino–Kuchitsu
diagonal PED requires cross terms to close the energy budget; it is not a normalized probability distribution.
Soft environment modes have large canceling terms. VEDA-balanced PED is not computed for embedded case F.

For example, the C12–O13 stretch has 87.1075% TED and 94.5051% diagonal PED at 1803.3402 cm^-1 (unscaled,
PEDRO mode 1430). This is a numerical extraction/check example, not validation of a protein/cofactor assignment.
