# Harmonic assignments for an embedded molecule

The 0.2.0b1 beta focuses on assigning modes from a completed ONIOM/QM/MM frequency calculation. Gaussian must
have written a Cartesian Hessian; the log's archive block or a formatted checkpoint supplies it. The full-model
Hessian retains coupling between the target molecule and its environment.

## First analysis

```powershell
pedro analyze model_freq.log --case F --anharmonic off
```

Confirm the subsystem atom numbers printed in the summary. The automatic choice is the largest bonded fragment
of the ONIOM high layer that is not tagged as peptide or water. Bond inference and fragment tags are heuristics;
use `--subsystem` to specify your intended target. CLI atom lists are **1-based Gaussian numbering**; Python
`Options(subsystem=[...])` uses **0-based indices**.

For example, if the target is atoms 1–14 in your own model:

```powershell
pedro analyze model_freq.log --subsystem 1-14 --case F --anharmonic off --scale 0.964
```

The scale is an explicit example, not a recommended universal correction. It affects displayed frequencies only.
The Hessian, eigenvectors and TED calculation use unscaled frequencies. No baseline correction, smoothing,
normalization or resampling is part of TED analysis.

## Reuse a coordinate set

To compare with a historical isolated-molecule assignment, retain its `.dd2` and a structure file with the same
atom numbering as that set. PEDRO maps the source molecule onto the selected embedded subsystem by graph
isomorphism, choosing the lowest all-atom RMSD among equivalent maps:

```powershell
pedro analyze model_freq.log --subsystem 1-14 --case F --anharmonic off --coords molecule.dd2 --coords-molecule molecule.fmu
```

Inspect the mapping, RMSD, coordinate-value differences and saved coordinate files. A matching graph is not
independent verification of an atom label or an isotope site. Supply `--atom-map` when a particular mapping is
required. Automatic coordinates can be used when no historical set is available; individual TED entries depend
on that choice, so use the same mapped set for comparisons whenever possible.

## Choose and interpret the case

| Case | Modes | Meaning of non-internal coordinates |
|---|---|---|
| F (default embedded workflow) | Full-model modes with global rigid motion projected out | Environment Cartesians also carry relative rigid motion of the subsystem; environment TED/amplitude equality is approximate |
| Fprime | All unprojected full-model modes | Six subsystem rigid-body coordinates plus environment Cartesians; environment equality is exact |
| P | Modes of the selected mobile Hessian block, without rigid-body projection | Six subsystem rigid-body coordinates plus mobile-environment Cartesians; frozen atoms do not move; environment equality is exact |

```powershell
pedro analyze model_freq.log --subsystem 1-14 --case Fprime --anharmonic off
pedro analyze model_freq.log --subsystem 1-14 --case P --mobile layer:H --anharmonic off
```

P defaults to only the subsystem mobile; `--mobile layer:H` includes the other high-layer atoms as well. These
cases answer different questions and need not yield identical frequencies or distributions. Unprojected
near-zero/negative rigid-motion eigenvalues in Fprime require interpretation before being called physical
instabilities. See [theory section 7](theory.md#7-a-molecule-inside-a-larger-hessian-basispy).

For each run, review frequency reconstruction (when applicable), masses, residual forces when available,
imaginary modes, coordinate conditioning, TED closure, and the case-specific environment check. TED is a signed
coordinate partition, not a probability; negative entries and entries above 100% are retained. Closure does not
validate a structure or a band assignment.

## Outputs and comparison

- `modes.csv` retains all modes, unscaled/scaled frequencies, IR intensities, subsystem TED and subsystem amplitude.
- `ted_subsystem.csv` contains the target's internal-coordinate TED. In Fprime/P it separately exports
  `subsystem_rigid_body_total` and `environment_total`; only the latter is environment TED. With P and only the
  subsystem mobile, environment TED is zero even when subsystem rigid-body TED is substantial.
- `group_totals.csv` separates chemical classes, subsystem rigid motion where present, and environment groups.
  Use `--groups` for residue groups and `--coord-classes` for named internal-coordinate classes.
- `manifest.json` records input hash, options, mass conventions, coordinate mapping, checks and package versions.
- `summary.txt` and `report.html` display the F approximation or Fprime/P identity deviation. The HTML table shows
  modes with at least 20% subsystem mass-weighted amplitude and states this selection; use CSVs for all modes.

Mode arrays and CSVs are ascending in unscaled frequency, with imaginary values negative. Text/HTML mode tables
and spectrum axes use descending wavenumber. Near-degenerate partners may rotate between calculations; inspect
cluster sums and mass-weighted mode overlaps instead of relying on a mode number alone. Every run writes a new
output folder and preserves earlier results.

The [614-atom benchmark](../validation/own/oniom_benzaldehyde_water/README.md) tests scale and numerical handling.
It is a weakly coupled, mechanically embedded water test model. A protein/cofactor application needs its own
subsystem, atom-map and assignment review; the benchmark is not a measured prediction of assignment changes.

## Optional anharmonic tools

Harmonic ONIOM assignments do not require GVPT2. Use `--anharmonic off` for the workflow above. When a compatible
Gaussian GVPT2 log is supplied, PEDRO can read Gaussian's anharmonic positions/intensities and label their harmonic
parents. The harmonic resonance screen gives candidates only; it does not supply anharmonic shifts or intensities.
