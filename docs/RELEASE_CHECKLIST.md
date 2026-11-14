# First public beta release record

Published version: **0.2.0b1**, November 14, 2026. The GitHub and Zenodo links below identify the actual records.
No PyPI upload is claimed. Publish only this `pedro/` tree, not the surrounding workspace.

## Decisions

| Decision | State |
|---|---|
| Main scope | Harmonic mode assignments for molecules embedded in large ONIOM/QM/MM models |
| Release status | Public beta, 0.2.0b1; F remains the default full-model analysis |
| Publication | Release first, develop publicly, choose the paper route after a convincing application demonstration |
| Public data | Calculations made for PEDRO and transcribed published values only; collaborator-owned research ONIOM logs are excluded, and VEDA's distributed examples and private research references stay external |
| Large benchmark | Owned 42 MB ONIOM log and derived PED/TED outputs: [dataset DOI](https://doi.org/10.5281/zenodo.24213306), identified by SHA-256 |
| Repository / Python distribution | `github.com/Fisicas/PEDRO` / `pedro-ped` |
| License and author | MIT, Michael Raymond Nelson |
| Logo | Option C, with light/dark/adaptive variants in `branding/` |
| Delivery | GitHub install command, private Windows environment, automatic synchronization to newer beta/stable tags on launch; manual pip/wheel route also available |

## Local preparation

- [x] Author, license holder, beta version and citation metadata; set the actual release date at publication.
- [x] Harmonic ONIOM workflow and explicit F/Fprime/P conventions.
- [x] Separate subsystem rigid-body and environment TED in Fprime/P CSVs and descriptions.
- [x] Display the approximate F or exact Fprime/P environment identity check in text, HTML and GUI.
- [x] Public Gaussian benzaldehyde/ethylene GVPT2 fixtures and reference tests; optional private fixtures remain optional.
- [x] ORCA reader, program-independent format and owned VEDA reference cases.
- [x] Reference metadata and Crossref lookup records in `docs/references_check/`; this is bibliographic checking,
      not a completed full-text novelty review.
- [x] README screenshots use owned public example calculations; reproduce with `scripts/make_screenshots.py`.
- [x] Portable UniMoVib cross-check recipe and recorded owned-model comparison.
- [x] Contribution and support instructions.
- [x] GitHub managed installer/launcher, explicit update command, offline launch and version pinning.
- [x] Tests, clean built-wheel installation, managed installer regression and GUI checks; see [beta validation](beta_validation.md).
- [x] Inspect the candidate file list and local changes; exclude raw private calculations, vendor programs and build artifacts.
- [x] Run the owned 42 MB ONIOM benchmark through full-model PED/TED export and Cartesian-energy checks;
      [record](../validation/own/oniom_benzaldehyde_water/ped_ted_check_20261113_v01.json). The water model establishes
      numerical handling and scale, with weak coupling; it is not evidence for a changed protein/cofactor assignment.

## Publication actions

- [x] Create/verify the GitHub repository and remote, then push the reviewed code and run CI.
      [Release CI run 43561441834](https://github.com/Fisicas/PEDRO/actions/runs/43561441834) passed for commit `b6344c4`.
- [x] Publish the [ONIOM benchmark archive](https://doi.org/10.5281/zenodo.24213306) and link the benchmark README and validation.
- [x] Set the actual release date in `CITATION.cff`, date the changelog, and create tag `v0.2.0b1`.
- [x] Publish the [GitHub beta release](https://github.com/Fisicas/PEDRO/releases/tag/v0.2.0b1).
- [x] Archive the [software release on Zenodo](https://doi.org/10.5281/zenodo.24213384) and add its actual DOI to the citation and README.

Follow the [archiving procedure](archiving.md) for the software and the separate owned benchmark dataset.

## Later work

- Set the GitHub social preview and visually check adaptive branding in dark mode if not already configured.
- Finish the full-text comparison with mode tracking, NOMODECO, GAR2PED and MOLVIB before claiming novelty.
- Demonstrate scientifically useful assignments in a realistically coupled application and select a paper venue.
- Reduced-dimensional GVPT2, hot-band/ladder positions, linear-bend coordinates and a Windows executable are future
  features, not prerequisites for this harmonic beta. Current limitations remain documented.
