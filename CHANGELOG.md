# Changelog

## 0.2.0b1 (14 November 2026): harmonic ONIOM public beta

- First public-beta scope: harmonic assignments for molecules embedded in large ONIOM/QM/MM models. F remains
  the default full-model convention; existing GVPT2 tools remain optional. Publication follows a convincing
  application demonstration and public development.
- Corrected the Fprime/P `ted_subsystem.csv` partition: `subsystem_rigid_body_total` is separate from
  `environment_total`. Per-mode descriptions now distinguish the two. Harmonic frequencies, eigenvectors,
  coordinate sets and TED values are unchanged.
- Case-specific environment TED/amplitude checks are displayed in text, HTML and GUI summaries. Documentation
  qualifies the approximate identity in F and the exact identity in Fprime/P.
- Added the harmonic ONIOM workflow, contribution/support guidance, portable UniMoVib cross-check recipe,
  public-beta release checklist and GUI/package verification.
- Restricted the public log exceptions to reviewed fixtures; the large ONIOM benchmark remains outside Git.
- Windows GitHub installer creates private version environments, updates from newer beta/stable tags on launch,
  retains the previous active version on update failure, and supports version pinning/offline launch. Manifests
  include pip's Git-install provenance when available.
- CSV consumers of Fprime/P must account for the new rigid-body column and the corrected environment-only total.
  The F and isolated CSV layouts are unchanged. The version DOI will be added after Zenodo archives the release.

## 0.2.0.dev4 (12 November 2026): more input programs, release metadata

- ORCA `.hess` reader (`io/orca.py`): Hessian, atoms and masses as ORCA used them, dipole derivatives, printed
  frequencies and IR intensities for the reconstruction check. Tested on a native ORCA 6.1.1 file (`tests/data`):
  frequencies within 0.0003 cm-1 and intensities within 1e-5 km/mol of ORCA's.
- `docs/input_formats.md`: what PEDRO reads from each program, the program-independent `.npz`/`.json` format, and how
  to add a reader.
- VEDA's size limit measured (docs/WHY_PEDRO.md): VEDA 4.99e stops without a message above 480 atoms.
- Author, license holder and citation metadata (Michael Raymond Nelson); repository URL.
- VEDA interchange: `.fmv` reader (VEDA's own Hessian file; unit factor measured against Gaussian's printed Hessian),
  `pedro.io.veda.write_fmu`, and `.dd2` export in VEDA's exact column layout (checked by loading it in VEDA 4.99e).
- Validation rebuilt on public data computed for PEDRO: VEDA 4.99e runs on our own Hessians
  (`validation/veda_reference/`, now tested by default), VEDA's own example files (statistics only), the methanol PED
  of Jamróz 2013 recomputed at the published level, a 614-atom ONIOM model, and GVPT2 examples.
- Fixed: the GVPT2 symmetry gauge labeled re-diagonalized variational states by their largest component, so two
  strongly mixed states (Darling–Dennison 1-1, benzaldehyde) could receive the same label and the parse check failed
  (78 cm⁻¹); labels are now assigned one-to-one. Results for previously tested molecules are unchanged.
- `docs/references.md`: every paper and program used, checked against Crossref.
- Public example calculations (`examples/`: ORCA acetone with two waters, Gaussian GVPT2 ethylene) and
  `scripts/make_screenshots.py`, which regenerates the README images from them.
- GVPT2 reference tests use neutral fixture names (`c2v_neutral_gvpt2.out`, `c2v_anion_gvpt2.log`,
  `water_dimer_gvpt2.out`); `scripts/make_branding.py` takes the preview caption and output folder as arguments.

## 0.2.0.dev3 (12 November 2026): reproducing published analyses

- Coordinate-set transfer by graph isomorphism (`mapping.py`, no NetworkX): a `.dd2` made for an isolated molecule is
  placed on the same molecule inside a QM/MM model automatically, choosing the isomorphism with the lowest RMSD.
- `masses='printed'` (Gaussian's 5-decimal masses exactly) and per-atom mass overrides for isotopologues.
- User-defined coordinate classes for group totals (`stretch` / `any` / `rest` rules).
- `metrics.py`: strongest mode (largest TED x IR), TED- and TED x IR-weighted centroids, class totals, mode following
  and one-to-one matching by mass-weighted overlap.
- Checked against an independent, unpublished analysis of an embedded molecule: identical atom maps, per-mode TED
  within 2e-8 pp, and identical band positions and isotope shifts at the stored precision.

## 0.2.0.dev2 (12 November 2026): resonances and overtones

- GVPT2 layer ported from dev18-dev21 and made to work from the Gaussian `.log` alone (no checkpoint): every
  fundamental, overtone and 1+1 combination with Gaussian's positions and intensities, Fermi and Darling-Dennison
  partners, variational projections, the dev18-dev21 quality audit (unchanged thresholds) and symmetry gauge, a
  per-state verdict and a TED-based harmonic-parentage label. Gaussian's mode numbering is mapped onto PEDRO's by
  matching harmonic energies (Hungarian assignment); all flags are reported in PEDRO's numbering.
- Large-model isolation rule: the layer runs only when GVPT2 data exist, fails closed with a reason (imaginary
  modes, foreign Hessian, reduced-dimensional runs), and never touches the harmonic results.
- Symmetry-element search rewritten: candidates only from exchangeable atom pairs, hashed de-duplication and
  probe-atom rejection. A 682-atom model is classified in about 0.1 s (dev21 effectively hung).
- Harmonic overtone and Fermi-candidate screen for jobs without GVPT2, with its recall measured against GVPT2.
- Automatic subsystem skips peptide and water fragments; all fragments are listed with formulas and atom numbers.
- Archive reader handles a closing `\@` wrapped onto its own line (some Gaussian 09 logs).
- Branding: logo C (1:2 Fermi-resonance figure) in light, dark and adaptive versions; data-driven social previews.
- GUI: "Resonances & overtones" tab; fragment table in the summary.

## 0.2.0.dev1 (12 November 2026), first release candidate of the cleaned core

New repository, ported from the frozen 0.1.0.dev21 research tree (which is not modified).

- Drop-in reading of Gaussian `.log`/`.out` files from the archive block (Hessian, geometry, dipole derivatives,
  forces), so ONIOM jobs need no checkpoint. Also `.fchk`, VEDA `.fmu`, and a program-independent `.npz`/`.json`.
- Analytic Wilson B-matrix rows for stretches, bends, torsions and out-of-plane angles (finite differences kept
  as a test oracle).
- TED for every mode, including imaginary modes; force-constant conventions only where the eigenvalue is positive.
- VEDA-compatible PED (inferred rule), with automatic comparison when a `.ved` file sits next to the input.
- Embedded-subsystem bases (cases F, Fprime, P) for a molecule inside a large QM/MM Hessian, with automatic
  subsystem detection from ONIOM layers, group totals, and the exact amplitude identity checked on every run.
- Partial-Hessian analysis.
- Near-degenerate mode clusters and cluster-summed TED.
- Desktop GUI with drag and drop; HTML report; every run writes a new folder and never overwrites.
- Not yet ported: the GVPT2 resonance/overtone layer (detected and reported as not available).
