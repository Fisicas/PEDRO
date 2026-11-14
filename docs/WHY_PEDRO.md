# Why PEDRO, and not the existing tools?

This page describes the scope of the 0.2.0b1 beta, including what PEDRO does **not** add.

## What you can already do without PEDRO

- **PED/TED of a free molecule:** VEDA 4 does this well, with automatic and optimizable coordinate sets. PEDRO
  reproduces VEDA's TED on identical input ([validation](validation.md)) but adds nothing new for a free molecule, and
  it does not yet optimize coordinate sets as VEDA's MIX/Optimize does.
- **Frequencies, partial Hessians, isotope substitution and atom-group energy shares from a Hessian:** a few dozen
  lines of NumPy. PEDRO packages them with checks, but they are not hard to write.
- **Picking out "the C=O stretch" band:** a crude diagonal measure, *k*Δ*r*²/**x**ᵀ**Hx**, often selects the same mode
  as a full TED. If that is all you need, you may not need a PED at all.

## What needs PEDRO, or something equivalent

**Mode character of a molecule inside a large QM/MM, ONIOM or partial Hessian, in that molecule's own internal
coordinates, with the environment's share in the same closed budget.** For each mode PEDRO gives how much is, for
example, C–O stretch, ring deformation, methyl bend, and environment, summing to 100%. In projected case F the
environment coordinates also carry relative rigid motion of the subsystem, so their TED total is only approximately
equal to the environment's mass-weighted amplitude share; PEDRO displays the deviation. Fprime and P separate the
subsystem rigid-body contribution and satisfy the exact environment identity ([theory §7](theory.md)).

This matters because an atom-group share is not a coordinate share: the fraction of a mode's energy on the atoms of a
C=O group can be several times larger than the C=O stretch contribution, and only the internal-coordinate breakdown
says what the stretch is coupled to (neighboring bends, side chains, the surrounding protein or solvent).

### What VEDA does with a large embedded model (tested)

VEDA 4.99e (`veda4xx.exe`, the largest build) was run by hand on our 614-atom ONIOM model (benzaldehyde in 200 UFF
waters, `validation/own/`) and on blocks cut from the same Hessian (13 November 2026):

| Input to VEDA | What VEDA did |
|---|---|
| Full 614-atom Hessian | Accepted the file, wrote only the input file name to its log, then stopped: no coordinates read, no error message, every button except *Input* stayed disabled. |
| 481 atoms (benzaldehyde plus nearest atoms) | Same as 614 atoms. |
| 480 atoms (benzaldehyde plus nearest atoms) | Read all 480 atoms and the Hessian. *Create DD2* then reported `CAUTION: shortage of coordinates 451 instead of 1434` and `error in DD2, coordinate No. 452` (the cut model is 187 separate fragments, which VEDA's automatic coordinates do not connect); *Calculate* stayed disabled. |
| Benzaldehyde's own 14-atom block of the full Hessian, environment deleted | Calculated, as a free molecule (3N − 6 = 36 modes). |

So the documented 480-atom limit is real and is enforced silently: an input with more atoms fails without telling
the user why. (VEDA behaved the same way, with messages of the same form, on a 682-atom protein QM/MM model.)

The last row is the usual workaround. Deleting the environment's rows and columns leaves a Hessian that is not
invariant to rigid translation and rotation, and VEDA projects six modes out of it as if it were a free molecule; how
those six are removed is a convention, and PEDRO's differs from VEDA's. On the identical 14-atom block with the
identical coordinate set the two programs therefore differ by up to 8.5 cm⁻¹ and 4 percentage points of TED in the
soft modes, while agreeing on the C=O stretch (87% in both). On a Hessian that is invariant to rigid motion, such as
the rev-GSVA effective Hessian of the same benzaldehyde, the two agree in every entry (TED within 0.5 pp, frequencies
within 0.03 cm⁻¹). In this test model the waters (UFF, mechanical embedding) barely couple to benzaldehyde: the full
analysis (case F), the partial Hessian (case P) and the deleted-environment block all give 87% C=O stretch. Whether
deleting the environment changes an assignment therefore depends on how strongly the molecule couples to its
surroundings, which is exactly what the full analysis measures and the workaround cannot.

A published QM/MM study reached the same limit from the other side: its authors state that VEDA could not be applied
to their QM/MM calculations and assigned the bands from separate implicit-solvent calculations (Le, Hastings, Gozem,
*ACS Omega* 2026).

## What PEDRO adds besides that

- **One dropped-in output file, no conversion**: Gaussian `.log`/`.fchk`, ORCA `.hess`, VEDA `.fmu`, or the
  program-independent format ([input formats](input_formats.md)), with checks VEDA does not report: recomputed versus
  printed frequencies, imaginary modes, residual forces, conditioning of the coordinate set.
- **VEDA's coordinate set transferred** from the isolated molecule onto the same molecule in the large model, by graph
  isomorphism.
- **GVPT2 quality control:** per-state flags and a point-group check of Gaussian's resonance couplings, which catch
  artifacts Gaussian does not flag (in one C2v test molecule, a symmetry-forbidden coupling produced a spurious 57%
  mixing).
- **Open, tested, citable code:** results can be regenerated exactly, and every number carries its convention.

## Closest existing work

Read in full where marked; others from abstracts.

| Tool / method | Large QM/MM models? | Relation to PEDRO |
|---|---|---|
| VEDA 4 (Jamróz 2013; read) | No (480 atoms; tested above) | Same TED definition for free molecules; PEDRO reproduces it. |
| GSVA / rev-GSVA in UniMoVib (Kraka group, *JCTC* 2018; *Theor. Chem. Acc.* 2021, read) | Yes | Effective subsystem Hessian with the environment relaxed: "intrinsic" vibrations of the embedded molecule. Removes mixing with the environment by construction; no TED of the computed modes. A complementary cross-check. |
| LModeA / LModeAGen (*J. Phys. Chem. A* 2022; read) | Yes, via GSVA | Local-mode decomposition (CNM) and local force constants; a different quantity from PED/TED. |
| Local normal-mode analysis for QM/MM (Huix-Rotllant, Schwinn, Ferré, *JCTC* 2020, *PCCP* 2021) | Yes, own QM/MM code | Localized modes on fragments; not a PED/TED. |
| Mode / intensity tracking (Reiher; CP2K) | Yes | Selected modes without a full Hessian. |
| PHVA, MBH, VSA (TAMkin; Vester & Olsen, *JCTC* 2024) | Yes | Frequencies with part of the system frozen; no internal-coordinate decomposition. PEDRO's partial-Hessian mode follows Vester & Olsen's advice not to project rigid-body motion. |
| NOMODECO (*J. Chem. Phys.* 2024), vibAnalysis | Small molecules | Free-molecule decompositions. |

**Cross-check with existing programs.** On our 614-atom ONIOM model (benzaldehyde in 200 UFF waters,
`validation/own/`), UniMoVib 1.5.4 (the Kraka group's MIT-licensed implementation of GSVA) reproduces PEDRO's harmonic
frequencies and IR intensities for all 1836 modes to its print precision (5×10⁻⁵), and PEDRO's implementation of the
rev-GSVA effective Hessian of benzaldehyde matches UniMoVib's to 6×10⁻¹¹ Eh/bohr² with the same pseudo-inverse cut-off.
VEDA 4.99e, run on our own Hessians with the same coordinate sets, agrees with PEDRO on every TED entry within 0.5
percentage point (its rounding limit) and every PED entry within 1 pp ([validation](validation.md)). ORCA 6.1.1's own frequencies and IR intensities are reproduced from its `.hess`
file (test suite).

The [portable UniMoVib recipe and recorded comparison](../validation/own/unimovib_crosscheck/README.md) specify the
inputs, mass convention and pseudo-inverse threshold. The owned large ONIOM log is available in the
[benchmark dataset archive](https://doi.org/10.5281/zenodo.24213306); see the
[benchmark description](../validation/own/oniom_benzaldehyde_water/README.md) for provenance and limitations.

**Not claimed:** that no other tool computes an embedded-subsystem PED/TED. We have found none, but the literature
check is not finished.

## The recommended workflow for a molecule in a large model

1. Full-Hessian TED (`--case F`) with the molecule's VEDA coordinate set: what the computed bands are made of,
   including mixing with the environment.
2. Partial Hessians (`--case P`): robustness against artificial coupling to force-field modes.
3. An intrinsic-vibration cross-check with an existing tool (rev-GSVA in UniMoVib), whose subsystem Hessian can also
   be analyzed in VEDA itself.

All three are cleanest at a true minimum; PEDRO reports the residual forces and imaginary modes so this can be judged.
