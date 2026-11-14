# Example calculations

Small public calculations for trying PEDRO and for the README screenshots (`scripts/make_screenshots.py`).

| File | Program | What it shows |
|---|---|---|
| `acetone_2h2o_b3lyp_orca.hess` (+ `.inp`) | ORCA 6.1.1, B3LYP-D3(BJ)/def2-SVP, Opt Freq | A molecule with an environment: PEDRO describes acetone in its own internal coordinates and gives the two hydrogen-bonded waters' share (case F). The optimized structure keeps one imaginary mode of −10.8 cm⁻¹ (a floppy water motion), which PEDRO reports. |
| `benzaldehyde_gvpt2_b3lyp_gaussian.log` (+ `.gjf`) | Gaussian 16, B3LYP/6-31+G(d), Opt=Tight Freq=Anharmonic | 702 GVPT2 states with 12 Fermi resonances (including the aldehyde C–H stretch with the overtone of the C–H in-plane bend) and Darling–Dennison couplings; used for the README screenshot and the social preview. User and host fields in the archive block replaced by placeholders. |
| `ethylene_gvpt2_b3lyp_gaussian.log` (+ `.gjf`) | Gaussian 16, B3LYP/6-31+G(d), Opt=Tight Freq=Anharmonic | GVPT2 states, a Fermi resonance and Darling–Dennison couplings, quality verdicts and the D2h symmetry gauge. The user and host fields in the log's archive block were replaced by placeholders; nothing else was changed. |

```bash
pedro analyze examples/acetone_2h2o_b3lyp_orca.hess --scale 0.965
pedro analyze examples/benzaldehyde_gvpt2_b3lyp_gaussian.log
```
