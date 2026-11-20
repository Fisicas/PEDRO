# How PEDRO computes PED and TED

This note derives every quantity PEDRO reports, starting from the Cartesian Hessian (and, for §10, Gaussian's GVPT2 output) that a quantum-chemistry
program writes. Each section names the module that implements it. Nothing here is new theory: it is the
Wilson GF method and the energy-distribution definitions of Morino–Kuchitsu and Keresztury–Jalsovszky, written
out so that every number in the output can be traced.

**Notation.** $`N`$ atoms, Cartesian displacements $`\mathbf{x}\in\mathbb{R}^{3N}`$ (bohr), atomic masses $`m_a`$ (u),
$`\mathbf{M}=\mathrm{diag}(m_1,m_1,m_1,\dots,m_N)`$. $`\mathbf{H}=\partial^2E/\partial\mathbf{x}\,\partial\mathbf{x}^{\mathsf T}`$
is the Cartesian Hessian ($`E_h\,\mathrm{bohr}^{-2}`$). Internal coordinates are $`s_i`$, $`i=1\dots n`$; normal
coordinates are $`Q_k`$, $`k=1\dots n`$.

---

## 1. Normal modes from the Cartesian Hessian (`harmonic.py`)

**Mass weighting.** With $`\tilde{\mathbf{H}}=\mathbf{M}^{-1/2}\mathbf{H}\mathbf{M}^{-1/2}`$ the harmonic equations
of motion become an ordinary symmetric eigenproblem.

**Removing rigid-body motion.** The six (five for a linear molecule) mass-weighted translation and rotation
vectors are

```math
\mathbf{d}^{(T_\alpha)}_a=\sqrt{m_a}\,\mathbf{e}_\alpha,\qquad
\mathbf{d}^{(R_\alpha)}_a=\sqrt{m_a}\,\mathbf{e}_\alpha\times\mathbf{r}_a ,
```

with $`\mathbf{r}_a`$ the position relative to the center of mass. They are orthonormalized (SVD) into
$`\mathbf{D}\in\mathbb{R}^{3N\times 6}`$, and $`\mathbf{V}`$ is an orthonormal basis of the complement
($`\mathbf{V}^{\mathsf T}\mathbf{D}=0`$). The vibrational problem is

```math
\left(\mathbf{V}^{\mathsf T}\tilde{\mathbf{H}}\mathbf{V}\right)\mathbf{W}=\mathbf{W}\boldsymbol{\Lambda},\qquad
\mathbf{L}_q=\mathbf{V}\mathbf{W},\qquad \mathbf{L}_x=\mathbf{M}^{-1/2}\mathbf{L}_q .
```

$`\mathbf{L}_q`$ has orthonormal columns (mass-weighted mode vectors); column $`k`$ of $`\mathbf{L}_x`$ is the Cartesian
displacement per unit normal coordinate, $`\mathbf{x}=\mathbf{L}_x\mathbf{Q}`$.

**Frequencies.**

```math
\tilde\nu_k=\frac{\mathrm{sign}(\lambda_k)\sqrt{|\lambda_k|}}{2\pi c}
=5140.4871\;\mathrm{sign}(\lambda_k)\sqrt{|\lambda_k|}\ \ \mathrm{cm^{-1}},\quad \lambda_k\ \text{in}\ E_h\,\mathrm{bohr^{-2}\,u^{-1}} .
```

Negative $`\lambda_k`$ (imaginary modes) are reported as negative wavenumbers. They are never replaced by
$`|\lambda_k|`$ and never dropped.

**IR intensities.** With dipole derivatives $`\partial\mu_\alpha/\partial x_i`$ (atomic units),

```math
\frac{\partial\mu_\alpha}{\partial Q_k}=\sum_i\frac{\partial\mu_\alpha}{\partial x_i}(\mathbf{L}_x)_{ik},\qquad
I_k=974.8801\sum_\alpha\left(\frac{\partial\mu_\alpha}{\partial Q_k}\right)^2\ \mathrm{km\,mol^{-1}} .
```

**Check against the source program.** When the log prints frequencies, PEDRO reports
$`\max_k|\tilde\nu_k^{\mathrm{PEDRO}}-\tilde\nu_k^{\mathrm{printed}}|`$. For the validation jobs this is at most
0.0003 cm⁻¹ for the Gaussian and ORCA free-molecule jobs and 0.030 cm⁻¹ for the 614-atom ONIOM model
([validation](validation.md)).

**Away from a stationary point.** The projected Hessian is the full curvature only at a stationary point. When the
gradient is not zero, rotations couple to it, and the curvilinear transformation to internal coordinates needs an
extra term $`\sum_i(\partial E/\partial s_i)\,\partial^2 s_i/\partial\mathbf{x}\,\partial\mathbf{x}^{\mathsf T}`$
that this analysis (like VEDA's) does not include. PEDRO reports the largest atomic force so that the user can
judge this; it never corrects for it silently.

## 2. Internal coordinates and the Wilson B matrix (`coordinates.py`)

A coordinate $`s_i(\mathbf{x})`$ is a linear combination of primitives. Its B-matrix row is the gradient
$`B_{ij}=\partial s_i/\partial x_j`$, evaluated analytically. In what follows $`\mathbf{u}_{ab}=\mathbf{x}_b-\mathbf{x}_a`$,
$`r_{ab}=|\mathbf{u}_{ab}|`$, $`\hat{\mathbf{u}}=\mathbf{u}/r`$.

**Stretch** $`r=r_{ab}`$:
$`\ \partial r/\partial\mathbf{x}_a=-\hat{\mathbf{u}}_{ab},\ \ \partial r/\partial\mathbf{x}_b=\hat{\mathbf{u}}_{ab}`$.

**Bend** $`\theta`$ at apex $`b`$, with $`\hat{\mathbf{u}}=\hat{\mathbf{u}}_{ba}`$, $`\hat{\mathbf{v}}=\hat{\mathbf{u}}_{bc}`$,
$`\cos\theta=\hat{\mathbf{u}}\cdot\hat{\mathbf{v}}`$:

```math
\frac{\partial\theta}{\partial\mathbf{x}_a}=\frac{\cos\theta\,\hat{\mathbf{u}}-\hat{\mathbf{v}}}{r_{ba}\sin\theta},\quad
\frac{\partial\theta}{\partial\mathbf{x}_c}=\frac{\cos\theta\,\hat{\mathbf{v}}-\hat{\mathbf{u}}}{r_{bc}\sin\theta},\quad
\frac{\partial\theta}{\partial\mathbf{x}_b}=-\frac{\partial\theta}{\partial\mathbf{x}_a}-\frac{\partial\theta}{\partial\mathbf{x}_c}.
```

**Torsion** $`\phi(a,b,c,d)\in(-\pi,\pi]`$ (IUPAC sign). With $`\mathbf{F}=\mathbf{x}_a-\mathbf{x}_b`$,
$`\mathbf{G}=\mathbf{x}_b-\mathbf{x}_c`$, $`\mathbf{H}=\mathbf{x}_d-\mathbf{x}_c`$, $`\mathbf{A}=\mathbf{F}\times\mathbf{G}`$,
$`\mathbf{B}=\mathbf{H}\times\mathbf{G}`$ (Blondel and Karplus, singularity-free form):

```math
\begin{aligned}
\partial\phi/\partial\mathbf{x}_a&=-\tfrac{|\mathbf{G}|}{A^2}\mathbf{A}, &
\partial\phi/\partial\mathbf{x}_d&=\tfrac{|\mathbf{G}|}{B^2}\mathbf{B},\\[0pt]
\partial\phi/\partial\mathbf{x}_b&=\tfrac{|\mathbf{G}|}{A^2}\mathbf{A}+\tfrac{\mathbf{F}\cdot\mathbf{G}}{A^2|\mathbf{G}|}\mathbf{A}-\tfrac{\mathbf{H}\cdot\mathbf{G}}{B^2|\mathbf{G}|}\mathbf{B}, &
\partial\phi/\partial\mathbf{x}_c&=\tfrac{\mathbf{H}\cdot\mathbf{G}}{B^2|\mathbf{G}|}\mathbf{B}-\tfrac{\mathbf{F}\cdot\mathbf{G}}{A^2|\mathbf{G}|}\mathbf{A}-\tfrac{|\mathbf{G}|}{B^2}\mathbf{B}.
\end{aligned}
```

**Out-of-plane** (VEDA `OUT A B C D`, $`D`$ central): the angle between the bond $`D\!\to\!A`$ and the plane $`(B,C,D)`$.
With $`\mathbf{u}=\mathbf{x}_A-\mathbf{x}_D`$, $`\mathbf{v}=\mathbf{x}_B-\mathbf{x}_D`$, $`\mathbf{w}=\mathbf{x}_C-\mathbf{x}_D`$,
$`\mathbf{n}=\mathbf{v}\times\mathbf{w}`$ and $`\sin\gamma=\mathbf{n}\cdot\mathbf{u}/(|\mathbf{n}||\mathbf{u}|)`$:

```math
\mathbf{g}\equiv\frac{\partial\sin\gamma}{\partial\mathbf{n}}=\frac{\mathbf{u}}{|\mathbf{n}||\mathbf{u}|}-\frac{\sin\gamma\,\mathbf{n}}{|\mathbf{n}|^2},\qquad
\frac{\partial\gamma}{\partial\mathbf{x}_A}=\frac{1}{\cos\gamma}\Big(\frac{\mathbf{n}}{|\mathbf{n}||\mathbf{u}|}-\frac{\sin\gamma\,\mathbf{u}}{|\mathbf{u}|^2}\Big),
```

```math
\frac{\partial\gamma}{\partial\mathbf{x}_B}=\frac{\mathbf{w}\times\mathbf{g}}{\cos\gamma},\qquad
\frac{\partial\gamma}{\partial\mathbf{x}_C}=\frac{\mathbf{g}\times\mathbf{v}}{\cos\gamma},\qquad
\frac{\partial\gamma}{\partial\mathbf{x}_D}=-\sum_{X=A,B,C}\frac{\partial\gamma}{\partial\mathbf{x}_X}.
```

Every analytic row is tested against central finite differences of the coordinate value (relative error
< 10⁻⁶, typically 10⁻⁹) and for invariance to rigid translation and rotation (`tests/test_core.py`).

**Coordinate sets.** PEDRO uses (i) a VEDA `.dd2` file (read and written, including multi-component coordinates),
(ii) a PEDRO JSON set, or (iii) an automatic set: bond stretches, then XH₂ symmetric/antisymmetric
combinations, bends, out-of-plane angles and torsions, selected greedily by Gram–Schmidt on the B rows until
$`3N-6`$ independent coordinates are found. TED values depend on the coordinate set, so the set used is always
written next to the results (`coordinates.dd2`, `coordinates.json`).

## 3. The mode matrix

For a square, non-redundant basis ($`n`$ coordinates for $`n`$ modes),

```math
\mathbf{L}=\mathbf{B}\,\mathbf{L}_x\qquad(\text{coordinates}\times\text{modes}),\qquad \mathbf{s}=\mathbf{L}\mathbf{Q}.
```

$`L_{ik}`$ is the change of coordinate $`i`$ per unit normal coordinate $`k`$. PEDRO requires $`\mathbf{L}`$ to be
invertible and reports its condition number; a singular $`\mathbf{L}`$ means the coordinates are not independent for
these modes, and the run stops.

## 4. Force constants, PED and TED (`distributions.py`)

**Internal force field.** The harmonic energy is $`V=\tfrac12\sum_k\lambda_kQ_k^2=\tfrac12\mathbf{Q}^{\mathsf T}\boldsymbol{\Lambda}\mathbf{Q}`$.
Substituting $`\mathbf{Q}=\mathbf{L}^{-1}\mathbf{s}`$,

```math
V=\tfrac12\,\mathbf{s}^{\mathsf T}\mathbf{F}\,\mathbf{s},\qquad \mathbf{F}=\mathbf{L}^{-\mathsf T}\boldsymbol{\Lambda}\mathbf{L}^{-1},
\qquad \mathbf{L}^{\mathsf T}\mathbf{F}\mathbf{L}=\boldsymbol{\Lambda}.
```

At a stationary point and with a complete basis this $`\mathbf{F}`$ is the internal-coordinate force field; the test
suite recovers a known valence force field exactly (to 10⁻¹⁰).

**Potential energy of one mode.** For mode $`k`$ alone, $`\mathbf{s}=\mathbf{L}_{:,k}Q_k`$, so

```math
\lambda_k=\sum_{i}F_{ii}L_{ik}^2+\sum_{i\neq j}F_{ij}L_{ik}L_{jk}.
```

*Morino–Kuchitsu PED* keeps the diagonal terms, $`P^{\mathrm{MK}}_{ki}=100\,F_{ii}L_{ik}^2/\lambda_k`$; with the pair
terms $`200\,F_{ij}L_{ik}L_{jk}/\lambda_k`$ ($`i<j`$) the mode closes to exactly 100 %.

*Row-partition PED* assigns each cross term half to each of its two coordinates:

```math
P_{ki}=\frac{100}{\lambda_k}\,L_{ik}\,(\mathbf{F}\mathbf{L})_{ik}.
```

Because $`\mathbf{F}\mathbf{L}=\mathbf{L}^{-\mathsf T}\boldsymbol{\Lambda}`$, $`(\mathbf{F}\mathbf{L})_{ik}=\lambda_k(\mathbf{L}^{-1})_{ki}`$ and

```math
\boxed{\;T_{ki}=100\,L_{ik}\,(\mathbf{L}^{-1})_{ki}\;}
```

— the force constants cancel. This is the Keresztury–Jalsovszky form, and it is what VEDA prints as TED.

**Why "total" energy.** The kinetic energy is $`\tfrac12\dot{\mathbf{Q}}^{\mathsf T}\dot{\mathbf{Q}}=\tfrac12\dot{\mathbf{s}}^{\mathsf T}\mathbf{L}^{-\mathsf T}\mathbf{L}^{-1}\dot{\mathbf{s}}`$.
Partitioning it the same way gives $`K_{ki}=100\,L_{ik}(\mathbf{L}^{-\mathsf T}\mathbf{L}^{-1}\mathbf{L})_{ik}=100\,L_{ik}(\mathbf{L}^{-1})_{ki}`$.
The potential and kinetic distributions are identical, so $`T`$ is also the distribution of the total energy.

**Closure.** $`\sum_iT_{ki}=100\,(\mathbf{L}^{-1}\mathbf{L})_{kk}=100`$ for every mode, and $`\sum_kT_{ki}=100\,(\mathbf{L}\mathbf{L}^{-1})_{ii}=100`$
for every coordinate. Individual entries are signed and may exceed 100; PEDRO never clips or renormalizes them.

**Invariances** (all tested):
1. *Row scaling.* Replacing $`s_i\to c_is_i`$ multiplies row $`i`$ of $`\mathbf{L}`$ by $`c_i`$ and column $`i`$ of
   $`\mathbf{L}^{-1}`$ by $`1/c_i`$; $`T`$ is unchanged. Units and the 1 Å scaling of angular rows therefore do not matter.
2. *Mode sign.* $`Q_k\to-Q_k`$ flips column $`k`$ of $`\mathbf{L}`$ and row $`k`$ of $`\mathbf{L}^{-1}`$; $`T`$ is unchanged.
3. *Group totals.* For a group $`g`$ of coordinates and any invertible change of coordinates $`\mathbf{C}`$ acting only
   inside $`g`$, $`\mathbf{L}_g\to\mathbf{C}\mathbf{L}_g`$ and $`(\mathbf{L}^{-1})_{:,g}\to(\mathbf{L}^{-1})_{:,g}\mathbf{C}^{-1}`$, so
   $`\sum_{i\in g}T_{ki}=100\,[(\mathbf{L}^{-1})_{k,g}\mathbf{L}_{g,k}]`$ is unchanged. A group total (for example "all C–O
   stretches", or "the environment") is a property of the subspace the group spans, not of how it is parameterized.

**Imaginary modes.** $`T`$ is defined for every mode because it uses only $`\mathbf{L}`$. Force-constant conventions
divide by $`\lambda_k`$ and are reported only for $`\lambda_k>0`$ (NaN otherwise).

## 5. VEDA-compatible PED

VEDA 4 prints, besides TED, a "PED" table whose absolute values sum to 100 both along each mode and along each
coordinate. VEDA's source code is not public. PEDRO reproduces that table by Sinkhorn balancing:

```math
\mathbf{X}^{(0)}=|\mathbf{T}|,\qquad \mathbf{X}\leftarrow\mathrm{diag}(\mathbf{X}\mathbf{1})^{-1}\mathbf{X},\ \ \mathbf{X}\leftarrow\mathbf{X}\,\mathrm{diag}(\mathbf{1}^{\mathsf T}\mathbf{X})^{-1}\quad\text{(repeat to }10^{-8}),
\qquad P^{\mathrm{VEDA}}_{ki}=100\,\mathrm{sign}(L_{ik})\,X_{ki}.
```

This rule was **inferred** from VEDA 4.99e output; it reproduces VEDA's printed PED within rounding for the
molecules in `docs/validation.md`. A different VEDA version could differ. The overall sign of a mode is arbitrary,
so comparisons align each mode's sign before comparing.

## 6. Near-degenerate modes

If two or more modes are (nearly) degenerate, any rotation of their eigenvectors is an equally valid solution, and
per-mode TED can change by tens of percentage points between programs. For a cluster $`C`$ and an orthogonal
rotation $`\mathbf{R}`$, $`\mathbf{L}_{:,C}\to\mathbf{L}_{:,C}\mathbf{R}`$ and $`(\mathbf{L}^{-1})_{C,:}\to\mathbf{R}^{\mathsf T}(\mathbf{L}^{-1})_{C,:}`$, so

```math
\sum_{k\in C}T_{ki}=100\,\big[\mathbf{L}_{:,C}(\mathbf{L}^{-1})_{C,:}\big]_{ii}
```

is invariant. PEDRO groups modes closer than 1 cm⁻¹ (adjustable) and reports cluster-summed TED alongside
per-mode TED.

## 7. A molecule inside a larger Hessian (`basis.py`)

Let $`S`$ be the subsystem (the molecule to describe, $`N_S`$ atoms) and $`E`$ the environment. PEDRO builds a square
basis from $`3N_S-6`$ internal coordinates of $`S`$ plus coordinates for everything else:

| Case | Modes | Basis | Size |
|---|---|---|---|
| `isolated` | projected | internals of the whole molecule | $`3N-6`$ |
| `F` (default for large/ONIOM inputs) | projected | internals of $`S`$ + Cartesians of all $`E`$ atoms | $`3N_S-6+3N_E=3N-6`$ |
| `Fprime` | unprojected | internals of $`S`$ + 6 rigid-body coordinates of $`S`$ + Cartesians of $`E`$ | $`3N`$ |
| `P` | partial Hessian (§8) | internals of $`S`$ + 6 rigid-body of $`S`$ + Cartesians of mobile $`E`$ | $`3N_M`$ |

The rigid-body coordinates of $`S`$ are the mass-weighted center-of-mass displacement and the linearized rotation
$`\boldsymbol\theta=\mathbf{I}_S^{-1}\sum_{a\in S}m_a\,\mathbf{r}_a\times\Delta\mathbf{x}_a`$ about the principal axes of $`S`$.
Both are linear in $`\mathbf{x}`$, so their B rows are constant.

**Exact identity for `Fprime` and `P`.** Here $`\mathbf{B}=\mathrm{diag}(\mathbf{B}_S,\mathbf{I}_E)`$ is block diagonal and
square, and the mode matrix $`\mathbf{L}_x`$ is square with $`\mathbf{L}_x^{-1}=\mathbf{L}_q^{\mathsf T}\mathbf{M}^{1/2}`$. Then
$`\mathbf{L}^{-1}=\mathbf{L}_q^{\mathsf T}\mathbf{M}^{1/2}\mathbf{B}^{-1}`$ and, for the environment rows,

```math
\sum_{i\in E}T_{ki}=100\sum_{i\in E}\big(m_i^{-1/2}(\mathbf{L}_q)_{ik}\big)\big(m_i^{1/2}(\mathbf{L}_q)_{ik}\big)
=100\sum_{a\in E}\lVert\mathbf{q}_{a,k}\rVert^2 ,
```

the share of mode $`k`$'s mass-weighted amplitude on the environment atoms. This equals the atom-group potential
energy fraction $`\sum_{i\in g}x_i(\mathbf{H}\mathbf{x})_i/\mathbf{x}^{\mathsf T}\mathbf{H}\mathbf{x}`$ (because
$`\mathbf{H}\mathbf{x}_k=\lambda_k\mathbf{M}\mathbf{x}_k`$ for an eigenvector), so the internal-coordinate TED of $`S`$ refines,
and never contradicts, the simpler atom-group partition. PEDRO checks the identity on every `Fprime`/`P` run
(`identity_max_dev_pp`, typically 10⁻¹² pp).

**Case `F`.** The basis is square but not block diagonal in the same sense, because the projected modes exclude
the global rigid motions and the rigid motion of $`S`$ relative to $`E`$ is carried by the environment Cartesians.
The identity then holds only approximately; PEDRO reports the deviation (`identity_F_max_dev_pp`). Use `Fprime`
or `P` to see the rigid-body share of $`S`$ separately.

Text, HTML and GUI summaries display this case-specific check. In `ted_subsystem.csv`, F exports
`environment_total` with its projected convention. Fprime/P export `subsystem_rigid_body_total` and
`environment_total` separately: the latter contains only the mobile environment Cartesians. When only $`S`$ is
mobile in P, `environment_total` is zero; the subsystem's rigid-body contribution may still be nonzero.

**Per-mode summaries.** `subsystem_ted_pct` $`=\sum_{i\in S\text{-internal}}T_{ki}`$; `subsystem_amplitude_pct`
$`=100\sum_{a\in S}\lVert\mathbf{q}_{a,k}\rVert^2`$. Group totals sum TED over coordinate classes of $`S`$ (ν(C–O), ν(C–C),
δ, τ, γ, …) and over environment layers or user-defined residue groups.

**Choosing the subsystem automatically.** With ONIOM layers: the largest covalently bonded fragment of the
high (QM) layer. Without layers: the whole structure if it is one molecule. The choice is printed and can be
overridden by an atom list.

## 8. Partial-Hessian analysis

Only the atoms of a mobile set $`M`$ move; all others are treated as infinitely heavy. The vibrational problem
is the mass-weighted block $`\tilde{\mathbf{H}}_{MM}`$, diagonalized without rigid-body projection (the frozen
atoms anchor the system), giving $`3N_M`$ modes (Head; Li and Jensen). This tests whether a band's character
depends on how much of the environment is allowed to move.

## 9. Units and constants

CODATA 2018: $`E_h=4.3597447222071\times10^{-18}`$ J, $`a_0=0.529177210903`$ Å, u $`=1.66053906660\times10^{-27}`$ kg,
$`c=2.99792458\times10^{10}`$ cm s⁻¹. Masses are read from the input; Gaussian's printed 5-decimal masses are
replaced by exact most-abundant-isotope masses only when they agree to within 6×10⁻⁶ u, and the manifest records
which masses were used. Frequency scale factors are for display only and never enter the analysis.

## 10. Resonances and overtones (`anharmonic/`)

### 10.1 What GVPT2 provides

Second-order vibrational perturbation theory expands the potential to quartic order in the normal coordinates,

```math
V=\tfrac12\sum_i\omega_iq_i^2+\tfrac16\sum_{ijk}\phi_{ijk}q_iq_jq_k+\tfrac1{24}\sum_{ijkl}\phi_{ijkl}q_iq_jq_kq_l ,
```

and gives term values $`E(\mathbf n)=\sum_i\omega_i(n_i+\tfrac12)+\sum_{i\le j}\chi_{ij}(n_i+\tfrac12)(n_j+\tfrac12)+\chi_0`$.
The transitions out of the ground state are therefore

```math
\nu_i=\omega_i+2\chi_{ii}+\tfrac12\sum_{j\ne i}\chi_{ij},\qquad
\nu_{2i}=2\nu_i+2\chi_{ii},\qquad
\nu_{i+j}=\nu_i+\nu_j+\chi_{ij}.
```

**Resonances.** The cubic contribution to $`\chi`$ contains denominators such as $`\omega_i-\omega_j-\omega_k`$ and
$`\omega_i-2\omega_j`$. When a fundamental nearly coincides with a combination or overtone (a Fermi resonance), or two
states of the same quantum number nearly coincide (Darling–Dennison resonances), these terms blow up. GVPT2 removes
the resonant terms from the perturbative sum (deperturbed VPT2, DVPT2) and treats them instead by diagonalizing a
small variational matrix that couples the near-degenerate states. Gaussian decides which terms to treat this way
with a frequency window and a test comparing the perturbative and variational treatments (the thresholds are
printed in the log, e.g. 200 cm⁻¹ and 1 cm⁻¹ for Fermi resonances). After the variational step each final state is
a mixture; Gaussian prints the projection of each deperturbed state on its final state when it is below 80 %.

PEDRO reads all of this from the log (`Anharmonic Infrared Spectroscopy`, `Resonance Analysis`, `Variational Matrix
Definition`, `Projection of DVPT2 states…`, the X-matrix tables and the cubic and quartic force constants). It never
recomputes or rescales Gaussian's positions or intensities.

### 10.2 Mapping Gaussian's states onto PEDRO's modes

Gaussian numbers anharmonic modes by descending harmonic frequency; PEDRO numbers them ascending. PEDRO matches
the printed harmonic energy $`E_\text{harm}`$ of every fundamental to its own frequencies, recomputed from the same
Hessian, by a one-to-one (Hungarian) assignment. The largest residual is reported (typically 5×10⁻⁴ cm⁻¹). The
layer refuses to run, with a stated reason, if the residual exceeds 0.5 cm⁻¹ (the GVPT2 step belongs to another
Hessian), if any harmonic frequency is imaginary (VPT2 is defined only at a minimum), or if Gaussian treated only a
subset of modes. Modes closer than 0.1 cm⁻¹ cannot be told apart by energy; states built on them carry
`MAPPING_AMBIGUOUS`.

### 10.3 Harmonic-parentage labels

Each state is labeled by the TED of its parent harmonic modes, for example
`ν₄₈ + ν₄₉ = [ν(C5–O6) 59 + ν(C1–O2) 15 + …] + [ν(C1–O2) 56 + …]`. This says which coordinates the parent modes
move. It is **not** an anharmonic energy distribution: after a resonance the final state is a mixture, and its
character is the projection-weighted mixture of the parents' characters, which PEDRO reports through the
resonance partners and projection percentage rather than by inventing new numbers.

### 10.4 Quality verdicts

Every state receives a verdict, clean, warning or failed, from flags computed per mode (thresholds are constants in
`anharmonic/quality.py`; every flag carries its numbers):

| Flag | Meaning |
|---|---|
| `ANOMALOUS_SHIFT` | deperturbed $`\nu/\omega`$ outside 0.93–1.01 for 400 < ω < 2500 cm⁻¹ (or off by > 15 % below 400 cm⁻¹): position unusable |
| `QUARTIC_OUTLIER` | a semi-diagonal quartic constant $`\lvert\phi_{iijj}\rvert`$ > 100 cm⁻¹ between two heavy-atom modes (typical values are a few cm⁻¹); usually a numerical-differentiation artifact |
| `RESONANCE_MIXED` | Gaussian's projection below 80 %: the state is a mixture |
| `CONTAMINATED` | more than max(5 cm⁻¹, 20 %) of the mode's shift, $`\tfrac12\sum_j\chi_{ij}`$ over flagged partners $`j`$, comes from modes VPT2 cannot treat |
| `COUPLED_TO_FLAGGED` | the variational matrix couples the mode to a state containing an anomalous mode |
| `PARENT_FLAGGED` | an overtone or combination whose parent mode is flagged |
| `FORBIDDEN_COUPLING` | Gaussian coupled the mode to a state that the point group forbids (§10.5) |
| `VPT2_BREAKDOWN` | a fundamental shifted by more than 25 %, or an overtone/combination below half its harmonic value |

A position flagged `ANOMALOUS_SHIFT` or `VPT2_BREAKDOWN` makes the verdict `failed`; any other flag gives `warning`.

### 10.5 Symmetry gauge

For molecules with an abelian point group (C₂ᵥ, Cₛ, C₂ₕ, D₂ₕ, …) a force constant whose symmetry product is not
totally symmetric must vanish. Gaussian run without symmetry prints such constants anyway, so their size measures
the numerical noise, and a resonance coupling that symmetry forbids is an artifact. PEDRO finds the point group from
the geometry, assigns irreducible representations to its own normal modes, classifies every printed resonance and
variational coupling (symmetry-forbidden, involves a flagged mode, or allowed and clean), and re-diagonalizes
Gaussian's printed variational matrix keeping only allowed, clean couplings. A full re-diagonalization must
reproduce Gaussian's final fundamentals (parse check, typically < 0.05 cm⁻¹). In one C2v test molecule this removes a 6.9 cm⁻¹
forbidden coupling between the B₂ and A₁ C=O stretches, which Gaussian had turned into a 57 % mixing.

The search for symmetry elements uses only directions that can be symmetry elements: principal axes, atom
directions, and differences and sums of position vectors of atoms that could be exchanged (same element, same
distance from the center of mass). Candidate operations are rejected after testing a few probe atoms, so a
614-atom model is classified in about 0.3 s. Above 200 atoms the gauge is skipped by default (`--symmetry force`
overrides this), because large embedded models are almost always C₁.

### 10.6 Without GVPT2: the harmonic resonance screen

GVPT2 on a 2000-mode QM/MM model is not feasible, so for harmonic jobs PEDRO lists where the subsystem's overtones
($`2\omega_i`$) and combinations ($`\omega_i+\omega_j`$) fall and which of them lie within Δ (default 10 cm⁻¹) of an
IR-active (≥ 1 km/mol) subsystem fundamental. Near-coincidence is necessary for a Fermi resonance but not
sufficient: the coupling $`\phi_{kij}`$ and the anharmonic shifts are unknown without an anharmonic calculation.

How good is this screen? On the public GVPT2 examples ([validation](validation.md) §5) the Δ = 10 cm⁻¹ screen
contained 7 of the 12 Fermi resonances Gaussian treated for benzaldehyde (58 %), with 11 % of its 61 candidates among
them, and missed ethylene's single resonance. A locality ranking (geometric-mean overlap of the three modes' amplitude distributions) was tested and
did not separate real resonances from coincidences, so it is not offered. Every row is labeled a candidate.

**Strong pumping (not yet implemented).** Hot-band and ladder transitions from a populated $`v=1`$ level follow
from the same anharmonic constants: $`\nu(1\to2)=\nu_i+2\chi_{ii}`$ and $`\nu(1_i\to1_i1_j)=\nu_j+\chi_{ij}`$. PEDRO keeps
the parsed X matrix so these positions can be added later.

## 11. Transferring a coordinate set, and band-level quantities (`mapping.py`, `metrics.py`)

### 11.1 Transfer by graph isomorphism

A coordinate set (for example VEDA's `.dd2`) is defined on the atom numbering of one structure (its `.fmu`). To use
it for the same molecule inside a larger model, PEDRO builds the covalent-bond graphs of both molecules (Cordero
radii × 1.25), refines atom colors jointly (element, then neighbor colors, Weisfeiler–Lehman), and enumerates
every element-preserving isomorphism by backtracking. Equivalent hydrogens of methyl and CH₂ groups make these
numerous (5,184 for one 31-atom test molecule with several methyl and CH₂ groups). Among them PEDRO keeps the one whose atoms superpose
with the lowest all-atom RMSD after optimal rotation (Kabsch). It then renumbers the coordinates and reports the
largest change of any coordinate value between the two structures. Different elements, bonding or atom counts stop
the run.

### 11.2 Coordinate classes

Group totals can use named classes instead of PEDRO's chemical defaults. Rules are applied in order and the first
match wins: `stretch a b` (a single stretch between exactly these atoms), `any atoms` (a coordinate touching any of
them), `rest`. By §4 a class total does not depend on how the coordinates inside the class are parameterized.

### 11.3 Band metrics

For a coordinate $`j`$ (for example a C–O stretch), a scaled window $`(\nu_\mathrm{lo},\nu_\mathrm{hi})`$ and modes $`k`$
in it with IR intensity $`I_k`$:

```math
k^\star=\arg\max_k T_{kj}I_k,\qquad
\bar\nu_{T}=\frac{\sum_kT_{kj}\tilde\nu_k}{\sum_kT_{kj}},\qquad
\bar\nu_{TI}=\frac{\sum_kT_{kj}I_k\tilde\nu_k}{\sum_kT_{kj}I_k}.
```

$`k^\star`$ is the "strongest mode" of that coordinate. The centroids average over all modes that carry it, so they do
not jump when near-degenerate modes exchange character.

### 11.4 Following modes between calculations

Isotopologues, partial Hessians and different charge states have different mode sets. PEDRO compares them with the
overlap $`|\langle\mathbf q_k^{(a)}|\mathbf q_l^{(b)}\rangle|`$ of mass-weighted mode vectors restricted to the
subsystem atoms and normalized there, each calculation with its own masses. `follow` gives the best partners of one
mode; `match` gives a one-to-one assignment that maximizes the total overlap (Hungarian algorithm) and flags pairs
below 0.8.

## 12. What PEDRO does not claim

- New PED theory: all definitions above are textbook.
- Source-level equivalence to VEDA: the VEDA-compatible PED rule is inferred from VEDA's output.
- Anharmonic energy distributions: GVPT2 labels are harmonic parentage only (§10.3); positions and intensities are Gaussian's.
- Fermi couplings for harmonic jobs: the screen lists candidates only (§10.6).
- Correctness away from stationary points: see §1; the force on every atom is reported.

## References

The full, publisher-checked list (with DOIs) is in [references.md](references.md).


1. E. B. Wilson, J. C. Decius, P. C. Cross, *Molecular Vibrations*, McGraw-Hill, 1955.
2. Y. Morino, K. Kuchitsu, *J. Chem. Phys.* **20**, 1809 (1952).
3. G. Keresztury, G. Jalsovszky, *J. Mol. Struct.* **10**, 304 (1971).
4. M. H. Jamróz, *Spectrochim. Acta A* **114**, 220 (2013) (VEDA).
5. P. Pulay, G. Fogarasi, F. Pang, J. E. Boggs, *J. Am. Chem. Soc.* **101**, 2550 (1979).
6. G. Fogarasi, X. Zhou, P. W. Taylor, P. Pulay, *J. Am. Chem. Soc.* **114**, 8191 (1992).
7. A. Blondel, M. Karplus, *J. Comput. Chem.* **17**, 1132 (1996).
8. R. Sinkhorn, *Ann. Math. Statist.* **35**, 876 (1964).
9. J. D. Head, *Int. J. Quantum Chem.* **65**, 827 (1997).
10. H. Li, J. H. Jensen, *Theor. Chem. Acc.* **107**, 211 (2002).
11. J. W. Ochterski, *Vibrational Analysis in Gaussian*, Gaussian, Inc. white paper (1999).
12. B. Cordero et al., *Dalton Trans.* 2832 (2008) (covalent radii).
13. E. Tiesinga et al., *Rev. Mod. Phys.* **93**, 025010 (2021) (CODATA 2018).
14. V. Barone, *J. Chem. Phys.* **122**, 014108 (2005) (GVPT2).
15. J. Bloino, V. Barone, *J. Chem. Phys.* **136**, 124108 (2012) (VPT2 transition properties and resonances).
16. J. M. L. Martin, T. J. Lee, P. R. Taylor, J.-P. François, *J. Chem. Phys.* **103**, 2589 (1995) (resonance test).
