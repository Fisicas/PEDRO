# Input formats

PEDRO needs, for one structure: the elements, the Cartesian coordinates, the masses and the Cartesian Hessian.
Dipole derivatives (for IR intensities), the gradient (for the stationarity check), the program's printed
frequencies (for the reconstruction check) and ONIOM layers (for the automatic choice of subsystem) are used when
present. Any quantum-chemistry program that writes a Cartesian Hessian can feed PEDRO; the readers below only differ
in where they look for these quantities.

| File | Program | Where PEDRO looks |
|---|---|---|
| `.log`, `.out` | Gaussian 09/16 | The archive block at the end of a `Freq` job (Hessian, geometry, dipole derivatives, forces); printed masses, frequencies, IR intensities and ONIOM layers from the body of the log; a `Freq=Anharmonic` GVPT2 section when present. No checkpoint needed. |
| `.fchk`, `.fch` | Gaussian (`formchk`) | `Cartesian Force Constants`, `Current cartesian coordinates`, `Real atomic weights`, `Dipole Derivatives`. A log with the same name next to it adds printed frequencies and ONIOM layers. |
| `.hess` | ORCA | `$hessian`, `$atoms` (element, mass, coordinates in bohr), `$dipole_derivatives`; `$vibrational_frequencies` and `$ir_spectrum` for the check. Masses are used as ORCA wrote them. |
| `.fmu` | VEDA 4 input | Same fields as `.fchk`. A `.dd2` with the same name is used as the coordinate set, and a `.ved` is compared with PEDRO entry by entry. |
| `.npz`, `.json` | anything | The program-independent format below. |

## Program-independent format (`pedro-hessian/1`)

For Psi4, PySCF, Q-Chem, CP2K, NWChem, xtb or any other program, write these arrays (atomic units except masses) to
an `.npz` file, or to a `.json` file with `"format": "pedro-hessian/1"`:

| Key | Shape | Units | Required |
|---|---|---|---|
| `elements` (or `Z`) | N | symbols (or atomic numbers) | yes |
| `xyz` | N × 3 | bohr | yes |
| `hessian` | 3N × 3N | Eh/bohr² | yes |
| `masses` | N | u | no (most-abundant isotopes otherwise) |
| `dipole_derivatives` | 3N × 3 | e (row = Cartesian coordinate, column = μx, μy, μz) | no |
| `gradient` | 3N | Eh/bohr | no |
| `layers` | N | `H`/`M`/`L` | no |

Example: whatever program produced them, once the arrays are in hand,

```python
import numpy as np

np.savez("my_molecule.npz",
         elements=["O", "H", "H"],
         xyz=xyz_bohr,                  # (N, 3), bohr
         hessian=hessian,               # (3N, 3N), Eh/bohr^2, rows/columns ordered x1, y1, z1, x2, ...
         masses=[15.999, 1.008, 1.008]) # optional
```

Then `pedro analyze my_molecule.npz`. Programs that return the Hessian as an (N, N, 3, 3) array (atom, atom,
Cartesian, Cartesian), such as PySCF, need `h.transpose(0, 2, 1, 3).reshape(3 * N, 3 * N)` first.

## Adding a reader

A reader is one function that returns a `pedro.job.Job` (see `src/pedro/io/orca.py`, about 100 lines). Register its
file extension in `src/pedro/io/__init__.py` and add a small real output file and a test that checks PEDRO's
recomputed frequencies against the program's printed ones (`tests/test_orca.py`).
