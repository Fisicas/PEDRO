"""The in-memory description of one vibrational calculation, independent of the file format it came from."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Job:
    source: str                      # absolute path of the file that was read
    sha256: str
    kind: str                        # gaussian_log | fchk | veda_fmu | npz | json
    elements: list[str]
    masses: np.ndarray               # u, one per atom
    mass_source: str                 # where the masses came from (file, printed log, default isotope table)
    xyz_bohr: np.ndarray             # (N, 3)
    hessian: np.ndarray              # (3N, 3N), Eh / bohr^2, same frame as xyz_bohr
    dipole_derivatives: np.ndarray | None = None   # (3N, 3), e; row = Cartesian coordinate, column = mu_x,y,z
    gradient: np.ndarray | None = None             # (3N,), Eh / bohr (sign as stored by the source)
    charge: int | None = None
    multiplicity: int | None = None
    printed_freqs: np.ndarray | None = None        # cm^-1 as printed by the program (validation only)
    printed_ir: np.ndarray | None = None           # km/mol as printed
    printed_masses: np.ndarray | None = None       # u, exactly as printed by the program (5 decimals in Gaussian)
    layers: list[str] | None = None                # ONIOM layer per atom (H/M/L) when known
    frozen: list[int] = field(default_factory=list)   # 0-based atoms Gaussian kept frozen
    route: str = ""
    program: str = ""
    title: str = ""
    termination: str = "unknown"                   # normal | error | incomplete | unknown
    anharmonic_present: bool = False               # a GVPT2 (Freq=Anharmonic) section exists
    notes: list[str] = field(default_factory=list)
    attachments: dict = field(default_factory=dict)  # sibling files found next to the source (e.g. VEDA .dd2/.ved)

    @property
    def natoms(self) -> int:
        return len(self.elements)

    @property
    def formula(self) -> str:
        counts: dict[str, int] = {}
        for e in self.elements:
            counts[e] = counts.get(e, 0) + 1
        order = (["C", "H"] if "C" in counts else []) + sorted(k for k in counts if k not in ("C", "H") or "C" not in counts)
        return "".join(f"{k}{counts[k] if counts[k] > 1 else ''}" for k in order if k in counts)

    def summary(self) -> dict:
        return {
            "source": self.source, "sha256": self.sha256, "kind": self.kind, "n_atoms": self.natoms,
            "formula": self.formula, "charge": self.charge, "multiplicity": self.multiplicity,
            "route": self.route, "program": self.program, "termination": self.termination,
            "mass_source": self.mass_source, "oniom_layers": sorted(set(self.layers)) if self.layers else None,
            "n_frozen": len(self.frozen), "dipole_derivatives": self.dipole_derivatives is not None,
            "anharmonic_present": self.anharmonic_present, "notes": list(self.notes),
            "attachments": {k: str(v) for k, v in self.attachments.items()},
        }
