"""PEDRO: Potential Energy Distribution with Resonances and Overtones.

Harmonic potential- and total-energy distributions (PED/TED) from Cartesian Hessians, for isolated molecules
(VEDA-compatible) and for a molecule embedded in a much larger Hessian (QM/MM, ONIOM, partial Hessians).

    import pedro
    res = pedro.analyze("my_freq_job.log")   # reads, chooses a basis, computes TED
    res.write()                                     # new results folder next to the input
"""

__version__ = "0.2.0b1"

from .io import load  # noqa: E402
from .pipeline import Options, analyze  # noqa: E402

__all__ = ["__version__", "load", "analyze", "Options"]
