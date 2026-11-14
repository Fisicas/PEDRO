"""Load any supported input into a Job. The format is chosen from the file name and contents."""

from __future__ import annotations

from pathlib import Path

from ..job import Job
from .common import ReadError

SUPPORTED = (".log", ".out", ".fchk", ".fch", ".fmu", ".fmv", ".hess", ".npz", ".json")


def load(path) -> Job:
    p = Path(path)
    if not p.exists():
        raise ReadError(f"{p}: file not found")
    suf = p.suffix.lower()
    if suf in (".log", ".out"):
        from . import gaussian_log
        job = gaussian_log.load(p)
    elif suf in (".fchk", ".fch"):
        from . import fchk
        job = fchk.load(p, "fchk")
        _attach_log(job, p)
    elif suf == ".fmu":
        from . import fchk
        job = fchk.load(p, "veda_fmu")
    elif suf == ".fmv":
        from . import veda
        job = veda.load_fmv(p)
    elif suf == ".hess":
        from . import orca
        job = orca.load(p)
    elif suf in (".npz", ".json"):
        from . import generic
        job = generic.load(p)
    else:
        raise ReadError(f"{p.name}: unsupported file type. PEDRO reads {', '.join(SUPPORTED)}")
    if job.masses is None:
        from ..constants import ISOTOPE_MASS
        import numpy as np
        job.masses = np.array([ISOTOPE_MASS[e] for e in job.elements])
        job.mass_source = "default most-abundant isotopes (file has no masses)"
    for ext in (".dd2", ".ved"):
        sib = p.with_suffix(ext)
        if sib.exists():
            job.attachments[ext[1:]] = sib
    return job


def _attach_log(job: Job, fchk_path: Path):
    """Use printed frequencies and ONIOM layers from a sibling log, when one exists, for validation."""
    for ext in (".log", ".out"):
        sib = fchk_path.with_suffix(ext)
        if sib.exists():
            from . import gaussian_log
            try:
                lg = gaussian_log.load(sib)
            except ReadError:
                return
            if lg.elements == job.elements:
                job.printed_freqs, job.printed_ir = lg.printed_freqs, lg.printed_ir
                job.layers = job.layers or lg.layers
                job.route, job.program, job.termination = lg.route, lg.program, lg.termination
                job.attachments["log"] = sib
            return


__all__ = ["load", "ReadError", "SUPPORTED"]
