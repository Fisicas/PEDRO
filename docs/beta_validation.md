# 0.2.0b1 beta preparation checks

Local checks performed on Windows 11, 13 November 2026. These are preparation results; no public GitHub CI run,
release tag or archive publication is implied.

| Environment | Python | NumPy | SciPy | PySide6 | pytest | Result |
|---|---|---|---|---|---|---|
| Existing science environment | 3.12.10 | 2.5.3 | 1.18.1 | 6.11.2 | 9.1.1 | 58 passed, 3 optional private GVPT2 references skipped |
| Fresh environment, built wheel and newly resolved GUI dependencies | 3.12.10 | 2.5.4 | 1.18.1 | 6.12.0 | 9.2.0 | 58 passed, same 3 skips |

Matplotlib was 3.11.2 in both environments. The six tests using VEDA's distributed examples ran against an existing
local VEDA installation; those third-party files are not redistributed. Without that installation, those tests
skip. Core reference tests use owned/synthetic inputs included in the repository.

The installer integration test uses an offline local Git repository: first install, no-op update, a newer
annotated beta tag, version pinning, a final release, unpinning, an unavailable version and offline fallback launch.
It verifies that failed updates preserve the active installation. Existing dependencies are reused only for this
offline installer test; the separate clean wheel installation resolved dependencies into a genuinely new venv.
GitHub bootstrap URLs require the actual repository and first version tag to be published before a live check.

The GUI was exercised with the public acetone/water example in P and the owned 614-atom ONIOM model in F. The
large-model check explicitly selected atoms 1–14 and transferred a public benzaldehyde coordinate set by graph
isomorphism; all four tabs, mode selection and report writing were exercised. Offscreen captures were inspected;
this is not a test of every interactive desktop operation or file-dialog path. README screenshots were regenerated
from the owned public examples, with local paths removed.

The 614-atom F result was compared with the preserved dev4 result using the same automatic coordinate labels.
Maximum differences were 2.30e-10 cm^-1 in frequency and 8.21e-9 percentage points in TED, consistent with numerical
roundoff between runs. No Hessian, mode-analysis or coordinate-construction algorithm was changed by the reporting
fixes. Signed values and ascending numerical mode order were retained; no spectral smoothing, baseline correction,
normalization, interpolation or resampling was introduced.

The numerical [validation report](validation.md) was regenerated into a new JSON record, preserving the previous
record. The [portable UniMoVib cross-check](../validation/own/unimovib_crosscheck/README.md) was rerun on the same
input hash with the recorded printed-mass convention and pseudo-inverse cutoff.

Wheel and source archives build successfully. Package contents and a candidate Git file list were inspected:
owned examples/reference inputs and the installer are included; vendor executables, checkpoints, the large ONIOM
log and external/private reference calculations are excluded. Original calculation inputs and prior validation
records were not modified. Full file hashes, test logs, package hashes and local GUI outputs are retained in the
maintainer's versioned preparation folder outside the repository.

Still pending: remote Windows/Linux CI, a live install/update from the public GitHub tags, the public benchmark
archive link, and a representative protein/cofactor application demonstration. The archive and application are
not asserted by these software checks. See the [release checklist](RELEASE_CHECKLIST.md).
