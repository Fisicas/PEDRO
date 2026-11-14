# Contributing and support

PEDRO is a research-software beta focused on harmonic mode assignments for molecules embedded in large ONIOM/QM/MM
models. Development will continue publicly through GitHub issues and pull requests. The maintainer is Michael
Raymond Nelson. Support is best effort; there is no guaranteed response time.

Use American English for documentation, user-facing text, release notes, communications and publication drafts.
Preserve exact published titles, quotations, program output and technical identifiers when reproducing sources.

For a bug report, include the PEDRO version, operating system, Python/NumPy/SciPy versions, the command or GUI
options, observed behavior, expected behavior and traceback. Attach a minimal input only if you have permission
to share it. A manifest can help, but inspect it for local paths and unpublished sample information first.
Describe units, atom numbering, coordinate set, case F/Fprime/P, mass choices and display scale explicitly.

Discuss a substantial method or file-reader change in an issue before implementing it. Keep numerical conventions
and limitations documented. Tests should check a known answer, independent implementation or scientific invariant;
passing closure alone does not establish a correct assignment. New fixtures must have clear provenance and
redistribution permission. Do not commit vendor executables, third-party example data or unpublished research logs.

From a checkout with Python 3.11 or later:

```powershell
python -m pip install -e ".[test,gui]"
python -m pytest -q
```

The core tests use owned/public or synthetic fixtures. Additional references are optional:
`PEDRO_VEDA_EXAMPLES` points to a local VEDA installation and `PEDRO_GVPT2_DATA` to private reference calculations.
Those skips are expected without the external inputs. GUI smoke tests skip when GUI dependencies are absent;
the CI GUI job installs those dependencies and checks the built wheel.

For a pull request, explain the problem, resulting behavior, numerical conventions affected and checks run.
Include focused regression coverage for scientific or output-accounting corrections. Preserve input files and
previous results. Contributions are under the repository's MIT license; identify borrowed code and its license.
Disclose substantial generated-code assistance in the pull-request description and review the code and its
scientific assumptions yourself.
