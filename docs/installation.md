# Installation and GitHub synchronization

The Windows managed installer is `Install-PEDRO.ps1`. It needs Windows PowerShell 5.1+, Git and Python 3.11+.
The README gives the bootstrap commands. A downloaded repository can use the same installer directly:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\Install-PEDRO.ps1 -Python C:/path/to/python.exe
```

The script reads version tags from `https://github.com/Fisicas/PEDRO.git` and selects the highest matching
`vMAJOR.MINOR.PATCH`, `aN`, `bN` or `rcN` version. The default beta channel includes pre-releases and final releases;
the stable channel includes final releases only. It follows published Git tags, not untagged development commits
or GitHub's stable-only `releases/latest` endpoint. Maintainers should create release tags only for reviewed builds.

Each installation uses `pip` with the exact selected Git revision and GUI extras in a new private environment.
Only after version and GUI dependency imports pass does the script switch `install-state.json` to that environment.
An existing application process continues using its original environment. A failed update retains the previous
active environment, and the launcher starts that copy if the update service is unavailable. Installation errors
are displayed by the launcher. No calculation input or result is uploaded; normal network operations read Git
tags/source from GitHub and retrieve dependencies from the configured Python package index.

`install-config.json` records the interpreter, repository, channel and optional version pin. The generated
`PEDRO.cmd` runs the updater before opening the GUI. `-Action Update` performs only the update; `-SkipUpdate`
starts the existing version without a network check. `-Version 0.2.0b1` pins the managed copy until `-Unpin` is used.
There is no background timer, service or scheduled task. Manual pip/source installations update only when asked.

The GUI entry point is also available in each private environment as `Scripts/pedro-gui.exe`; the CLI is
`Scripts/pedro.exe`. Find the active environment in `install-state.json`. Direct entry points bypass the updater.
The run manifest includes pip's `direct_url.json`, whose `vcs_info.commit_id` identifies the installed Git commit.
Dependencies and numerical conventions are recorded separately in the manifest; a version pin alone does not
guarantee identical results across different numerical libraries or coordinate choices.

For manual source installs, use `python -m pip install ".[gui]"` in the checkout. For a particular public tag:

```powershell
python -m pip install --upgrade "pedro-ped[gui] @ git+https://github.com/Fisicas/PEDRO@v0.2.0b1"
```

The managed installer does not modify or pull into a user's development checkout. Its environments and update
state live under `%LOCALAPPDATA%\PEDRO` by default; research files should be kept in their own folders. Previous
version environments are retained for diagnosis or deliberate reuse. The `-NoDeps` and `-NoBuildIsolation` switches
are intended for offline packaging/testing with dependencies already available; ordinary users should omit them.
