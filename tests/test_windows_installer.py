"""Offline installation/update/pinning against a local Git remote, including an annotated beta tag."""

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import pytest


@pytest.mark.skipif(sys.platform != "win32" or shutil.which("git") is None or shutil.which("powershell") is None,
                    reason="Windows PowerShell and Git required")
def test_managed_install_update_pin_and_offline_launch(tmp_path):
    if not all(__import__("importlib.util", fromlist=["find_spec"]).find_spec(n) for n in ("PySide6", "matplotlib", "setuptools")):
        pytest.skip("offline installer integration uses existing GUI/build dependencies")
    remote = tmp_path / "remote"
    remote.mkdir()
    env = dict(os.environ)
    # Dependencies are deliberately reused for this offline test only; the installer defaults to an isolated venv.
    helper = tmp_path / "dependency_path"
    helper.mkdir()
    # Append dependency access after the new venv's own packages, including when CI has PEDRO installed.
    (helper / "sitecustomize.py").write_text(
        f"import site\nsite.addsitedir({str(Path(sys.prefix) / 'Lib' / 'site-packages')!r})\n", encoding="utf-8")
    env["PYTHONPATH"] = str(helper)
    env["PIP_NO_INDEX"] = "1"
    marker = tmp_path / "launched.txt"
    env["PEDRO_TEST_LAUNCH_MARKER"] = str(marker)

    def git(*args):
        result = subprocess.run(["git", *args], cwd=remote, env=env, capture_output=True, text=True, encoding="utf-8")
        assert result.returncode == 0, result.stdout + result.stderr
        return result.stdout.strip()

    git("init")
    git("config", "user.name", "PEDRO installer test")
    git("config", "user.email", "pedro-test@example.invalid")
    package = remote / "pedro"
    package.mkdir()
    (package / "__main__.py").write_text("import os\nfrom pathlib import Path\nPath(os.environ['PEDRO_TEST_LAUNCH_MARKER']).write_text('launched')\n", encoding="utf-8")

    def publish(version, annotated=False):
        (package / "__init__.py").write_text(f'__version__ = "{version}"\n', encoding="utf-8")
        (remote / "pyproject.toml").write_text(
            '[build-system]\nrequires=["setuptools>=68"]\nbuild-backend="setuptools.build_meta"\n'
            f'[project]\nname="pedro-ped"\nversion="{version}"\n[project.optional-dependencies]\ngui=[]\n', encoding="utf-8")
        git("add", ".")
        git("commit", "-m", version)
        git("tag", *( ["-a", f"v{version}", "-m", version] if annotated else [f"v{version}"] ))
        return git("rev-parse", "HEAD")

    first = publish("0.2.0b1")
    install = tmp_path / "managed install with spaces"
    script = Path(__file__).resolve().parents[1] / "Install-PEDRO.ps1"

    def manage(action, *args, success=True):
        result = subprocess.run(["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", str(script),
                                 "-Action", action, "-InstallDir", str(install), "-NoDeps", "-NoBuildIsolation", *args],
                                cwd=tmp_path, env=env, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=150)
        assert (result.returncode == 0) == success, result.stdout + result.stderr
        return json.loads((install / "install-state.json").read_text(encoding="utf-8-sig"))

    state = manage("Install", "-Python", sys.executable, "-Repository", remote.as_uri())
    assert state["Revision"] == first and state["Version"] == "0.2.0b1"
    old_environment = state["Environment"]
    assert (install / "PEDRO.cmd").is_file()
    # An unchanged tag must not create another environment.
    assert manage("Update")["Environment"] == old_environment
    second = publish("0.2.0b2", annotated=True)
    state = manage("Update")
    assert state["Revision"] == second and state["Version"] == "0.2.0b2"
    assert Path(old_environment).is_dir()
    assert manage("Update", "-Version", "0.2.0b1")["Revision"] == first
    stable = publish("0.2.0")
    assert manage("Update")["Revision"] == first  # pin survives a new release
    assert manage("Update", "-Unpin")["Revision"] == stable
    # An unavailable requested version cannot replace the working installation.
    assert manage("Update", "-Version", "9.0.0", success=False)["Revision"] == stable
    # Make the remote unavailable; Launch still starts the installed package without changing the pointer.
    unavailable = tmp_path / "remote unavailable"
    remote.rename(unavailable)
    state = manage("Launch")
    assert state["Revision"] == stable
    for _ in range(100):
        if marker.is_file():
            break
        time.sleep(0.05)
    assert marker.read_text() == "launched"
