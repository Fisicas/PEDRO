"""Run Qt in a separate process so its plotting backend cannot affect numerical tests."""

import importlib.util
import os
from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None or importlib.util.find_spec("matplotlib") is None,
                    reason="GUI extras not installed")
def test_installed_gui_embedded_partition_summary(tmp_path):
    root = Path(__file__).resolve().parents[1]
    env = dict(os.environ, QT_QPA_PLATFORM="offscreen")
    process = subprocess.run([sys.executable, "-X", "utf8", str(root / "scripts" / "smoke_gui.py"),
                              "--input", str(root / "examples" / "acetone_2h2o_b3lyp_orca.hess"),
                              "--case", "P", "--expected-modes", "30", "--outdir", str(tmp_path / "gui")],
                             cwd=tmp_path, env=env, capture_output=True, text=True, encoding="utf-8", timeout=150)
    assert process.returncode == 0, process.stdout + process.stderr
    assert '"status": "passed"' in process.stdout
