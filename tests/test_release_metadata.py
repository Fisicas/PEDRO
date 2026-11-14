"""Prevent version drift between the installed package, citation and public install commands."""

from pathlib import Path
import re
import tomllib

import pedro


def test_beta_metadata_and_install_tag_agree():
    root = Path(__file__).resolve().parents[1]
    citation = (root / "CITATION.cff").read_text(encoding="utf-8")
    assert re.search(r'^version: "([^"]+)"$', citation, re.M).group(1) == pedro.__version__
    assert f"@v{pedro.__version__}" in (root / "README.md").read_text(encoding="utf-8")
    project = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
    assert "Development Status :: 4 - Beta" in project["project"]["classifiers"]
    assert "ONIOM" in project["project"]["description"]
