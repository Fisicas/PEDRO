"""Regenerate the README screenshots from the public example calculations in examples/.

    python scripts/make_screenshots.py

Opens the desktop window, analyzes each example, switches to the tab shown in the README and saves
docs/images/gui_<name>.png. Output folders are written to a temporary directory, never next to the examples.
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from PySide6 import QtCore, QtGui, QtWidgets  # noqa: E402

from pedro.gui.app import Main  # noqa: E402

SHOTS = [  # (output name, input, tab index, scale, scaled window)
    ("modes", ROOT / "examples" / "acetone_2h2o_b3lyp_orca.hess", 1, 0.965, (1150, 1800)),
    ("summary", ROOT / "examples" / "acetone_2h2o_b3lyp_orca.hess", 0, 0.965, (1550, 1800)),
    ("resonances", ROOT / "examples" / "benzaldehyde_gvpt2_b3lyp_gaussian.log", 3, 1.0, (400, 4000)),
]


def main():
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
    app.setStyle("Fusion")
    if os.environ.get("QT_QPA_PLATFORM") == "offscreen":
        import matplotlib
        font = Path(matplotlib.get_data_path()) / "fonts" / "ttf" / "DejaVuSans.ttf"
        QtGui.QFontDatabase.addApplicationFont(str(font))
        app.setFont(QtGui.QFont("DejaVu Sans", 9))
    elif sys.platform == "win32":
        app.setFont(QtGui.QFont("Segoe UI", 9))
    tmp = tempfile.TemporaryDirectory()
    queue = list(SHOTS)
    w = Main()
    w.resize(1400, 860)
    w.show()

    def next_shot():
        if not queue:
            app.quit()
            return
        name, path, tab, scale, (lo, hi) = queue[0]
        w.scale.setValue(scale)
        w.wmin.setValue(lo)
        w.wmax.setValue(hi)
        w.outdir.setText(tmp.name)
        w.res = None
        w.start(str(path))
        QtCore.QTimer.singleShot(500, wait)

    def wait():
        # poll from the GUI thread: the analysis runs in a worker thread, and the window must only be
        # touched from here
        if w.thread is not None or w.res is None:
            QtCore.QTimer.singleShot(300, wait)
            return
        name, path, tab, scale, win = queue.pop(0)
        w.tabs.setCurrentIndex(tab)

        def snap():
            hide_paths()
            app.processEvents()
            out = ROOT / "docs" / "images" / f"gui_{name}.png"
            w.grab().save(str(out))
            print("saved", out, flush=True)
            QtCore.QTimer.singleShot(300, next_shot)

        QtCore.QTimer.singleShot(800, snap)

    def hide_paths():
        """Keep local folder names out of the published images."""
        local = [str(Path(tmp.name)) + "\\", str(Path(tmp.name)) + "/", str(Path(tmp.name)), str(Path(tmp.name).resolve()), str(ROOT / "examples") + "\\", str(ROOT / "examples") + "/"]

        def clean(t):
            for a in local:
                t = t.replace(a, "")
            return t.replace(str(Path(tmp.name)), "")

        w.log.setPlainText(clean(w.log.toPlainText()))
        w.summary.setHtml(clean(w.summary.toHtml()))
        w.statusBar().showMessage(clean(w.statusBar().currentMessage()))
        w.outdir.setText("")

    QtCore.QTimer.singleShot(300, next_shot)
    app.exec()
    tmp.cleanup()


if __name__ == "__main__":
    main()
