"""Exercise the installed GUI's worker, reports and tabs without opening a visible window.

Run from a checkout after installing PEDRO, for example:
python scripts/smoke_gui.py --input examples/acetone_2h2o_b3lyp_orca.hess --outdir /new/output
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--case", default="auto", choices=["auto", "isolated", "F", "Fprime", "P"])
    ap.add_argument("--coords", type=Path)
    ap.add_argument("--subsystem")
    ap.add_argument("--expected-modes", type=int)
    ap.add_argument("--screenshot", type=Path)
    ap.add_argument("--screenshot-tab", type=int, default=0, choices=range(4), help="tab captured after all tabs are checked")
    ap.add_argument("--timeout", type=int, default=120, help="seconds")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=False)
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6 import QtCore, QtGui, QtWidgets
    import matplotlib
    from pedro import __version__
    from pedro.gui.app import Main

    failures = []
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(["pedro-gui-smoke"])
    app.setStyle("Fusion")
    # Offscreen Windows Qt may not enumerate system fonts in a sandbox. Use the GUI dependency's own font.
    font = Path(matplotlib.get_data_path()) / "fonts" / "ttf" / "DejaVuSans.ttf"
    QtGui.QFontDatabase.addApplicationFont(str(font))
    app.setFont(QtGui.QFont("DejaVu Sans", 9))

    class SmokeWindow(Main):
        def fail(self, message):
            # A failure dialog would block unattended CI; record the original worker error instead.
            failures.append(message)

    window = SmokeWindow()
    window.resize(1400, 860)
    window.case.setCurrentText(args.case)
    window.outdir.setText(str(args.outdir.resolve()))
    if args.coords:
        window.coords.setText(str(args.coords.resolve()))
    if args.subsystem:
        window.subsys.setText(args.subsystem)
    window.show()

    def unhandled(kind, value, traceback):
        failures.append(f"{kind.__name__}: {value}")
        sys.__excepthook__(kind, value, traceback)

    sys.excepthook = unhandled

    def inspect():
        if window.thread is not None:
            QtCore.QTimer.singleShot(50, inspect)
            return
        if failures or window.res is None:
            if not failures:
                failures.append("GUI returned no result")
            app.exit(1)
            return
        result = window.res
        try:
            if args.expected_modes is not None:
                assert result.modes.n == args.expected_modes, result.modes.n
            summary = window.summary.toPlainText()
            if result.basis.case == "F":
                assert "case F; approximate" in summary
            elif result.basis.case in ("Fprime", "P"):
                assert "exact identity" in summary
            for index in range(window.tabs.count()):
                window.tabs.setCurrentIndex(index)
                app.processEvents()
            if window.model.rowCount():
                window.tabs.setCurrentIndex(1)
                window.table.selectRow(0)
                app.processEvents()
            for name in ("modes.csv", "ted_subsystem.csv", "manifest.json", "report.html"):
                assert (Path(window.out) / name).is_file(), name
            if args.screenshot:
                window.tabs.setCurrentIndex(args.screenshot_tab)
                app.processEvents()
                args.screenshot.parent.mkdir(parents=True, exist_ok=True)
                assert not args.screenshot.exists(), "screenshot already exists"
                assert window.grab().save(str(args.screenshot))
            print(json.dumps({"version": __version__, "case": result.basis.case,
                              "modes": result.modes.n, "output": window.out,
                              "tabs_checked": window.tabs.count(), "status": "passed"}), flush=True)
            app.exit(0)
        except Exception as error:
            failures.append(str(error))
            app.exit(1)

    def start():
        window.start(str(args.input.resolve()))
        QtCore.QTimer.singleShot(50, inspect)

    def timeout():
        failures.append("GUI smoke check timed out")
        app.exit(2)

    QtCore.QTimer.singleShot(0, start)
    QtCore.QTimer.singleShot(args.timeout * 1000, timeout)
    code = app.exec()
    if window.thread is not None:
        window.thread.quit()
        window.thread.wait(30000)
    window.close()
    if failures:
        print("\n".join(failures), file=sys.stderr)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
