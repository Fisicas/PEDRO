"""PEDRO desktop window: drop a Gaussian .log/.fchk, VEDA .fmu or .npz/.json Hessian and get the full analysis.

Run with `pedro` or `pedro gui` (needs the optional 'gui' extras: PySide6 and matplotlib).
"""

from __future__ import annotations

import os
import sys
import traceback
from pathlib import Path

import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets
from PySide6.QtCore import Qt

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from .. import __version__  # noqa: E402
from ..io import SUPPORTED  # noqa: E402
from ..pipeline import Options, analyze  # noqa: E402
from ..report import environment_identity_check  # noqa: E402

ACCENT = "#0f6e74"
ASSETS = Path(__file__).resolve().parent / "assets"


# ----------------------------------------------------------------------------- background worker
class Worker(QtCore.QObject):
    progress = QtCore.Signal(str)
    finished = QtCore.Signal(object, object)   # result, output folder
    failed = QtCore.Signal(str)

    def __init__(self, path, options, out_dir):
        super().__init__()
        self.path, self.options, self.out_dir = path, options, out_dir

    @QtCore.Slot()
    def run(self):
        try:
            res = analyze(self.path, self.options, progress=self.progress.emit)
            self.progress.emit("Writing tables, coordinate files and the HTML report…")
            out = res.write(self.out_dir or None)
            self.progress.emit(f"Results saved in {out}")
            self.finished.emit(res, out)
        except Exception as e:  # shown to the user, with the full traceback in the log
            self.progress.emit(traceback.format_exc())
            self.failed.emit(f"{type(e).__name__}: {e}")


# ----------------------------------------------------------------------------- widgets
class DropZone(QtWidgets.QFrame):
    fileDropped = QtCore.Signal(str)

    def __init__(self):
        super().__init__()
        self.setAcceptDrops(True)
        self.setMinimumHeight(170)
        self.setObjectName("drop")
        lay = QtWidgets.QVBoxLayout(self)
        lay.setAlignment(Qt.AlignCenter)
        logo = ASSETS / "logo.svg"
        if logo.exists():
            from PySide6.QtSvgWidgets import QSvgWidget
            w = QSvgWidget(str(logo))
            w.setFixedSize(72, 72)
            lay.addWidget(w, alignment=Qt.AlignCenter)
        self.label = QtWidgets.QLabel("Drop a frequency calculation here\n"
                                      "Gaussian .log / .out / .fchk · VEDA .fmu · .npz / .json Hessian")
        self.label.setAlignment(Qt.AlignCenter)
        self.label.setStyleSheet("color:#44505e;font-size:13px")
        lay.addWidget(self.label)
        btn = QtWidgets.QPushButton("Choose file…")
        btn.clicked.connect(self.browse)
        lay.addWidget(btn, alignment=Qt.AlignCenter)
        self._idle()

    def _idle(self):
        self.setStyleSheet("#drop{border:2px dashed #9aa6b2;border-radius:12px;background:#f7f9fb}")

    def browse(self):
        pat = " ".join("*" + s for s in SUPPORTED)
        f, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open a frequency calculation", "", f"Hessian files ({pat})")
        if f:
            self.fileDropped.emit(f)

    def dragEnterEvent(self, e):
        if e.mimeData().hasUrls():
            e.acceptProposedAction()
            self.setStyleSheet(f"#drop{{border:2px solid {ACCENT};border-radius:12px;background:#e8f4f4}}")

    def dragLeaveEvent(self, e):
        self._idle()

    def dropEvent(self, e):
        self._idle()
        urls = [u.toLocalFile() for u in e.mimeData().urls() if u.isLocalFile()]
        if urls:
            self.fileDropped.emit(urls[0])


class Plot(FigureCanvasQTAgg):
    def __init__(self, h=3.0):
        self.fig = Figure(figsize=(6, h), tight_layout=True)
        super().__init__(self.fig)


class ModeModel(QtCore.QAbstractTableModel):
    HEAD = ["#", "cm⁻¹", "scaled", "IR", "S-TED %", "S-amp %", "Description (TED %)"]

    def __init__(self):
        super().__init__()
        self.res, self.rows = None, []

    def set(self, res, rows):
        self.beginResetModel()
        self.res, self.rows = res, list(rows)
        self.endResetModel()

    def rowCount(self, parent=QtCore.QModelIndex()):
        return len(self.rows)

    def columnCount(self, parent=QtCore.QModelIndex()):
        return len(self.HEAD)

    def headerData(self, s, o, role=Qt.DisplayRole):
        if role == Qt.DisplayRole and o == Qt.Horizontal:
            return self.HEAD[s]

    def data(self, idx, role=Qt.DisplayRole):
        k = self.rows[idx.row()]
        r, c = self.res, idx.column()
        if role == Qt.DisplayRole:
            ir = r.modes.ir_km_mol
            return [str(k + 1), f"{r.modes.freqs[k]:.2f}", f"{r.scaled[k]:.2f}", "" if ir is None else f"{ir[k]:.1f}",
                    f"{r.subsystem_ted[k]:.0f}", f"{r.subsystem_amplitude[k]:.0f}", r.descriptions[k]][c]
        if role == Qt.ForegroundRole and r.modes.freqs[k] < 0:
            return QtGui.QColor("#a33")
        if role == Qt.TextAlignmentRole and c < 6:
            return int(Qt.AlignRight | Qt.AlignVCenter)


# ----------------------------------------------------------------------------- main window
class Main(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"PEDRO {__version__} — energy distributions from Hessians")
        icon = ASSETS / "logo.svg"
        if icon.exists():
            self.setWindowIcon(QtGui.QIcon(str(icon)))
        self.resize(1280, 820)
        self.res = None
        self.out = None
        self.thread = None

        left = QtWidgets.QWidget()
        lv = QtWidgets.QVBoxLayout(left)
        self.drop = DropZone()
        self.drop.fileDropped.connect(self.start)
        lv.addWidget(self.drop)
        lv.addWidget(self._options())
        self.log = QtWidgets.QPlainTextEdit(readOnly=True)
        self.log.setPlaceholderText("Progress appears here.")
        self.log.setStyleSheet("font-family:Consolas,monospace;font-size:11px")
        lv.addWidget(self.log, 1)
        self.open_btn = QtWidgets.QPushButton("Open results folder")
        self.open_btn.setEnabled(False)
        self.open_btn.clicked.connect(lambda: self.out and os.startfile(self.out))
        self.report_btn = QtWidgets.QPushButton("Open HTML report")
        self.report_btn.setEnabled(False)
        self.report_btn.clicked.connect(lambda: self.out and os.startfile(Path(self.out) / "report.html"))
        hb = QtWidgets.QHBoxLayout()
        hb.addWidget(self.open_btn)
        hb.addWidget(self.report_btn)
        lv.addLayout(hb)
        left.setMaximumWidth(430)

        self.tabs = QtWidgets.QTabWidget()
        self.summary = QtWidgets.QTextBrowser()
        self.summary.setOpenExternalLinks(True)
        self.summary.setHtml(self._welcome())
        self.tabs.addTab(self.summary, "Summary")
        self.tabs.addTab(self._modes_tab(), "Modes")
        self.spec = Plot(4.2)
        self.spec.mpl_connect("pick_event", self._picked)
        self.tabs.addTab(self.spec, "Spectrum")
        self.tabs.addTab(self._anh_tab(), "Resonances && overtones")

        split = QtWidgets.QSplitter()
        split.addWidget(left)
        split.addWidget(self.tabs)
        split.setStretchFactor(1, 1)
        self.setCentralWidget(split)
        self.statusBar().showMessage("Ready. Drop a file to start.")

    # ---------------------------------------------------------------- panels
    def _options(self):
        box = QtWidgets.QGroupBox("Options (defaults work for most jobs)")
        f = QtWidgets.QFormLayout(box)
        self.coords = QtWidgets.QLineEdit(placeholderText="automatic (or a .dd2 next to the input)")
        cb = QtWidgets.QToolButton(text="…")
        cb.clicked.connect(lambda: self._pick(self.coords, "Coordinate set (*.dd2 *.json)"))
        f.addRow("Coordinates", self._row(self.coords, cb))
        self.amap = QtWidgets.QLineEdit(placeholderText="only to place a .dd2 on a larger model")
        mb = QtWidgets.QToolButton(text="…")
        mb.clicked.connect(lambda: self._pick(self.amap, "Atom map (*.csv)"))
        f.addRow("Atom map", self._row(self.amap, mb))
        self.subsys = QtWidgets.QLineEdit(placeholderText="automatic, e.g. 1-31,40")
        f.addRow("Subsystem atoms", self.subsys)
        self.case = QtWidgets.QComboBox()
        self.case.addItems(["auto", "isolated", "F", "Fprime", "P"])
        self.case.setToolTip("auto: isolated molecule, or case F (subsystem internals + environment Cartesians) "
                             "for large/ONIOM models.\nFprime: unprojected check.\nP: partial Hessian.")
        f.addRow("Basis", self.case)
        self.mobile = QtWidgets.QLineEdit(placeholderText="case P only: layer:H or atom list")
        f.addRow("Mobile atoms", self.mobile)
        self.scale = QtWidgets.QDoubleSpinBox(decimals=4, minimum=0.5, maximum=1.5, singleStep=0.001, value=1.0)
        self.scale.setToolTip("Display only; the analysis always uses unscaled frequencies.")
        f.addRow("Frequency scale", self.scale)
        self.outdir = QtWidgets.QLineEdit(placeholderText="new folder next to the input")
        ob = QtWidgets.QToolButton(text="…")
        ob.clicked.connect(lambda: self._pick_dir(self.outdir))
        f.addRow("Output", self._row(self.outdir, ob))
        return box

    def _row(self, *ws):
        w = QtWidgets.QWidget()
        h = QtWidgets.QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        for x in ws:
            h.addWidget(x)
        return w

    def _pick(self, edit, filt):
        f, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Choose file", "", filt)
        if f:
            edit.setText(f)

    def _pick_dir(self, edit):
        d = QtWidgets.QFileDialog.getExistingDirectory(self, "Output folder (a new subfolder is created inside)")
        if d:
            edit.setText(d)

    def _modes_tab(self):
        w = QtWidgets.QWidget()
        v = QtWidgets.QVBoxLayout(w)
        bar = QtWidgets.QHBoxLayout()
        self.only_s = QtWidgets.QCheckBox("Only modes with ≥ 20% amplitude on the subsystem", checked=True)
        self.only_s.toggled.connect(self._fill_modes)
        self.wmin = QtWidgets.QDoubleSpinBox(minimum=-5000, maximum=5000, value=-5000, decimals=0, prefix="from ")
        self.wmax = QtWidgets.QDoubleSpinBox(minimum=-5000, maximum=5000, value=5000, decimals=0, prefix="to ", suffix=" cm⁻¹")
        for s in (self.wmin, self.wmax):
            s.valueChanged.connect(self._fill_modes)
        bar.addWidget(self.only_s)
        bar.addStretch()
        bar.addWidget(QtWidgets.QLabel("Scaled window"))
        bar.addWidget(self.wmin)
        bar.addWidget(self.wmax)
        v.addLayout(bar)
        split = QtWidgets.QSplitter(Qt.Vertical)
        self.model = ModeModel()
        self.table = QtWidgets.QTableView()
        self.table.setModel(self.model)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.verticalHeader().setVisible(False)
        self.table.setAlternatingRowColors(True)
        self.table.selectionModel().currentRowChanged.connect(lambda cur, _: self._show_mode(cur.row()))
        split.addWidget(self.table)
        self.detail = Plot(3.2)
        split.addWidget(self.detail)
        split.setSizes([420, 300])
        v.addWidget(split)
        return w

    def _anh_tab(self):
        w = QtWidgets.QWidget()
        v = QtWidgets.QVBoxLayout(w)
        self.anh_info = QtWidgets.QLabel("Run an analysis to see overtones, combination bands and resonances.")
        self.anh_info.setWordWrap(True)
        v.addWidget(self.anh_info)
        bar = QtWidgets.QHBoxLayout()
        self.anh_strong = QtWidgets.QCheckBox("Only fundamentals and bands with I(anh) ≥ 1 km/mol", checked=True)
        self.anh_strong.toggled.connect(self._fill_anh)
        self.anh_kind = QtWidgets.QComboBox()
        self.anh_kind.addItems(["all kinds", "fundamental", "overtone", "combination"])
        self.anh_kind.currentIndexChanged.connect(self._fill_anh)
        self.anh_verdict = QtWidgets.QComboBox()
        self.anh_verdict.addItems(["any verdict", "clean", "warning", "failed"])
        self.anh_verdict.currentIndexChanged.connect(self._fill_anh)
        for x in (self.anh_strong, self.anh_kind, self.anh_verdict):
            bar.addWidget(x)
        bar.addStretch()
        v.addLayout(bar)
        self.anh_table = QtWidgets.QTableWidget()
        self.anh_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.anh_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.anh_table.setAlternatingRowColors(True)
        self.anh_table.verticalHeader().setVisible(False)
        self.anh_table.setWordWrap(True)
        v.addWidget(self.anh_table, 1)
        return w

    def _fill_anh(self):
        r = self.res
        if r is None:
            return
        a, sc = r.anharmonic or {}, r.resonance_screen
        t = self.anh_table
        t.clearContents()
        col = {"clean": "#1d7a3a", "warning": "#a86b00", "failed": "#a33"}
        if a.get("states"):
            self.anh_info.setText(
                f"<b>GVPT2</b> ({a.get('pt2_model')}): {a['counts']['fundamental']} fundamentals, "
                f"{a['counts']['overtone']} overtones, {a['counts']['combination']} combination bands · verdicts "
                f"{a['verdicts']} · symmetry gauge: {a['symmetry'].get('group') or a['symmetry'].get('status')}.<br>"
                f"<span style='color:#5b6675'>{a['caveat']}</span>")
            for x in (self.anh_strong, self.anh_kind, self.anh_verdict):
                x.setEnabled(True)
            rows = [s for s in a["states"]
                    if (not self.anh_strong.isChecked() or s["kind"] == "fundamental" or (s["i_anharm_km_mol"] or 0) >= 1)
                    and self.anh_kind.currentText() in ("all kinds", s["kind"])
                    and self.anh_verdict.currentText() in ("any verdict", s["verdict"])]
            rows.sort(key=lambda s: -s["e_anharm_cm1"])
            head = ["State", "E(harm)", "E(anh)", "I(anh)", "Verdict", "Parentage (harmonic TED of parents)",
                    "Resonances", "Flags"]
            t.setColumnCount(len(head))
            t.setHorizontalHeaderLabels(head)
            t.setRowCount(len(rows))
            for i, s in enumerate(rows):
                vals = [s["pedro_quanta"], f"{s['e_harm_cm1']:.1f}", f"{s['e_anharm_cm1']:.1f}",
                        f"{s['i_anharm_km_mol']:.2f}", s["verdict"], s["parentage"], "; ".join(s["resonances"]),
                        "; ".join(s["flags"])]
                for j, val in enumerate(vals):
                    it = QtWidgets.QTableWidgetItem(val)
                    if j == 4:
                        it.setForeground(QtGui.QColor(col[s["verdict"]]))
                    t.setItem(i, j, it)
        elif sc:
            self.anh_info.setText(
                f"<b>No GVPT2 data</b> ({a.get('reason')}). Harmonic screen: {sc['n_candidates']} Fermi-resonance "
                f"candidates for {sc['n_fundamentals_screened']} IR-active fundamentals of {sc['n_modes_screened']} subsystem modes (|ω<sub>k</sub> − ω<sub>i</sub> − ω<sub>j</sub>| ≤ "
                f"{sc['delta_cm1']:g} cm⁻¹; scaled frequencies shown).<br><span style='color:#5b6675'>{sc['caveat']}</span>")
            for x in (self.anh_strong, self.anh_kind, self.anh_verdict):
                x.setEnabled(False)
            sf = r.options.scale
            head = ["Fundamental", "scaled cm⁻¹", "IR", "Partner", "scaled cm⁻¹", "Δ (cm⁻¹)",
                    "Fundamental character", "Partner character"]
            t.setColumnCount(len(head))
            t.setHorizontalHeaderLabels(head)
            t.setRowCount(len(sc["candidates"]))
            for i, c in enumerate(sc["candidates"]):
                vals = [f"ν{c['fundamental']}", f"{c['omega_k_cm1'] * sf:.1f}",
                        "" if c["ir_k_km_mol"] is None else f"{c['ir_k_km_mol']:.1f}", c["partner"],
                        f"{c['partner_harmonic_cm1'] * sf:.1f}", f"{c['delta_cm1']:+.1f}",
                        c["fundamental_character"], c["partner_character"]]
                for j, val in enumerate(vals):
                    t.setItem(i, j, QtWidgets.QTableWidgetItem(val))
        else:
            self.anh_info.setText(f"Resonance/overtone layer: {a.get('status')} ({a.get('reason')}).")
            t.setRowCount(0)
        t.resizeColumnsToContents()
        for j in range(t.columnCount()):
            t.setColumnWidth(j, min(t.columnWidth(j), 420))

    def _welcome(self):
        return f"""<h2 style='color:{ACCENT}'>PEDRO {__version__}</h2>
        <p>Public beta: harmonic mode assignments for molecules embedded in large ONIOM/QM/MM models.</p>
        <p>Drop a completed frequency calculation on the left. PEDRO reads the Cartesian Hessian, recomputes the normal modes,
        checks them against the frequencies the program printed, builds an internal-coordinate basis and reports the
        total energy distribution (TED) of every mode. For an isolated molecule it also gives the VEDA-compatible PED.
        For a large or ONIOM model it describes the embedded molecule in its own internal coordinates and reports the
        environment's share.</p>
        <p>Nothing is clipped, renormalized or dropped. Imaginary modes are kept and shown in red.</p>
        <p>Every run writes a new folder (tables, coordinate files, manifest, HTML report). Existing results are never
        overwritten.</p>"""

    # ---------------------------------------------------------------- run
    def _parse_atoms(self, s):
        from ..cli import _atoms
        return _atoms(s.strip()) if s.strip() else None

    def start(self, path):
        if self.thread is not None:
            QtWidgets.QMessageBox.information(self, "PEDRO", "An analysis is already running.")
            return
        try:
            mob = self.mobile.text().strip()
            opts = Options(coords=self.coords.text().strip() or None, atom_map=self.amap.text().strip() or None,
                           subsystem=self._parse_atoms(self.subsys.text()), case=self.case.currentText(),
                           mobile=mob if mob.startswith("layer:") else self._parse_atoms(mob), scale=self.scale.value())
        except ValueError as e:
            QtWidgets.QMessageBox.warning(self, "PEDRO", f"Could not read an atom list: {e}")
            return
        out = None
        if self.outdir.text().strip():
            from datetime import date
            base = Path(self.outdir.text().strip())
            for v in range(1, 1000):
                out = base / f"{Path(path).stem}_PEDRO_{date.today():%Y%m%d}_v{v:02d}"
                if not out.exists():
                    break
        self.log.clear()
        self.log.appendPlainText(f"Starting {path}")
        self.statusBar().showMessage("Working…")
        self.drop.label.setText(f"Analyzing {Path(path).name}…")
        self.open_btn.setEnabled(False)
        self.report_btn.setEnabled(False)
        self.thread = QtCore.QThread()
        self.worker = Worker(path, opts, out)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.progress.connect(self.log.appendPlainText)
        self.worker.finished.connect(self.done)
        self.worker.failed.connect(self.fail)
        self.worker.finished.connect(self.thread.quit)
        self.worker.failed.connect(self.thread.quit)
        self.thread.finished.connect(self._cleanup)
        self.thread.start()

    def _cleanup(self):
        self.thread.deleteLater()
        self.thread = None

    def fail(self, msg):
        self.statusBar().showMessage("Failed")
        self.drop.label.setText("Drop a frequency calculation here")
        QtWidgets.QMessageBox.critical(self, "PEDRO could not finish", msg)

    def done(self, res, out):
        self.res, self.out = res, str(out)
        self.drop.label.setText(f"{Path(res.job.source).name} — done. Drop another file to start again.")
        self.open_btn.setEnabled(True)
        self.report_btn.setEnabled(True)
        self.statusBar().showMessage(f"Done in {res.timings['total_s']:.1f} s · results in {out}")
        self.only_s.setChecked(res.basis.case != "isolated")
        self.summary.setHtml(self._summary_html(res, out))
        self._fill_modes()
        self._draw_spectrum()
        self._fill_anh()
        self.tabs.setCurrentIndex(0)

    # ---------------------------------------------------------------- views
    def _badge(self, ok, text):
        col = {True: "#1d7a3a", False: "#a33", None: "#a86b00"}[ok]
        sym = {True: "✔", False: "✖", None: "!"}[ok]
        return f"<span style='color:{col};font-weight:600'>{sym}</span> {text}"

    def _summary_html(self, res, out):
        j, c = res.job, res.checks
        rec = c.get("reconstruction")
        items = []
        if rec:
            items.append(self._badge(rec["max_abs_cm1"] < 0.5, f"Frequencies reproduce the printed values "
                                     f"(max difference {rec['max_abs_cm1']:.4f} cm⁻¹ over {rec['n_printed']} modes)"))
        else:
            items.append(self._badge(None, "No printed frequencies to compare with"))
        items.append(self._badge(c["n_imaginary"] == 0, f"{c['n_imaginary']} imaginary mode(s) — kept and flagged"))
        st = c.get("stationarity")
        if st:
            items.append(self._badge(st["max_atom_force_Eh_bohr"] < 4.5e-4,
                                     f"Largest atomic force {st['max_atom_force_Eh_bohr']:.1e} Eh/bohr "
                                     f"(subsystem {st['max_atom_force_subsystem_Eh_bohr']:.1e})"))
        items.append(self._badge(c["ted_closure_max_pp"] < 1e-6, f"TED closes to 100% for every mode "
                                 f"(max error {c['ted_closure_max_pp']:.0e} pp)"))
        identity = environment_identity_check(c)
        if identity:
            # F is an approximate convention, so a nonzero deviation is information, not a failed identity.
            exact = "identity_max_dev_pp" in c
            ok = c["identity_max_dev_pp"] < 1e-6 if exact else None
            items.append(self._badge(ok, f"{identity[0]}: {identity[1]}"))
        items.append(self._badge(c["condition_L"] < 1e8, f"Coordinate basis well conditioned (cond L = {c['condition_L']:.3g})"))
        v = res.validation
        if v and v.get("status") == "compared":
            items.append(self._badge(v["ted_fraction_within_1pp"] > 0.98,
                                     f"VEDA {v['veda_version']}: {100 * v['ted_fraction_within_1pp']:.2f}% of TED entries "
                                     f"within 1 pp; PED {100 * v.get('ped_fraction_within_1pp', 0):.2f}%"))
        a = res.anharmonic or {}
        if a.get("states"):
            items.append(self._badge(a["status"] == "ok", f"GVPT2: {a['n_states']} states read and labeled "
                                     f"({a['verdicts']['clean']} clean, {a['verdicts']['warning']} warning, "
                                     f"{a['verdicts']['failed']} failed); mapping residual "
                                     f"{a['mapping']['max_residual_cm1']:.3f} cm⁻¹"))
        elif a.get("status") == "failed":
            items.append(self._badge(False, f"GVPT2 layer failed closed: {a.get('reason')} (harmonic results unaffected)"))
        if res.resonance_screen:
            items.append(self._badge(None, f"No GVPT2 data: harmonic screen lists {res.resonance_screen['n_candidates']} "
                                     "Fermi-resonance candidates (see Resonances & overtones)"))
        notes = "".join(f"<li>{n}</li>" for n in j.notes)
        return f"""<h2 style='color:{ACCENT};margin-bottom:2px'>{Path(j.source).name}</h2>
        <p style='color:#5b6675'>{j.natoms} atoms · {j.formula} · charge {j.charge} · multiplicity {j.multiplicity}<br>
        {j.route}<br>{j.program}</p>
        <h3>Checks</h3><p>{'<br>'.join(items)}</p>
        <h3>Analysis</h3><p>Basis case <b>{c['case']}</b>: {res.basis.n_internal} internal coordinates on
        {len(res.basis.subsystem)} subsystem atoms ({c['subsystem_choice']}); {len(res.basis.labels)} coordinates,
        {res.modes.n} modes.<br>Masses: {j.mass_source}.<br>Time: {res.timings['total_s']:.1f} s.</p>
        {'<ul>' + notes + '</ul>' if notes else ''}
        {self._fragments_html(c)}
        <p>Results: <a href='file:///{out}'>{out}</a></p>"""

    def _fragments_html(self, c):
        fr = c.get("fragments")
        if not fr or len(fr) < 2:
            return ""
        rows = "".join(f"<tr><td>{f['n']}</td><td>{f['formula']}</td><td>{f['ranges']}</td><td>{f['tag']}</td></tr>"
                       for f in fr)
        return ("<h3>Bonded fragments</h3><p>To describe another fragment, paste its atom numbers into "
                "<i>Subsystem atoms</i> and drop the file again.</p><table cellpadding='3'><tr><th>atoms</th>"
                f"<th>formula</th><th>atom numbers</th><th>type</th></tr>{rows}</table>")

    def _fill_modes(self):
        r = self.res
        if r is None:
            return
        s = r.scaled
        keep = (s >= self.wmin.value()) & (s <= self.wmax.value())
        if self.only_s.isChecked():
            keep &= r.subsystem_amplitude >= 20
        rows = np.where(keep)[0][::-1]
        self.model.set(r, rows)
        self.table.resizeColumnsToContents()
        if len(rows):
            self.table.selectRow(0)

    def _show_mode(self, row):
        r = self.res
        if r is None or row < 0 or row >= len(self.model.rows):
            return
        k = self.model.rows[row]
        fig = self.detail.fig
        fig.clear()
        ax1, ax2 = fig.subplots(1, 2, gridspec_kw={"width_ratios": [3, 2]})
        ni = r.basis.n_internal
        t = r.ted[k, :ni]
        idx = np.argsort(-np.abs(t))[:12][::-1]
        ax1.barh(range(idx.size), t[idx], color=[ACCENT if v >= 0 else "#c0504d" for v in t[idx]])
        ax1.set_yticks(range(idx.size), [r.basis.labels[i] for i in idx], fontsize=8)
        ax1.axvline(0, color="k", lw=0.5)
        ax1.set_xlabel("TED %")
        ax1.set_title(f"Mode {k + 1}: {r.modes.freqs[k]:.1f} cm⁻¹ (scaled {r.scaled[k]:.1f})", fontsize=9)
        g = {name: v[k] for name, v in r.group_totals.items()}
        names = sorted(g, key=lambda n: -abs(g[n]))[:10]
        ax2.barh(range(len(names))[::-1], [g[n] for n in names], color="#7a8fa6")
        ax2.set_yticks(range(len(names))[::-1], names, fontsize=8)
        ax2.axvline(0, color="k", lw=0.5)
        ax2.set_xlabel("group total %")
        self.detail.draw_idle()

    def _draw_spectrum(self):
        r = self.res
        fig = self.spec.fig
        fig.clear()
        ax = fig.add_subplot(111)
        f = r.scaled
        ir = r.modes.ir_km_mol if r.modes.ir_km_mol is not None else np.ones(r.modes.n)
        ax.vlines(f, 0, ir, color="#c5ccd4", lw=0.7)
        sc = ax.scatter(f, ir, c=np.clip(r.subsystem_ted, 0, 100), cmap="viridis", vmin=0, vmax=100, s=16, picker=5,
                        zorder=3)
        fig.colorbar(sc, ax=ax, label="subsystem TED %")
        ax.set_xlabel("scaled wavenumber / cm⁻¹")
        ax.set_ylabel("IR intensity / km mol⁻¹")
        ax.invert_xaxis()
        ax.set_title("Click a point to open that mode", fontsize=9)
        self.spec.draw_idle()

    def _picked(self, ev):
        if self.res is None or not len(ev.ind):
            return
        k = int(ev.ind[0])
        self.only_s.setChecked(False)
        self.wmin.setValue(-5000)
        self.wmax.setValue(5000)
        rows = list(self.model.rows)
        if k in rows:
            self.tabs.setCurrentIndex(1)
            self.table.selectRow(rows.index(k))


def main(argv=None):
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv if argv is None else argv)
    app.setStyle("Fusion")
    app.setApplicationName("PEDRO")
    if sys.platform == "win32":
        app.setFont(QtGui.QFont("Segoe UI", 9))
    w = Main()
    w.show()
    args = (sys.argv if argv is None else argv)[1:]
    files = [a for a in args if Path(a).suffix.lower() in SUPPORTED and Path(a).exists()]
    if files:
        QtCore.QTimer.singleShot(200, lambda: w.start(files[0]))
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
