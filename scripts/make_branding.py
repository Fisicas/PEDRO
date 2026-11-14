"""PEDRO branding: logo (option C, 1:2 Fermi-resonance figure) in light, dark and adaptive versions, wordmark
lockups, and GitHub social previews built from a real calculation.

    python scripts/make_branding.py [GVPT2.log "caption"] [--out DIR]

The optional Gaussian Freq=Anharmonic log supplies the social-preview spectrum: each fundamental drawn as a stacked
bar of its TED shares by coordinate class, GVPT2 overtones and combination bands as thin lines, and Fermi-resonant
fundamentals marked. Without it a synthetic spectrum is used. SVGs go to branding/, PNGs are rendered with Qt.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "branding"
FONT = "Segoe UI, Helvetica Neue, Arial, sans-serif"

THEMES = {
    "light": {"bg": "#f7f9fb", "ink": "#1b2430", "mute": "#5b6675", "tick": "#7a8fa6", "grid": "#d9dee5",
              "a": "#0f6e74", "b": "#e09f3e",
              "cls": {"ν(C–O)": "#c0504d", "ν(C–C)": "#0f6e74", "ν(C–H)": "#6a5acd", "δ (bend)": "#2f6fd6",
                      "γ/τ": "#7a8fa6"}, "line": "#d08a1e"},
    "dark": {"bg": "#0e1621", "ink": "#e6edf3", "mute": "#9aa7b5", "tick": "#8aa0b8", "grid": "#25303d",
             "a": "#3fb6bc", "b": "#f0b450",
             "cls": {"ν(C–O)": "#ff7a6b", "ν(C–C)": "#3fb6bc", "ν(C–H)": "#a99be6", "δ (bend)": "#6fa8ff",
                     "γ/τ": "#8aa0b8"}, "line": "#f0b450"},
}


def svg(w, h, body, defs="", style=""):
    st = f"<style>{style}</style>" if style else ""
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">{st}'
            f"<defs>{defs}</defs>{body}</svg>\n")


# ----------------------------------------------------------------------------- the mark (1:2 Lissajous)
def lobes():
    """Two lobes of x = sin t, y = sin 2t as smooth paths (right lobe t in [0, pi], left lobe t in [pi, 2 pi])."""
    def path(t0, t1):
        t = np.linspace(t0, t1, 181)
        x, y = 256 + 150 * np.sin(t), 256 - 110 * np.sin(2 * t)
        return "M" + " L".join(f"{a:.2f},{b:.2f}" for a, b in zip(x, y))
    return path(0, math.pi), path(math.pi, 2 * math.pi)


def mark(theme=None, prefix="m"):
    """theme=None gives the adaptive version (colors from CSS classes switched by prefers-color-scheme)."""
    right, left = lobes()
    if theme is None:
        ca, cb, ring, tick = "var(--pa)", "var(--pb)", "var(--ring)", "var(--tick)"
    else:
        T = THEMES[theme]
        ca, cb, ring, tick = T["a"], T["b"], T["ink"], T["tick"]
    defs = (f'<linearGradient id="{prefix}R" gradientUnits="userSpaceOnUse" x1="0" y1="146" x2="0" y2="366">'
            f'<stop offset="0" style="stop-color:{ca}"/><stop offset="1" style="stop-color:{cb}"/></linearGradient>'
            f'<linearGradient id="{prefix}L" gradientUnits="userSpaceOnUse" x1="0" y1="146" x2="0" y2="366">'
            f'<stop offset="0" style="stop-color:{cb}"/><stop offset="1" style="stop-color:{ca}"/></linearGradient>')
    body = (f'<circle cx="256" cy="256" r="206" fill="none" style="stroke:{ring}" stroke-width="12"/>'
            + "".join(f'<line x1="{a}" y1="{b}" x2="{c}" y2="{d}" style="stroke:{tick}" stroke-width="8" stroke-linecap="round"/>'
                      for a, b, c, d in ((70, 256, 100, 256), (412, 256, 442, 256), (256, 70, 256, 100), (256, 412, 256, 442)))
            + f'<path d="{left}" fill="none" stroke="url(#{prefix}L)" stroke-width="16" stroke-linecap="round"/>'
            + f'<path d="{right}" fill="none" stroke="url(#{prefix}R)" stroke-width="16" stroke-linecap="round"/>')
    return defs, body


ADAPTIVE_CSS = (":root{--ring:#1b2430;--tick:#7a8fa6;--pa:#0f6e74;--pb:#e09f3e;--word:#1b2430;--sub:#5b6675}"
                "@media (prefers-color-scheme: dark){:root{--ring:#e6edf3;--tick:#8aa0b8;--pa:#3fb6bc;--pb:#f0b450;"
                "--word:#e6edf3;--sub:#9aa7b5}}")


def logo(theme):
    d, b = mark(theme)
    return svg(512, 512, b, d, ADAPTIVE_CSS if theme is None else "")


def lockup(theme):
    d, b = mark(theme)
    word = "var(--word)" if theme is None else THEMES[theme]["ink"]
    sub = "var(--sub)" if theme is None else THEMES[theme]["mute"]
    body = (f'<g transform="scale(0.5)">{b}</g>'
            f'<text x="276" y="160" font-family="{FONT}" font-weight="700" font-size="104" letter-spacing="6" '
            f'style="fill:{word}">PEDRO</text>'
            f'<text x="280" y="204" font-family="{FONT}" font-size="25" style="fill:{sub}">energy distributions · '
            f'resonances · overtones</text>')
    return svg(820, 256, body, d, ADAPTIVE_CSS if theme is None else "")


# ----------------------------------------------------------------------------- social preview data
def spectrum_data(log, caption=None):
    if log:
        sys.path.insert(0, str(ROOT / "src"))
        import pedro
        r = pedro.analyze(log)
        a = r.anharmonic
        ni = r.basis.n_internal
        cls_of = []
        for g in r.basis.groups[:ni]:
            cls_of.append("γ/τ" if g.startswith(("γ", "τ")) else g if g in ("ν(C–O)", "ν(C–C)", "ν(C–H)", "δ (bend)")
                          else "δ (bend)" if g.startswith("δ") else "γ/τ")
        fund, extra = [], []
        for s in a["states"]:
            if not 980 <= s["e_anharm_cm1"] <= 1790:
                continue
            if s["kind"] == "fundamental":
                k = int(next(iter(s["pedro_modes"]))) - 1
                t = np.abs(r.ted[k, :ni])
                shares = {}
                for c, v in zip(cls_of, t):
                    shares[c] = shares.get(c, 0) + v
                tot = sum(shares.values())
                fund.append({"x": s["e_anharm_cm1"], "I": s["i_anharm_km_mol"], "fermi": any(x.startswith("Fermi") for x in s["resonances"]),
                             "shares": {c: v / tot for c, v in shares.items()}, "label": r.descriptions[k]})
            elif s["verdict"] != "failed" and s["i_anharm_km_mol"] >= 0.5:
                extra.append({"x": s["e_anharm_cm1"], "I": s["i_anharm_km_mol"]})
        name = caption or Path(log).stem
        return fund, extra, name, a["n_states"]
    rng = np.random.default_rng(3)
    fund = [{"x": x, "I": rng.uniform(1, 200), "fermi": rng.random() < 0.3,
             "shares": dict(zip(THEMES["light"]["cls"], rng.dirichlet(np.ones(5)))), "label": ""}
            for x in np.linspace(1000, 1750, 22)]
    extra = [{"x": x, "I": rng.uniform(0.5, 20)} for x in rng.uniform(1000, 1780, 60)]
    return fund, extra, "synthetic example", 0


def social(theme, data):
    T = THEMES[theme]
    fund, extra, name, nstates = data
    W, H = 1280, 640
    x0, x1, base, top = 640, 1210, 470, 120
    lo, hi = 980, 1790
    X = lambda v: x0 + (hi - v) / (hi - lo) * (x1 - x0)
    imax = max(f["I"] for f in fund)
    Y = lambda i: (base - top) * math.sqrt(min(max(i, 0), imax) / imax)   # extra bands capped at the tallest fundamental
    d, mk = mark(theme, "s")
    b = [f'<rect width="{W}" height="{H}" fill="{T["bg"]}"/>']
    # axis and grid
    for v in range(1000, 1800, 100):
        b.append(f'<line x1="{X(v):.1f}" y1="{top - 10}" x2="{X(v):.1f}" y2="{base}" stroke="{T["grid"]}" stroke-width="1"/>')
        b.append(f'<text x="{X(v):.1f}" y="{base + 26}" font-family="{FONT}" font-size="16" fill="{T["mute"]}" text-anchor="middle">{v}</text>')
    b.append(f'<line x1="{x0}" y1="{base}" x2="{x1}" y2="{base}" stroke="{T["mute"]}" stroke-width="1.5"/>')
    b.append(f'<text x="{(x0 + x1) / 2}" y="{base + 50}" font-family="{FONT}" font-size="16" fill="{T["mute"]}" '
             f'text-anchor="middle">wavenumber / cm⁻¹ (anharmonic)</text>')
    for e in extra:
        b.append(f'<line x1="{X(e["x"]):.1f}" y1="{base}" x2="{X(e["x"]):.1f}" y2="{base - Y(e["I"]):.1f}" '
                 f'stroke="{T["line"]}" stroke-width="2" stroke-opacity="0.75"/>')
    order = [c for c in T["cls"] if any(f["shares"].get(c, 0) > 0.02 for f in fund)]
    for f in fund:
        h, y = Y(f["I"]), base
        for c in order:
            seg = h * f["shares"].get(c, 0)
            if seg <= 0:
                continue
            b.append(f'<rect x="{X(f["x"]) - 3.5:.1f}" y="{y - seg:.1f}" width="7" height="{seg:.1f}" fill="{T["cls"][c]}"/>')
            y -= seg
        if f["fermi"]:
            b.append(f'<circle cx="{X(f["x"]):.1f}" cy="{base - h - 12:.1f}" r="4.5" fill="none" stroke="{T["line"]}" stroke-width="2"/>')
    # legend
    lx, ly = x0, 560
    for c in order:
        b.append(f'<rect x="{lx}" y="{ly - 13}" width="14" height="14" fill="{T["cls"][c]}"/>')
        b.append(f'<text x="{lx + 20}" y="{ly}" font-family="{FONT}" font-size="16" fill="{T["ink"]}">{c}</text>')
        lx += 20 + 9.2 * len(c) + 22
    b.append(f'<line x1="{x0}" y1="{ly + 22}" x2="{x0}" y2="{ly + 36}" stroke="{T["line"]}" stroke-width="3"/>'
             f'<text x="{x0 + 10}" y="{ly + 35}" font-family="{FONT}" font-size="16" fill="{T["ink"]}">overtones and combination bands</text>'
             f'<circle cx="{x0 + 283}" cy="{ly + 29}" r="4.5" fill="none" stroke="{T["line"]}" stroke-width="2"/>'
             f'<text x="{x0 + 294}" y="{ly + 35}" font-family="{FONT}" font-size="16" fill="{T["ink"]}">Fermi resonance</text>')
    b.append(f'<text x="{x0}" y="78" font-family="{FONT}" font-size="17" fill="{T["mute"]}">{name}: each band split by its TED</text>')
    # left column
    b.append(f'<g transform="translate(64,58) scale(0.25)">{mk}</g>')
    b.append(f'<text x="208" y="148" font-family="{FONT}" font-weight="700" font-size="92" letter-spacing="5" fill="{T["ink"]}">PEDRO</text>')
    b.append(f'<text x="66" y="236" font-family="{FONT}" font-size="30" fill="{T["ink"]}">Potential- and total-energy distributions</text>')
    b.append(f'<text x="66" y="274" font-family="{FONT}" font-size="30" fill="{T["ink"]}">with resonances and overtones</text>')
    stats = [("±0.5 pp", "TED agreement with VEDA on identical input"),
             ("614 atoms", "ONIOM Hessian analyzed in about 4 s"),
             (f"{nstates:,} states" if nstates else "GVPT2", "labeled from one Gaussian .log")]
    for k, (big, small) in enumerate(stats):
        y = 360 + 78 * k
        b.append(f'<rect x="66" y="{y - 36}" width="5" height="48" fill="{[T["a"], T["b"], T["cls"]["ν(C–O)"]][k]}"/>')
        b.append(f'<text x="86" y="{y - 4}" font-family="{FONT}" font-weight="700" font-size="30" fill="{T["ink"]}">{big}</text>')
        b.append(f'<text x="86" y="{y + 20}" font-family="{FONT}" font-size="18" fill="{T["mute"]}">{small}</text>')
    return svg(W, H, "".join(b), d)


def render_png(svg_path, png_path, w, h):
    from PySide6 import QtCore, QtGui, QtSvg
    QtGui.QGuiApplication.instance() or QtGui.QGuiApplication([])
    r = QtSvg.QSvgRenderer(str(svg_path))
    img = QtGui.QImage(w, h, QtGui.QImage.Format_ARGB32)
    img.fill(QtCore.Qt.transparent)
    p = QtGui.QPainter(img)
    p.setRenderHint(QtGui.QPainter.Antialiasing)
    r.render(p)
    p.end()
    img.save(str(png_path))


def main():
    global OUT
    args = sys.argv[1:]
    if "--out" in args:
        i = args.index("--out")
        OUT = Path(args[i + 1])
        del args[i:i + 2]
    OUT.mkdir(exist_ok=True, parents=True)
    data = spectrum_data(args[0] if args else None, args[1] if len(args) > 1 else None)
    (OUT / "logo.svg").write_text(logo(None), encoding="utf-8")
    (OUT / "lockup.svg").write_text(lockup(None), encoding="utf-8")
    for th in ("light", "dark"):
        (OUT / f"logo_{th}.svg").write_text(logo(th), encoding="utf-8")
        (OUT / f"lockup_{th}.svg").write_text(lockup(th), encoding="utf-8")
        (OUT / f"social_preview_{th}.svg").write_text(social(th, data), encoding="utf-8")
        render_png(OUT / f"logo_{th}.svg", OUT / f"logo_{th}.png", 512, 512)
        render_png(OUT / f"lockup_{th}.svg", OUT / f"lockup_{th}.png", 1640, 512)
        render_png(OUT / f"social_preview_{th}.svg", OUT / f"social_preview_{th}.png", 1280, 640)
    print("written to", OUT)


if __name__ == "__main__":
    main()
