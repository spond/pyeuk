"""
Graphical report generation for the PyEuk cluster-count sweep.

This module turns the output of ``CyclosporaClusterFinder.cluster_sweep()`` -- either
the returned dict or the ``*_SWEEP.json`` it writes -- into a single, self-contained
HTML string. There is no JavaScript anywhere: charts are inline SVG, and the optional
distance heatmap is an embedded PNG.

Design constraints (deliberate, see the module tests):

* **Dependency-light.** Only ``numpy`` is required at import time; ``Pillow`` (PIL) is
  imported lazily and only for the heatmap -- if it is missing the heatmap is skipped
  with a note rather than crashing. The module has no relative imports and never pulls
  in ``numba`` or the heavy package siblings, so it can be exercised standalone.
* **Portable bytes.** The document declares ``<meta charset="utf-8">`` and the returned
  string is guaranteed to contain **zero raw non-ASCII bytes** -- every non-ASCII code
  point is emitted as an HTML entity / numeric character reference. This renders
  correctly even when a bare server sends no charset header.
* **Two themes.** ``studio`` uses Fraunces/Inter via a Google Fonts link (fine for
  standalone viewing); ``galaxy`` uses system fonts and the Galaxy brand palette and
  contains **no external assets at all** (no ``http(s)://`` references), so it renders
  embedded inside Galaxy with nothing to block.
* **Three flavors.** ``dashboard`` (default), ``clinical``, ``narrative``. A confident
  cohort is reported as a single green number; a fuzzy one as an amber range. The stable
  cores are drawn on the tree itself as numbered bars.

Public entry point::

    html = render(sweep, dist_df=None, flavor="dashboard", theme="studio")
"""

import base64
import colorsys
import io
import json

import numpy as np

__all__ = [
    "render",
    "dendrogram_svg",
    "sweep_curve_svg",
    "heatmap_img",
    "cores_rows",
    "FLAVORS",
    "THEMES",
]

FLAVORS = ("dashboard", "clinical", "narrative")
THEMES = ("studio", "galaxy")


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def _esc(text):
    """HTML-escape a plain-text string for safe embedding as element text."""
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _ascii_only(html):
    """Return ``html`` with every non-ASCII code point replaced by a numeric
    character reference, so the emitted document has zero raw non-ASCII bytes."""
    return html.encode("ascii", "xmlcharrefreplace").decode("ascii")


def _load_sweep(sweep):
    """Accept a sweep result dict or a path to a ``*_SWEEP.json`` file."""
    if isinstance(sweep, (str, bytes)):
        with open(sweep, "r") as fh:
            return json.load(fh)
    if isinstance(sweep, dict):
        return sweep
    raise TypeError("sweep must be a dict or a path to a *_SWEEP.json file")


def _load_matrix(dist_df, leaf_order):
    """Return a tree-ordered square distance array, or ``None`` if no usable
    matrix is available.

    ``dist_df`` may be ``None``, a path to a CSV/TSV distance matrix, or a
    pandas-like object exposing ``.index`` and ``.values``. Ids are matched to
    the tree ``leaf_order``; a trailing ``B`` batch/background suffix on the
    matrix ids is tolerated when it improves the overlap.
    """
    if dist_df is None:
        return None

    ids, M = _matrix_ids_values(dist_df)
    if M is None or not len(ids):
        return None

    M = np.asarray(M, dtype=float)
    M = (M + M.T) / 2.0
    np.fill_diagonal(M, 0.0)

    order = _match_order(ids, leaf_order)
    if len(order) < 2:
        return None
    return M[np.ix_(order, order)]


def _matrix_ids_values(dist_df):
    """Extract (ids, 2-D array) from a path or a pandas-like DataFrame."""
    # pandas-like: has .index and .values
    if hasattr(dist_df, "index") and hasattr(dist_df, "values"):
        ids = [str(x) for x in list(dist_df.index)]
        return ids, np.asarray(dist_df.values, dtype=float)

    if isinstance(dist_df, (str, bytes)):
        path = dist_df.decode() if isinstance(dist_df, bytes) else dist_df
        with open(path, "r") as fh:
            text = fh.read()
        rows = [ln for ln in text.splitlines() if ln.strip()]
        if not rows:
            return [], None
        sep = "\t" if ("\t" in rows[0] or path.lower().endswith((".tsv", ".tab"))) else ","
        header = rows[0].split(sep)
        ids = [h.strip() for h in header[1:]]
        data = []
        for ln in rows[1:]:
            cells = ln.split(sep)[1:]
            data.append([float(c) if c.strip() != "" else 0.0 for c in cells])
        return ids, np.asarray(data, dtype=float)

    raise TypeError("dist_df must be None, a file path, or a DataFrame-like object")


def _match_order(ids, leaf_order):
    """Return matrix-row indices reordered to follow ``leaf_order``.

    Prefer exact id matches; if those cover fewer than half the leaves, fall
    back to normalising a single trailing ``B`` suffix on the matrix ids.
    """
    idx = {s: i for i, s in enumerate(ids)}
    exact = [idx[s] for s in leaf_order if s in idx]
    if len(exact) >= 0.5 * len(leaf_order):
        return exact
    norm = {}
    for i, s in enumerate(ids):
        norm.setdefault(s, i)
        if s.endswith("B"):
            norm.setdefault(s[:-1], i)
    return [norm[s] for s in leaf_order if s in norm]


# ---------------------------------------------------------------------------
# Component renderers (pure inline SVG / embedded PNG)
# ---------------------------------------------------------------------------
def dendrogram_svg(tree, w=320, h=150, cores=None, corelabels=False):
    """Confidence tree: leaves on x, merge height on y, branch style by support.

    If ``cores`` (a list of specimen-id lists) is given, each stable core is
    drawn as a bar beneath its leaves -- the reproduced groups, shown on the
    tree itself; ``corelabels`` numbers them to match the cores table.
    """
    nodes = tree.get("nodes", [])
    hmax = tree.get("hmax", 0.0) or 0.0
    order = tree.get("leaf_order", [])
    n = len(order)
    mb = 22 if (cores and corelabels) else 16 if cores else 6
    ml, mr, mt = 6, 6, 8
    pw, ph = w - ml - mr, h - mt - mb
    X = (lambda p: ml + (p / (n - 1)) * pw) if n > 1 else (lambda p: ml + pw / 2)
    Y = (lambda hh: mt + (1 - hh / hmax) * ph) if hmax > 0 else (lambda hh: mt + ph)
    seg = []
    for nd in nodes:
        s = nd.get("support", 0.0)
        if s >= 0.75:
            st = "stroke:var(--ink);stroke-width:1.7;opacity:1"
        elif s >= 0.45:
            st = "stroke:var(--ink2);stroke-width:1.2;opacity:.72"
        else:
            st = "stroke:var(--muted);stroke-width:1;opacity:.4;stroke-dasharray:2.5 2.5"
        xa, xb, yh = X(nd["ya"]), X(nd["yb"]), Y(nd["h"])
        seg.append(f'<line x1="{xa:.1f}" x2="{xb:.1f}" y1="{yh:.1f}" y2="{yh:.1f}" style="{st}"/>')
        seg.append(f'<line x1="{xa:.1f}" x2="{xa:.1f}" y1="{yh:.1f}" y2="{Y(nd["hca"]):.1f}" style="{st}"/>')
        seg.append(f'<line x1="{xb:.1f}" x2="{xb:.1f}" y1="{yh:.1f}" y2="{Y(nd["hcb"]):.1f}" style="{st}"/>')
    if cores and n > 0:
        pos = {s: i for i, s in enumerate(order)}
        lw = (pw / (n - 1)) if n > 1 else pw
        by = Y(0) + 5
        bh = 4
        for ci, core in enumerate(sorted(cores, key=len, reverse=True)):
            idxs = sorted(pos[s] for s in core if s in pos)
            if not idxs:
                continue
            runs = []
            a = p = idxs[0]
            for q in idxs[1:]:
                if q == p + 1:
                    p = q
                else:
                    runs.append((a, p))
                    a = p = q
            runs.append((a, p))
            op = 0.92 if ci % 2 == 0 else 0.55  # alternate tint so adjacent cores separate
            for x0i, x1i in runs:
                x0 = X(x0i) - lw * 0.42
                x1 = X(x1i) + lw * 0.42
                seg.append(
                    f'<rect x="{x0:.1f}" y="{by:.1f}" width="{max(1.2, x1 - x0):.1f}" '
                    f'height="{bh}" rx="1.6" fill="var(--blue)" opacity="{op}"/>'
                )
            if corelabels:
                cx = (X(idxs[0]) + X(idxs[-1])) / 2
                seg.append(
                    f'<text x="{cx:.1f}" y="{by + bh + 7:.1f}" text-anchor="middle" '
                    f'style="font-family:var(--fm);font-size:7px;fill:var(--blue);font-weight:600">{ci + 1}</text>'
                )
    return (
        f'<svg viewBox="0 0 {w} {h}" width="100%" style="max-width:{w}px" role="img" '
        f'aria-label="Confidence tree with stable cores marked beneath the leaves.">'
        f'{"".join(seg)}</svg>'
    )


def sweep_curve_svg(r, w=360, h=180):
    """Silhouette and bootstrap stability vs k, with the supported-range band."""
    sw = list(r.get("sweep", []))
    ks = [x["k"] for x in sw if x.get("k") is not None]
    ml, mr, mt, mb = 34, 10, 10, 22
    pw, ph = w - ml - mr, h - mt - mb
    if len(ks) < 2:
        return (
            f'<svg viewBox="0 0 {w} {h}" width="100%" style="max-width:{w}px" role="img" '
            f'aria-label="Sweep curve unavailable.">'
            f'<text x="{w / 2:.0f}" y="{h / 2:.0f}" text-anchor="middle" class="ax">'
            f'sweep curve unavailable</text></svg>'
        )
    kmin, kmax = min(ks), max(ks)
    span = (kmax - kmin) or 1
    X = lambda k: ml + (k - kmin) / span * pw
    Y = lambda v: mt + (1 - v) * ph

    def line(key, col):
        pts = [(X(x["k"]), Y(x[key])) for x in sw if x.get(key) is not None]
        if not pts:
            return ""
        d = "M" + " L".join(f"{a:.1f},{b:.1f}" for a, b in pts)
        return f'<path d="{d}" fill="none" style="stroke:{col};stroke-width:1.8;opacity:.9"/>'

    lo, hi = r.get("count_range", [kmin, kmax])
    lo = max(kmin, min(lo, kmax))
    hi = max(kmin, min(hi, kmax))
    band = (
        f'<rect x="{X(lo):.1f}" y="{mt}" width="{max(0.0, X(hi) - X(lo)):.1f}" '
        f'height="{ph}" fill="var(--blue)" opacity="0.08"/>'
    )
    grid = "".join(
        f'<line x1="{ml}" x2="{w - mr}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" '
        f'stroke="var(--line)" stroke-width=".5"/>'
        for v in (0, 0.5, 1)
    )
    ax = (
        f'<text x="{ml - 4}" y="{Y(1) + 3:.0f}" text-anchor="end" class="ax">1.0</text>'
        f'<text x="{ml - 4}" y="{Y(0) + 3:.0f}" text-anchor="end" class="ax">0</text>'
        f'<text x="{X(lo):.0f}" y="{h - 6}" class="ax" fill="var(--blue)">{lo}</text>'
        f'<text x="{X(hi):.0f}" y="{h - 6}" class="ax" fill="var(--blue)">{hi}</text>'
        f'<text x="{w - mr}" y="{h - 6}" text-anchor="end" class="ax">k &rarr;</text>'
    )
    return (
        f'<svg viewBox="0 0 {w} {h}" width="100%" style="max-width:{w}px" role="img" '
        f'aria-label="Silhouette and stability across candidate cluster counts.">'
        f'{band}{grid}{line("silhouette", "var(--blue)")}{line("stability", "var(--good)")}{ax}</svg>'
    )


def heatmap_img(D, w=240):
    """Tree-ordered pairwise distance heatmap as an embedded PNG.

    Returns an ``<img>`` tag on success, or a small note ``<div>`` when the
    matrix is unavailable or Pillow is not installed -- never raises for those.
    """
    if D is None:
        return (
            '<div class="hm-note">Distance heatmap omitted: no distance matrix '
            "was provided (pass --matrix / dist_df).</div>"
        )
    try:
        from PIL import Image
    except ImportError:
        return (
            '<div class="hm-note">Distance heatmap omitted: the optional Pillow '
            "dependency is not installed (pip install &#39;pyeuk[report]&#39;).</div>"
        )
    D = np.asarray(D, dtype=float)
    n = D.shape[0]
    if n < 1:
        return '<div class="hm-note">Distance heatmap omitted: empty matrix.</div>'
    vmax = float(np.percentile(D[D > 0], 95)) if (D > 0).any() else 1.0
    img = Image.new("RGB", (n, n))
    px = img.load()
    for i in range(n):
        for j in range(n):
            v = min(D[i, j] / vmax, 1.0) if vmax > 0 else 0.0
            rr, gg, bb = colorsys.hls_to_rgb(196 / 360, (18 + v * 70) / 100, 0.40)
            px[j, i] = (int(rr * 255), int(gg * 255), int(bb * 255))
    cell = max(2, 480 // n)  # upscale for crisp blocks; NEAREST keeps cells sharp
    img = img.resize((n * cell, n * cell), Image.NEAREST)
    buf = io.BytesIO()
    img.save(buf, "PNG", optimize=True)
    b64 = base64.b64encode(buf.getvalue()).decode("ascii")
    return (
        f'<img src="data:image/png;base64,{b64}" width="100%" '
        f'style="max-width:{w}px;border-radius:3px;display:block" '
        f'alt="Tree-ordered pairwise distance heatmap; dark blocks are related specimens."/>'
    )


def cores_rows(r, top=8):
    """Return ((rank, size, members-preview) rows, total-core-count)."""
    cs = sorted(r.get("stable_cores", []), key=len, reverse=True)
    out = []
    for i, c in enumerate(cs[:top]):
        preview = ", ".join(_esc(m) for m in c[:6]) + (" &hellip;" if len(c) > 6 else "")
        out.append((i + 1, len(c), preview))
    return out, len(cs)


# ---------------------------------------------------------------------------
# Themes: palette + font stacks. `studio` links Google Fonts; `galaxy` has
# zero external assets and uses system fonts + the Galaxy brand palette.
# ---------------------------------------------------------------------------
_STUDIO_FONTS = (
    '<link rel="preconnect" href="https://fonts.googleapis.com">'
    '<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>'
    '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?'
    "family=Fraunces:opsz,wght@9..144,500;9..144,600&family=Inter:wght@400;500;600;700&"
    'family=IBM+Plex+Mono:wght@400;500;600&display=swap">'
)

_STUDIO_PALETTE = (
    ":root{--ground:#f4f6f7;--panel:#fff;--panel2:#eef2f4;--ink:#111a1d;--ink2:#46555b;"
    "--muted:#8a999e;--line:#d6dee1;--line2:#c4cfd3;--blue:#2a78d6;--good:#1a7f4b;"
    "--good-s:#dcefe4;--warn:#9a5b0c;--warn-s:#f6ead4;"
    "--mast:#16211f;--mast-fg:#f4f6f7;--mast-accent:#7fb8ee;"
    "--shadow:0 1px 2px rgba(17,26,29,.05),0 10px 26px -18px rgba(17,26,29,.4);"
    '--fd:"Fraunces",Georgia,serif;--fb:"Inter",-apple-system,Segoe UI,Helvetica,Arial,sans-serif;'
    '--fm:"IBM Plex Mono",ui-monospace,Menlo,monospace}'
    '@media(prefers-color-scheme:dark){:root:not([data-theme="light"]){--ground:#0b1214;'
    "--panel:#121b1e;--panel2:#18242a;--ink:#e7eef0;--ink2:#aab9bf;--muted:#6b7b81;"
    "--line:#213036;--line2:#2d4149;--blue:#3987e5;--good:#5fbe8e;--good-s:#12301f;"
    "--warn:#d3a34a;--warn-s:#32240f;"
    "--shadow:0 1px 2px rgba(0,0,0,.5),0 12px 30px -20px rgba(0,0,0,.9)}}"
    ':root[data-theme="dark"]{--ground:#0b1214;--panel:#121b1e;--panel2:#18242a;'
    "--ink:#e7eef0;--ink2:#aab9bf;--muted:#6b7b81;--line:#213036;--line2:#2d4149;"
    "--blue:#3987e5;--good:#5fbe8e;--good-s:#12301f;--warn:#d3a34a;--warn-s:#32240f;"
    "--shadow:0 1px 2px rgba(0,0,0,.5),0 12px 30px -20px rgba(0,0,0,.9)}"
)

_GALAXY_PALETTE = (
    ":root{--ground:#f8f9fa;--panel:#ffffff;--panel2:#eef1f4;--ink:#2c3143;--ink2:#4b5563;"
    "--muted:#868e96;--line:#dee2e6;--line2:#ced4da;--brand:#25537b;--blue:#2077b3;"
    "--good:#3f9b3f;--good-s:#e4f3e4;--warn:#c96a02;--warn-s:#ffedd6;--danger:#e31a1e;"
    "--mast:#2c3143;--mast-fg:#ffffff;--mast-accent:#8ec5ef;"
    "--shadow:0 1px 2px rgba(44,49,67,.06),0 8px 22px -16px rgba(44,49,67,.32);"
    '--fb:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Helvetica Neue",Arial,'
    '"Noto Sans","Liberation Sans",sans-serif;--fd:var(--fb);'
    '--fm:SFMono-Regular,Menlo,Monaco,Consolas,"Liberation Mono","Courier New",monospace}'
    '@media(prefers-color-scheme:dark){:root:not([data-theme="light"]){--ground:#1a1e2a;'
    "--panel:#232838;--panel2:#2c3143;--ink:#e9ecf1;--ink2:#aeb6c2;--muted:#7b8494;"
    "--line:#3a4155;--line2:#454d63;--brand:#5a90c0;--blue:#4a9fd6;--good:#66cc66;"
    "--good-s:#1c3320;--warn:#fe9a3a;--warn-s:#3a2a12;"
    "--shadow:0 1px 2px rgba(0,0,0,.5),0 12px 30px -20px rgba(0,0,0,.9)}}"
    ':root[data-theme="dark"]{--ground:#1a1e2a;--panel:#232838;--panel2:#2c3143;'
    "--ink:#e9ecf1;--ink2:#aeb6c2;--muted:#7b8494;--line:#3a4155;--line2:#454d63;"
    "--brand:#5a90c0;--blue:#4a9fd6;--good:#66cc66;--good-s:#1c3320;--warn:#fe9a3a;"
    "--warn-s:#3a2a12;--shadow:0 1px 2px rgba(0,0,0,.5),0 12px 30px -20px rgba(0,0,0,.9)}"
)

# Shared base rules (theme-independent). `--brand` falls back to `--ink` under the
# studio palette, which does not define it.
_BASE_CSS = (
    "*{box-sizing:border-box}"
    "body{margin:0;background:var(--ground);color:var(--ink);font-family:var(--fb);"
    "font-size:15px;line-height:1.5;-webkit-font-smoothing:antialiased}"
    ".ax{font-family:var(--fm);font-size:9px;fill:var(--muted)}"
    ".mono{font-family:var(--fm);font-variant-numeric:tabular-nums}"
    ".rk-headline{font-family:var(--fm);font-size:11.5px;color:var(--ink2);margin:0}"
    ".hm-note{font-size:12px;color:var(--muted);font-style:italic;padding:10px 0}"
    ".legend{display:flex;gap:14px;flex-wrap:wrap;font-size:11px;color:var(--ink2)}"
    ".lg{display:inline-flex;align-items:center;gap:6px}"
    ".lb{width:22px;border-top:2px solid var(--ink)}"
    ".lb.mid{border-top:1.3px solid var(--ink2);opacity:.72}"
    ".lb.fuz{border-top:1.3px dashed var(--muted);opacity:.5}"
    ".dot{width:9px;height:9px;border-radius:50%}"
    ".cbar{width:16px;height:5px;border-radius:2px;background:var(--blue);display:inline-block}"
    ".pyeuk-report.verdict-good{--accent:var(--good)}"
    ".pyeuk-report.verdict-warn{--accent:var(--warn)}"
)

_LEGEND = (
    '<div class="legend">'
    '<span class="lg"><span class="lb"></span>confident split</span>'
    '<span class="lg"><span class="lb mid"></span>moderate</span>'
    '<span class="lg"><span class="lb fuz"></span>uncertain</span>'
    '<span class="lg"><span class="cbar"></span>stable core</span>'
    '<span class="lg"><span class="dot" style="background:var(--blue)"></span>silhouette</span>'
    '<span class="lg"><span class="dot" style="background:var(--good)"></span>stability</span>'
    "</div>"
)


def _theme_head(theme):
    """Return (fonts-link, palette-css) for a theme name."""
    if theme == "galaxy":
        return "", _GALAXY_PALETTE
    if theme == "studio":
        return _STUDIO_FONTS, _STUDIO_PALETTE
    raise ValueError(f"unknown theme {theme!r}; expected one of {THEMES}")


# ---------------------------------------------------------------------------
# Shared data model for the flavors
# ---------------------------------------------------------------------------
def _pct(x):
    try:
        return int(round(float(x) * 100))
    except (TypeError, ValueError):
        return 0


def _provenance(r):
    lm = r.get("linkage_method", "ward")
    return f"{r.get('n', '?')} specimens &middot; wIBS + {_esc(lm).capitalize()}"


def _facts(r, D):
    """Assemble everything the flavor templates need from a sweep dict + matrix."""
    lo, hi = r.get("count_range", [None, None])
    sel = r.get("naive_selectors", {}) or {}
    conf = bool(r.get("confident"))
    point = r.get("point_estimate")
    cores, ncores = cores_rows(r, top=8)
    big = str(point) if (conf and point is not None) else f"{lo}&ndash;{hi}"
    pr = _pct(r.get("pairs_resolved", 0))
    if conf:
        blurb = (
            f"Independent selectors agree at {point} (knee {sel.get('knee', '?')}, "
            f"silhouette {sel.get('silhouette', '?')}, gap {sel.get('gap', '?')}) and the "
            f"tree has a dominant fully-supported split. {r.get('n_stable_cores', 0)} stable "
            f"cores; {pr}% of specimen pairs resolved."
        )
    else:
        spread_hi = max(sel.values()) if sel else "?"
        blurb = (
            f"Independent selectors scatter ({sel.get('knee', '?')}&ndash;{spread_hi}); the "
            f"tree resolves {r.get('count_at_solid_support', '?')} groups at full support and "
            f"{r.get('count_at_moderate_support', '?')} including moderate splits. "
            f"{r.get('n_stable_cores', 0)} stable cores; {pr}% of pairs resolved."
        )
    return dict(
        r=r,
        D=D,
        lo=lo,
        hi=hi,
        sel=sel,
        conf=conf,
        point=point,
        big=big,
        cls="good" if conf else "warn",
        vclass="verdict-good" if conf else "verdict-warn",
        head="the count is determined" if conf else "the count is not determined",
        blurb=blurb,
        cores=cores,
        ncores=ncores,
        prov=_provenance(r),
        pr=pr,
        headline=r.get("headline", ""),
    )


def _core_table(cores):
    rows = "".join(
        f'<tr><td class="mono">{i}</td><td class="mono">{n}</td><td>{m}</td></tr>'
        for i, n, m in cores
    )
    return (
        '<table class="cores"><thead><tr><th>core</th><th>n</th><th>members</th></tr>'
        f"</thead><tbody>{rows}</tbody></table>"
    )


def _document(title, theme, body_class, body, extra_css=""):
    """Wrap a flavor body in a complete, portable HTML document."""
    fonts, palette = _theme_head(theme)
    html = (
        "<!doctype html>\n"
        '<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f"<title>{_esc(title)}</title>\n"
        f"{fonts}"
        f"<style>{palette}{_BASE_CSS}{extra_css}</style>\n"
        "</head>\n"
        f'<body>\n<div class="pyeuk-report {body_class}">\n{body}\n</div>\n</body>\n</html>\n'
    )
    return _ascii_only(html)


# ---------------------------------------------------------------------------
# Flavor: DASHBOARD (default)
# ---------------------------------------------------------------------------
_DASHBOARD_CSS = (
    ".top{background:var(--mast);color:var(--mast-fg);padding:12px 22px;display:flex;"
    "justify-content:space-between;align-items:center;flex-wrap:wrap;gap:6px}"
    ".top .t{font-weight:600;font-size:15px;letter-spacing:.01em}"
    ".top .t b{color:var(--mast-accent);font-weight:600}"
    ".top .p{font-size:11.5px;color:rgba(255,255,255,.72);font-family:var(--fm)}"
    ".wrap{max-width:1060px;margin:0 auto;padding:18px 22px 50px;display:flex;"
    "flex-direction:column;gap:14px}"
    ".tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:1px;"
    "background:var(--line);border:1px solid var(--line);border-radius:4px;overflow:hidden}"
    ".tile{background:var(--panel);padding:14px 16px}"
    ".tn{font-family:var(--fm);font-size:25px;font-weight:600;font-variant-numeric:tabular-nums;"
    "line-height:1}"
    ".tk{font-size:12px;color:var(--muted);margin-top:5px}"
    ".grid{display:grid;grid-template-columns:1.3fr 1fr;gap:14px}"
    "@media(max-width:760px){.grid{grid-template-columns:1fr}}"
    ".card{background:var(--panel);border:1px solid var(--line);border-radius:4px;"
    "box-shadow:var(--shadow);padding:15px 17px;display:flex;flex-direction:column;gap:8px}"
    ".card h3{margin:0;font-size:13px;font-weight:600;color:var(--brand,var(--ink));"
    "display:flex;justify-content:space-between;align-items:center;gap:8px}"
    ".chip{font-size:10px;font-weight:600;padding:2px 9px;border-radius:3px;"
    "background:var(--warn-s);color:var(--warn);border:1px solid var(--warn)}"
    ".chip.g{background:var(--good-s);color:var(--good);border-color:var(--good)}"
    ".desc{font-size:12px;color:var(--muted)}"
    ".cores{width:100%;border-collapse:collapse;font-size:12px}"
    ".cores th{text-align:left;font-family:var(--fm);font-size:9.5px;text-transform:uppercase;"
    "letter-spacing:.06em;color:var(--muted);padding:3px 6px 3px 0;border-bottom:2px solid var(--line)}"
    ".cores td{padding:3px 6px 3px 0;border-bottom:1px solid var(--line)}"
    ".hmrow{display:grid;grid-template-columns:240px 1fr;gap:18px;align-items:center}"
    "@media(max-width:520px){.hmrow{grid-template-columns:1fr}}"
)


def _flavor_dashboard(f, theme):
    r = f["r"]
    tiles = [
        ("specimens", str(r.get("n", "?")), "ink"),
        ("count", f["big"], f["cls"]),
        ("determined?", ("yes" if f["conf"] else "no"), f["cls"]),
        ("stable cores", str(r.get("n_stable_cores", 0)), "good"),
        ("pairs resolved", f"{f['pr']}%", "good"),
        ("strong splits", str(r.get("strong_splits", 0)), "ink"),
    ]
    thtml = "".join(
        f'<div class="tile"><div class="tn" style="color:var(--{col})">{v}</div>'
        f'<div class="tk">{k}</div></div>'
        for k, v, col in tiles
    )
    chip = "count determined" if f["conf"] else "count not determined"
    chipc = "g" if f["conf"] else ""
    body = (
        f'<div class="top"><div class="t">PyEuk &middot; <b>run dashboard</b></div>'
        f'<div class="p">{f["prov"]}</div></div>'
        f'<div class="wrap">'
        f'<p class="rk-headline">{_esc(f["headline"])}</p>'
        f'<div class="tiles">{thtml}</div>'
        f'<div class="grid">'
        f'<div class="card"><h3>Confidence tree <span class="chip {chipc}">{chip}</span></h3>'
        f'{dendrogram_svg(r["tree"], w=560, h=246, cores=r.get("stable_cores"), corelabels=True)}'
        f"{_LEGEND}<div class=\"desc\">Solid branches reproduce across bootstraps; the numbered "
        f"blue bars beneath the leaves are the stable cores (numbers match the table).</div></div>"
        f'<div class="card"><h3>Count sweep</h3>{sweep_curve_svg(r, w=380, h=210)}'
        f'<div class="desc">Band = supported range. Silhouette and stability rising together '
        f"with no plateau means no single count is picked.</div></div>"
        f"</div>"
        f'<div class="grid">'
        f'<div class="card"><h3>Stable cores <span class="chip g">{f["ncores"]} reproduced</span>'
        f'</h3>{_core_table(f["cores"])}</div>'
        f'<div class="card"><h3>Distance structure</h3><div class="hmrow">'
        f'<div>{heatmap_img(f["D"], w=240)}</div>'
        f'<div class="desc">Pairwise wIBS distance, tree-ordered. Dark diagonal blocks are the '
        f"stable cores.</div></div></div>"
        f"</div>"
        f"</div>"
    )
    return _document("PyEuk Run Dashboard", theme, f["vclass"], body, _DASHBOARD_CSS)


# ---------------------------------------------------------------------------
# Flavor: CLINICAL
# ---------------------------------------------------------------------------
_CLINICAL_CSS = (
    ".wrap{max-width:820px;margin:0 auto;padding:44px 30px 70px}"
    ".masthead{display:flex;justify-content:space-between;align-items:flex-end;"
    "border-bottom:2.5px solid var(--ink);padding-bottom:12px;gap:12px;flex-wrap:wrap}"
    ".masthead h1{font-family:var(--fd);font-weight:600;font-size:25px;margin:0}"
    ".masthead .meta{font-family:var(--fm);font-size:10.5px;color:var(--muted);"
    "text-align:right;line-height:1.5}"
    ".prov{font-family:var(--fm);font-size:11px;color:var(--ink2);margin:8px 0 6px}"
    ".rk-headline{margin:0 0 20px}"
    ".verdict{border:1px solid var(--accent);background:var(--accent-s);border-radius:4px;"
    "padding:16px 20px;display:flex;gap:20px;align-items:center;margin-bottom:22px;flex-wrap:wrap}"
    ".verdict .big{font-family:var(--fm);font-size:44px;font-weight:600;line-height:1;"
    "font-variant-numeric:tabular-nums;white-space:nowrap;color:var(--accent)}"
    ".verdict .lab{font-size:13px;color:var(--ink2)}.verdict .lab b{color:var(--ink)}"
    ".sect{font-family:var(--fm);font-size:11px;letter-spacing:.14em;text-transform:uppercase;"
    "color:var(--blue);margin:24px 0 10px;border-bottom:1px solid var(--line);padding-bottom:5px}"
    ".two{display:grid;grid-template-columns:1fr 1fr;gap:22px;align-items:start}"
    "@media(max-width:640px){.two{grid-template-columns:1fr}}"
    ".cap{font-size:12px;color:var(--muted);margin-top:6px}.legend{margin-top:8px}"
    ".cores{width:100%;border-collapse:collapse;font-size:12.5px;margin-top:4px}"
    ".cores th{text-align:left;font-family:var(--fm);font-size:10px;text-transform:uppercase;"
    "letter-spacing:.08em;color:var(--muted);border-bottom:1px solid var(--line);padding:4px 8px 4px 0}"
    ".cores td{padding:4px 8px 4px 0;border-bottom:1px solid var(--line)}"
    ".foot{border-top:1px solid var(--line);margin-top:26px;padding-top:12px;font-size:11px;"
    "color:var(--muted);line-height:1.6}"
    # bind the confidence-tinted surface for the verdict box
    ".pyeuk-report.verdict-good{--accent-s:var(--good-s)}"
    ".pyeuk-report.verdict-warn{--accent-s:var(--warn-s)}"
)


def _flavor_clinical(f, theme):
    r = f["r"]
    band = f["big"]
    body = (
        f'<div class="wrap">'
        f'<div class="masthead"><h1>Clustering Report</h1>'
        f'<div class="meta">PyEuk sweep diagnostic<br>{_esc(r.get("linkage_method", "ward"))} linkage'
        f"<br>unsupervised</div></div>"
        f'<div class="prov">{f["prov"]}</div>'
        f'<p class="rk-headline">{_esc(f["headline"])}</p>'
        f'<div class="verdict"><div class="big">{f["big"]}</div>'
        f'<div class="lab"><b>clusters &mdash; {f["head"]}.</b> {f["blurb"]}</div></div>'
        f'<div class="sect">Confidence tree with stable cores</div>'
        f'{dendrogram_svg(r["tree"], w=760, h=210, cores=r.get("stable_cores"))}'
        f'<div class="cap">Ward tree; solid branches are splits the data reproduces, '
        f"faded/dashed are uncertain. Blue bars beneath the leaves are the {f['ncores']} stable "
        f"cores &mdash; groups reproduced in &ge;90% of bootstraps.</div>"
        f"{_LEGEND}"
        f'<div class="sect">Count sweep &amp; distance structure</div>'
        f'<div class="two">'
        f'<div>{sweep_curve_svg(r, w=380, h=180)}<div class="cap">Silhouette and bootstrap '
        f"stability across k; band marks {band}.</div></div>"
        f'<div>{heatmap_img(f["D"], w=240)}<div class="cap">Pairwise wIBS distance, tree-ordered '
        f"&mdash; dark blocks are related.</div></div>"
        f"</div>"
        f'<div class="sect">Stable cores</div>'
        f'{_core_table(f["cores"])}<div class="cap">{f["ncores"]} total; top {len(f["cores"])} shown.</div>'
        f'<div class="foot">Method: Ward linkage on the PyEuk wIBS distance matrix. Branch support '
        f"= 1 &minus; mean cross-cluster co-assignment over bootstrap resamples marginalised across "
        f"resolution. Count selectors: merge-gap knee, silhouette, Tibshirani gap statistic. "
        f'Unsupervised &mdash; no ground truth. Generated by <span class="mono">pyeuk report</span>.</div>'
        f"</div>"
    )
    return _document("PyEuk Clustering Report", theme, f["vclass"], body, _CLINICAL_CSS)


# ---------------------------------------------------------------------------
# Flavor: NARRATIVE
# ---------------------------------------------------------------------------
_NARRATIVE_CSS = (
    ".wrap{max-width:760px;margin:0 auto;padding:64px 28px 90px;display:flex;"
    "flex-direction:column;gap:30px}"
    ".eyebrow{font-family:var(--fm);font-size:11px;letter-spacing:.18em;text-transform:uppercase;"
    "color:var(--blue)}"
    ".hero h1{font-family:var(--fd);font-weight:600;font-size:clamp(30px,6vw,50px);line-height:1.04;"
    "margin:.3em 0 0;letter-spacing:-.02em;text-wrap:balance}"
    ".hero .sub{font-size:19px;color:var(--ink2);margin-top:16px;max-width:60ch}"
    ".hero .sub b{color:var(--ink)}"
    ".rk-headline{margin:12px 0 0}"
    ".herotree{background:var(--panel);border:1px solid var(--line);border-radius:8px;"
    "box-shadow:var(--shadow);padding:22px 18px 14px;margin-top:8px}"
    ".herotree .legend{margin-top:12px;justify-content:center;font-size:12px}"
    ".pull{font-family:var(--fd);font-weight:500;font-size:25px;line-height:1.25;color:var(--ink);"
    "border-left:3px solid var(--blue);padding-left:20px}"
    "h2{font-family:var(--fd);font-weight:600;font-size:24px;margin:14px 0 0}"
    "p{margin:12px 0;font-size:16.5px;color:var(--ink2)}p b{color:var(--ink)}"
    ".split{display:grid;grid-template-columns:1fr 1fr;gap:22px;align-items:center;"
    "background:var(--panel2);border-radius:8px;padding:18px}"
    "@media(max-width:600px){.split{grid-template-columns:1fr}}"
    ".num{font-family:var(--fm);font-size:40px;font-weight:600;color:var(--ink);"
    "font-variant-numeric:tabular-nums}"
    ".corechips{display:flex;flex-wrap:wrap;gap:6px;margin-top:10px}"
    ".cc{font-family:var(--fm);font-size:11px;background:var(--good-s);color:var(--good);"
    "border:1px solid var(--good);border-radius:12px;padding:3px 9px}"
    ".foot{border-top:1px solid var(--line);padding-top:16px;color:var(--muted);font-size:12.5px}"
)


def _flavor_narrative(f, theme):
    r = f["r"]
    sel = f["sel"]
    point = f["point"]
    if f["conf"]:
        h1 = (
            f"The data falls cleanly into {point} groups &mdash; and PyEuk says so with confidence."
        )
        sub = (
            "On this cohort the count is not a judgment call. Independent methods agree, one split "
            "dominates the tree, and the report says so plainly &mdash; a single number, earned."
        )
        why = (
            f"The merge-gap knee, silhouette, and gap statistic all land on about {point}, and the "
            f"tree has one dominant, fully-supported split with everything below it far shorter. "
            f"When the evidence converges like this, a single number is the honest answer."
        )
        report = f"<b>{point} clusters, with confidence.</b>"
        h2 = "Why a single number?"
    else:
        spread_hi = max(sel.values()) if sel else "?"
        h1 = (
            f"The data supports {f['lo']} to {f['hi']} groups &mdash; and says so, instead of "
            f"guessing one."
        )
        sub = (
            "This cohort does not fall into a single obvious number of clusters. PyEuk does not "
            "paper over that with one cut &mdash; it reports the <b>range the data supports</b>, "
            "the <b>structure it is sure of</b>, and draws its <b>own uncertainty</b> into the tree."
        )
        why = (
            f"Independent ways of choosing a count disagree: the knee says <b>{sel.get('knee', '?')}</b>, "
            f"silhouette <b>{sel.get('silhouette', '?')}</b>, the gap statistic <b>{sel.get('gap', '?')}</b>. "
            f"When honest methods scatter from {sel.get('knee', '?')} to {spread_hi}, the count is "
            f"genuinely undetermined."
        )
        report = (
            f"<b>{f['lo']}&ndash;{f['hi']} clusters, count not determined; "
            f"{r.get('n_stable_cores', 0)} stable cores.</b>"
        )
        h2 = "Why not a single number?"
    chips = "".join(f'<span class="cc">core {i} &middot; {n}</span>' for i, n, _ in f["cores"])
    body = (
        f'<div class="wrap">'
        f'<div class="hero"><div class="eyebrow">PyEuk &middot; {f["prov"]}</div>'
        f'<h1>{h1}</h1><p class="sub">{sub}</p>'
        f'<p class="rk-headline">{_esc(f["headline"])}</p></div>'
        f'<div class="herotree">{dendrogram_svg(r["tree"], w=680, h=250, cores=r.get("stable_cores"))}'
        f"{_LEGEND}</div>"
        f'<p class="pull">Solid branches are splits the data reproduces. Faded ones are splits it '
        f"doesn't. The blue bars are the groups you can act on.</p>"
        f"<h2>{h2}</h2><p>{why}</p>"
        f'<div class="split"><div><div class="num">{r.get("n_stable_cores", 0)}</div>'
        f'<div style="font-size:13px;color:var(--ink2)">stable cores &mdash; groups that cluster '
        f"together in &ge;90% of bootstrap resamples. The findings you can act on.</div>"
        f'<div class="corechips">{chips}</div></div>'
        f'<div><div class="num">{f["pr"]}%</div>'
        f'<div style="font-size:13px;color:var(--ink2)">of specimen pairs are resolved decisively. '
        f"The uncertainty is about <i>how many</i> groups, not <i>which</i> specimens are related.</div>"
        f"</div></div>"
        f"<h2>What to report</h2><p>{report} The same tool, reading the same tree, tells you which "
        f"kind of answer your data actually supports.</p>"
        f'<div class="foot">Ward linkage on the PyEuk wIBS distance matrix. Branch support = 1 '
        f"&minus; mean cross-cluster co-assignment over bootstrap resamples. Selectors: merge-gap "
        f"knee, silhouette, Tibshirani gap statistic. Unsupervised. "
        f'<span style="font-family:var(--fm)">pyeuk report --flavor narrative</span>.</div>'
        f"</div>"
    )
    return _document("What PyEuk Found", theme, f["vclass"], body, _NARRATIVE_CSS)


_FLAVOR_FUNCS = {
    "dashboard": _flavor_dashboard,
    "clinical": _flavor_clinical,
    "narrative": _flavor_narrative,
}


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------
def render(sweep, dist_df=None, flavor="dashboard", theme="studio"):
    """Render a self-contained HTML report from a cluster sweep.

    Parameters
    ----------
    sweep : dict or str
        A ``cluster_sweep()`` result dict, or a path to the ``*_SWEEP.json``
        it writes.
    dist_df : None, str, or DataFrame-like, optional
        Distance matrix for the tree-ordered heatmap. May be a path to a
        CSV/TSV, or a pandas DataFrame (ids on the index). If ``None`` -- or if
        Pillow is not installed -- the heatmap is skipped with a note.
    flavor : {"dashboard", "clinical", "narrative"}
        Report layout. ``dashboard`` is the default.
    theme : {"studio", "galaxy"}
        ``studio`` links Google Fonts (Fraunces/Inter); ``galaxy`` uses system
        fonts and the Galaxy brand palette with no external assets.

    Returns
    -------
    str
        A complete HTML document with zero raw non-ASCII bytes and no JavaScript.
    """
    if flavor not in _FLAVOR_FUNCS:
        raise ValueError(f"unknown flavor {flavor!r}; expected one of {FLAVORS}")
    if theme not in THEMES:
        raise ValueError(f"unknown theme {theme!r}; expected one of {THEMES}")
    r = _load_sweep(sweep)
    tree = r.get("tree", {})
    D = _load_matrix(dist_df, tree.get("leaf_order", []))
    f = _facts(r, D)
    return _FLAVOR_FUNCS[flavor](f, theme)
