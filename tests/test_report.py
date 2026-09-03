"""
Unit tests for the graphical report renderer (``pyeuk.report``).

The report turns a ``cluster_sweep()`` result into a single self-contained HTML
document -- confident cohorts as a green single number, fuzzy ones as an amber
range, with the stable cores drawn on the tree. These tests pin the guarantees
that matter for portability and for embedding inside Galaxy: every flavor renders
for both a confident and a fuzzy cohort, the verdict colour class tracks
confidence, the bytes are pure ASCII (HTML entities, never raw non-ASCII), and the
``galaxy`` theme carries no external assets.

``pyeuk.report`` is loaded standalone via importlib rather than ``import pyeuk`` so
that the suite runs in a lean environment (numpy [+ optional Pillow] only) without
dragging in numba or the heavy package siblings through the package ``__init__``.
"""

import importlib.util
import os
import unittest

_REPORT_PATH = os.path.join(os.path.dirname(__file__), "..", "pyeuk", "report.py")
_spec = importlib.util.spec_from_file_location("pyeuk_report_standalone", _REPORT_PATH)
report = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(report)


def _tree(n):
    """A small but well-formed confidence tree: n leaves, n-1 internal nodes,
    mixing high- and low-support branches so both branch styles render."""
    leaves = [f"S{i:02d}" for i in range(n)]
    nodes = []
    for i in range(n - 1):
        nodes.append(
            {
                "ya": float(i),
                "yb": float(i + 1),
                "h": round(0.12 * (i + 1), 3),
                "hca": 0.0,
                "hcb": 0.0,
                # alternate strong / weak so solid and dashed branches both appear
                "support": 0.95 if i % 2 == 0 else 0.30,
            }
        )
    return {
        "leaf_order": leaves,
        "hmax": round(0.12 * (n - 1), 3),
        "nodes": nodes,
        "newick": "(" + ",".join(leaves) + ");",
    }


def _sweep_rows(k_hi=8):
    return [
        {"k": k, "clusters": k, "silhouette": round(0.9 - 0.05 * k, 3),
         "rel_gap": 0.2, "singletons": 0, "stability": round(0.95 - 0.03 * k, 3)}
        for k in range(2, k_hi + 1)
    ]


def _confident_sweep():
    tree = _tree(12)
    cores = [tree["leaf_order"][:7], tree["leaf_order"][7:11]]  # one >6 to test the ellipsis
    return {
        "n": 12,
        "linkage_method": "ward",
        "count_range": [3, 3],
        "point_estimate": 3,
        "confident": True,
        "count_at_solid_support": 3,
        "count_at_moderate_support": 3,
        "naive_selectors": {"knee": 3, "silhouette": 3, "stability": 3, "gap": 3},
        "pairs_resolved": 0.91,
        "strong_splits": 3,
        "n_stable_cores": len(cores),
        "stable_cores": cores,
        "sweep": _sweep_rows(),
        "tree": tree,
        "headline": "3 clusters (confident: count selectors agree)",
    }


def _fuzzy_sweep():
    tree = _tree(14)
    cores = [tree["leaf_order"][:4], tree["leaf_order"][6:9]]
    return {
        "n": 14,
        "linkage_method": "ward",
        "count_range": [2, 5],
        "point_estimate": None,
        "confident": False,
        "count_at_solid_support": 2,
        "count_at_moderate_support": 5,
        "naive_selectors": {"knee": 2, "silhouette": 6, "stability": 6, "gap": 8},
        "pairs_resolved": 0.74,
        "strong_splits": 2,
        "n_stable_cores": len(cores),
        "stable_cores": cores,
        "sweep": _sweep_rows(),
        "tree": tree,
        "headline": "2-5 clusters (count not determined; selectors scatter 2-8)",
    }


def _identity_matrix_df(sweep):
    """A trivial DataFrame-like distance matrix over the leaves (for the heatmap
    path). Kept dependency-light: a tiny shim exposing .index and .values, so the
    test needs neither pandas nor a temp file."""
    import numpy as np

    ids = list(sweep["tree"]["leaf_order"])
    n = len(ids)
    M = np.abs(np.subtract.outer(np.arange(n), np.arange(n))) / float(n)

    class _DF:
        def __init__(self, index, values):
            self.index = index
            self.values = values

    return _DF(ids, M)


class TestReportRender(unittest.TestCase):

    def test_all_flavors_render_confident_and_fuzzy(self):
        for sweep in (_confident_sweep(), _fuzzy_sweep()):
            for flavor in report.FLAVORS:
                for theme in report.THEMES:
                    html = report.render(sweep, flavor=flavor, theme=theme)
                    self.assertTrue(html, f"empty html for {flavor}/{theme}")
                    self.assertGreater(len(html), 2000)
                    self.assertIn("<!doctype html>", html)
                    self.assertIn('charset="utf-8"', html)
                    # the canonical sweep headline is embedded verbatim
                    self.assertIn(sweep["headline"], html)

    def test_verdict_colour_class_tracks_confidence(self):
        for flavor in report.FLAVORS:
            conf = report.render(_confident_sweep(), flavor=flavor)
            fuzz = report.render(_fuzzy_sweep(), flavor=flavor)
            # confident -> green verdict; fuzzy -> amber verdict
            self.assertIn("verdict-good", conf)
            self.assertNotIn("pyeuk-report verdict-warn", conf)
            self.assertIn("verdict-warn", fuzz)
            self.assertNotIn("pyeuk-report verdict-good", fuzz)

    def test_confident_shows_single_number_fuzzy_shows_range(self):
        conf = report.render(_confident_sweep(), flavor="dashboard")
        fuzz = report.render(_fuzzy_sweep(), flavor="dashboard")
        self.assertIn("count not determined", fuzz)
        self.assertIn("&ndash;", fuzz)          # the amber range uses an en-dash entity
        self.assertIn("count determined", conf)

    def test_zero_raw_non_ascii_bytes(self):
        for sweep in (_confident_sweep(), _fuzzy_sweep()):
            for flavor in report.FLAVORS:
                for theme in report.THEMES:
                    html = report.render(sweep, flavor=flavor, theme=theme)
                    raw = html.encode("utf-8")
                    self.assertFalse(
                        any(b > 127 for b in raw),
                        f"{flavor}/{theme} contains raw non-ASCII bytes",
                    )
                    # a stronger equivalent guarantee: the string is pure ASCII
                    html.encode("ascii")

    def test_galaxy_theme_has_no_external_assets(self):
        for sweep in (_confident_sweep(), _fuzzy_sweep()):
            for flavor in report.FLAVORS:
                html = report.render(sweep, flavor=flavor, theme="galaxy")
                self.assertNotIn("http://", html)
                self.assertNotIn("https://", html)
                self.assertNotIn("fonts.googleapis.com", html)

    def test_studio_theme_links_google_fonts(self):
        html = report.render(_confident_sweep(), flavor="dashboard", theme="studio")
        self.assertIn("https://fonts.googleapis.com", html)

    def test_heatmap_embeds_png_when_matrix_given(self):
        sweep = _confident_sweep()
        df = _identity_matrix_df(sweep)
        html = report.render(sweep, dist_df=df, flavor="dashboard")
        try:
            import PIL  # noqa: F401
        except ImportError:
            self.assertIn("Pillow", html)  # graceful note, no crash
        else:
            self.assertIn("data:image/png;base64,", html)

    def test_missing_matrix_is_skipped_gracefully(self):
        html = report.render(_confident_sweep(), dist_df=None, flavor="dashboard")
        self.assertIn("no distance matrix", html)
        self.assertNotIn("Traceback", html)

    def test_accepts_json_path_and_bad_args_raise(self):
        import json
        import tempfile

        sweep = _fuzzy_sweep()
        with tempfile.NamedTemporaryFile("w", suffix="_SWEEP.json", delete=False) as fh:
            json.dump(sweep, fh)
            path = fh.name
        try:
            html = report.render(path, flavor="narrative", theme="galaxy")
            self.assertIn(sweep["headline"], html)
        finally:
            os.unlink(path)
        with self.assertRaises(ValueError):
            report.render(sweep, flavor="nope")
        with self.assertRaises(ValueError):
            report.render(sweep, theme="nope")


if __name__ == "__main__":
    unittest.main()
