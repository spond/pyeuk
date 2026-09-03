"""
Unit tests for the cluster-count SWEEP diagnostic (CyclosporaClusterFinder.cluster_sweep).

The sweep is the default clustering answer: rather than commit to one k, it reports the range
of counts the data supports, whether that count is determined (do the independent selectors
agree?), and a per-branch confidence tree. These tests pin the two behaviours that matter --
a well-separated cohort yields a confident single number, a diffuse one yields a range with
its confident sub-structure -- plus determinism and the emitted artifacts.
"""

import glob
import os
import tempfile
import unittest

import numpy as np
import pandas as pd

from pyeuk.clustering import CyclosporaClusterFinder


def _separated_blocks(k=4, per=12, seed=1):
    """k tight, WELL-separated blocks: within ~0.02, between ~0.6. Every count selector
    should land on k, so the sweep must report k with confidence."""
    rng = np.random.default_rng(seed)
    labels = [f"C{b}" for b in range(k) for _ in range(per)]
    n = len(labels)
    D = np.full((n, n), 0.6)
    for i in range(n):
        for j in range(n):
            if i == j:
                D[i, j] = 0.0
            elif labels[i] == labels[j]:
                D[i, j] = abs(rng.normal(0.02, 0.004))
    D = (D + D.T) / 2.0
    np.fill_diagonal(D, 0.0)
    ids = [f"{labels[i]}_{i}" for i in range(n)]
    return pd.DataFrame(D, index=ids, columns=ids), k


def _diffuse_cloud(n=44, seed=3):
    """Points spread uniformly with no dominant grouping: the count is genuinely
    undetermined, so the sweep must decline a single number and report a range."""
    rng = np.random.default_rng(seed)
    X = rng.uniform(0, 1, (n, 2))
    D = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
    D = (D + D.T) / 2.0
    np.fill_diagonal(D, 0.0)
    ids = [f"s{i}" for i in range(n)]
    return pd.DataFrame(D, index=ids, columns=ids)


class TestClusterSweep(unittest.TestCase):

    def test_confident_on_separated_blocks(self):
        df, k = _separated_blocks(4, 12)
        r = CyclosporaClusterFinder().cluster_sweep(df, n_boot=40, seed=0)
        self.assertTrue(r["confident"])
        self.assertEqual(r["point_estimate"], k)
        self.assertEqual(r["count_range"], [k, k])

    def test_confidence_tree_well_formed(self):
        df, _ = _separated_blocks(3, 10)
        r = CyclosporaClusterFinder().cluster_sweep(df, n_boot=30, seed=0)
        nodes = r["tree"]["nodes"]
        self.assertEqual(len(nodes), df.shape[0] - 1)          # n-1 internal nodes
        self.assertEqual(len(r["tree"]["leaf_order"]), df.shape[0])
        for nd in nodes:
            self.assertGreaterEqual(nd["support"], 0.0)
            self.assertLessEqual(nd["support"], 1.0)
        self.assertTrue(r["tree"]["newick"].strip().endswith(";"))

    def test_deterministic_for_fixed_seed(self):
        df, _ = _separated_blocks(3, 10)
        f = CyclosporaClusterFinder()
        r1 = f.cluster_sweep(df, n_boot=30, seed=0)
        r2 = f.cluster_sweep(df, n_boot=30, seed=0)
        self.assertEqual(r1["count_range"], r2["count_range"])
        self.assertEqual([n["support"] for n in r1["tree"]["nodes"]],
                         [n["support"] for n in r2["tree"]["nodes"]])

    def test_fuzzy_cohort_reports_range_not_a_number(self):
        df = _diffuse_cloud()
        r = CyclosporaClusterFinder().cluster_sweep(df, n_boot=40, seed=0, do_gap=False)
        self.assertFalse(r["confident"])
        self.assertIsNone(r["point_estimate"])
        lo, hi = r["count_range"]
        self.assertLessEqual(lo, hi)
        self.assertIn("stable_cores", r)

    def test_writes_sweep_artifacts(self):
        df, _ = _separated_blocks(3, 8)
        with tempfile.TemporaryDirectory() as d:
            CyclosporaClusterFinder().cluster_sweep(df, n_boot=20, seed=0, do_gap=False, output_dir=d)
            self.assertTrue(glob.glob(os.path.join(d, "*_SWEEP.json")))
            self.assertTrue(glob.glob(os.path.join(d, "*_confidence_tree.nwk")))
            self.assertTrue(glob.glob(os.path.join(d, "*_RESULTING_CLUSTERS_*.txt")))


if __name__ == "__main__":
    unittest.main()
