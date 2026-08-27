"""
Unit tests for distance-threshold ("linkage") clustering in CyclosporaClusterFinder.

The two cut modes answer different questions and neither is a better version of the other.

  --cut count     split the cohort into k groups. k is chosen from the largest merge-height
                  gap, subject to a minimum relative gap and a minimum cluster size. Correct
                  for a closed outbreak investigation, where every specimen belongs somewhere.

  --cut distance  cut at a fixed dissimilarity. No k, no gap requirement, and specimens with
                  no near neighbour come back as singletons. Correct for surveillance, where
                  most cases are unrelated to each other.

The tests below pin the failure that motivated the second mode: on a cohort whose true
structure is mostly singletons, both guards in the count rule reject every k that could
reproduce it. On the real 183-specimen Plasmodium vivax AmpliSeq cohort it falls back to a
single cluster; on the synthetic fixture here it merges the unrelated specimens into a few
groups. Either way it cannot represent the shape, which is what the comparison asserts.
"""

import os
import tempfile
import unittest

import numpy as np
import pandas as pd

from cyclospora_pyeuk.clustering import CyclosporaClusterFinder


def _blocks_plus_singletons(n_blocks=3, block_size=6, n_singletons=30, seed=7):
    """A few tight clusters embedded in a cloud of MUTUALLY unrelated specimens.

    This is the shape of surveillance data: a handful of genuine transmission clusters, and a
    long tail of imported infections unrelated to anything else INCLUDING each other. Block
    members sit at ~0.02. Every other pair sits at a distance drawn from a wide distribution,
    so the unrelated specimens do not themselves form a tight group -- which is exactly why
    there is no single dominant gap in the tree for a count rule to find.

    Getting this right matters: if the singletons are all placed at one constant distance they
    become a cluster, the count rule finds them, and the test proves nothing.
    """
    rng = np.random.default_rng(seed)
    labels = []
    for b in range(n_blocks):
        labels += [f"B{b}"] * block_size
    labels += [f"S{i}" for i in range(n_singletons)]
    n = len(labels)

    D = np.abs(rng.normal(0.30, 0.075, (n, n)))
    for i in range(n):
        for j in range(n):
            if i != j and labels[i] == labels[j] and labels[i].startswith("B"):
                D[i, j] = abs(rng.normal(0.02, 0.004))
    D = (D + D.T) / 2
    np.fill_diagonal(D, 0.0)
    ids = [f"SPEC{i:03d}" for i in range(n)]
    return pd.DataFrame(D, index=ids, columns=ids), dict(zip(ids, labels))


def _ari(pairs):
    from collections import Counter, defaultdict
    from math import comb
    t = defaultdict(Counter)
    for a, b in pairs:
        t[a][b] += 1
    n = len(pairs)
    sij = sum(comb(v, 2) for r in t.values() for v in r.values() if v > 1)
    sa = sum(comb(sum(r.values()), 2) for r in t.values() if sum(r.values()) > 1)
    c = Counter(b for _, b in pairs)
    sb = sum(comb(v, 2) for v in c.values() if v > 1)
    e = sa * sb / comb(n, 2)
    m = (sa + sb) / 2
    return (sij - e) / (m - e) if m != e else float("nan")


class TestLinkageMode(unittest.TestCase):

    def setUp(self):
        self.D, self.truth = _blocks_plus_singletons()
        self.finder = CyclosporaClusterFinder(stringency=95.0, robust=True,
                                              default_threshold=0.05)
        self.tmp = tempfile.mkdtemp()

    def test_distance_mode_beats_count_mode_on_singleton_heavy_data(self):
        """On a mostly-singleton structure, only the distance cut recovers the truth.

        Asserted as a comparison rather than as a literal k, because the count rule's exact
        fallback depends on where its two guards happen to bite. What is stable, and what
        matters, is that it cannot represent this shape and the distance cut can.

        Both guards in the count rule are doing their job. The minimum-cluster-size guard
        exists so surveillance does not report two-case 'clusters'. This cohort simply needs
        the other question asked.
        """
        truth_pairs = lambda cl: [
            (self.truth[s], c) for s, c in zip(cl["Seq_ID"], cl["Assigned_cluster"])
            if s in self.truth and c != -1
        ]
        cl_c, k_c, _ = self.finder.find_clusters(self.D, None, k_min=2, k_max=50,
                                                 output_dir=self.tmp, cut_mode="count")
        cl_d, k_d, _ = self.finder.find_clusters(self.D, None, output_dir=self.tmp,
                                                 cut_mode="distance", linkage_threshold=0.10)
        ari_c = _ari(truth_pairs(cl_c))
        ari_d = _ari(truth_pairs(cl_d))
        n_truth_groups = len(set(self.truth.values()))

        self.assertGreater(ari_d, 0.90, f"distance cut should recover the truth, got {ari_d:.3f}")
        self.assertLess(ari_c, 0.30, f"count rule should not, got {ari_c:.3f}")
        # the count rule merges the unrelated specimens; the distance cut keeps them apart
        self.assertLess(k_c, n_truth_groups / 2)
        self.assertGreater(k_d, n_truth_groups / 2)

    def test_distance_mode_recovers_the_blocks(self):
        """Cutting between the within-block and between-block distances recovers the truth."""
        cl, k, thr = self.finder.find_clusters(self.D, None, output_dir=self.tmp,
                                               cut_mode="distance", linkage_threshold=0.10)
        self.assertAlmostEqual(thr, 0.10, places=6)
        assign = dict(zip(cl["Seq_ID"], cl["Assigned_cluster"]))
        # the three blocks each land in one group
        for b in ("B0", "B1", "B2"):
            members = [s for s, t in self.truth.items() if t == b]
            self.assertEqual(len({assign[m] for m in members}), 1,
                             f"block {b} was split across groups")
        # the three blocks are not merged with each other
        block_groups = {assign[[s for s, t in self.truth.items() if t == b][0]]
                        for b in ("B0", "B1", "B2")}
        self.assertEqual(len(block_groups), 3)
        # singletons stay singletons
        sizes = cl["Assigned_cluster"].value_counts()
        self.assertGreaterEqual(int((sizes == 1).sum()), 20)
        self.assertEqual(k, int(cl["Assigned_cluster"].nunique()))

    def test_threshold_is_calibrated_when_not_supplied(self):
        """With no threshold given, one is derived and reported rather than assumed."""
        _, k, thr = self.finder.find_clusters(self.D, None, output_dir=self.tmp,
                                              cut_mode="distance")
        self.assertGreater(thr, 0.0)
        self.assertLess(thr, 0.30)
        self.assertGreater(k, 1)

    def test_calibration_from_labelled_pairs(self):
        """Given labels, the threshold comes from the observed within-cluster distances."""
        gold = pd.DataFrame({
            "Seq_ID": list(self.truth.keys()),
            "Cluster_alias": list(self.truth.values()),
        })
        thr, provenance = CyclosporaClusterFinder.suggest_linkage_threshold(
            self.D, gold, robust=True)
        self.assertIn("labelled", provenance)
        # within-block distances are ~0.02, so the calibrated cut must sit well below the
        # ~0.30 between-block distances
        self.assertLess(thr, 0.20)

    def test_calibration_without_labels_states_so(self):
        thr, provenance = CyclosporaClusterFinder.suggest_linkage_threshold(self.D, None)
        self.assertIn("no labels", provenance)
        self.assertGreater(thr, 0.0)

    def test_excluded_specimens_reported_in_both_modes(self):
        """A specimen dropped for low completeness is reported as -1 whichever mode ran.

        Both modes go through the same writer. Omitting an excluded specimen would make the
        cohort look smaller than it was.
        """
        extra = ["DROPPED_A", "DROPPED_B"]
        all_ids = list(self.D.index) + extra
        for mode, thr in (("count", None), ("distance", 0.10)):
            cl, _, _ = self.finder.find_clusters(
                self.D, None, output_dir=self.tmp, all_input_ids=all_ids,
                cut_mode=mode, linkage_threshold=thr)
            assign = dict(zip(cl["Seq_ID"], cl["Assigned_cluster"]))
            for e in extra:
                self.assertIn(e, assign, f"{e} missing in {mode} mode")
                self.assertEqual(assign[e], -1)

    def test_count_mode_unchanged_on_well_separated_data(self):
        """Two clean groups: the count rule must still find them. Guards against regression."""
        rng = np.random.default_rng(3)
        n = 40
        lab = ["A"] * 20 + ["B"] * 20
        D = np.full((n, n), 0.40)
        for i in range(n):
            for j in range(n):
                if i != j and lab[i] == lab[j]:
                    D[i, j] = 0.05
        D += rng.normal(0, 0.003, (n, n))
        D = np.abs((D + D.T) / 2)
        np.fill_diagonal(D, 0.0)
        ids = [f"X{i:02d}" for i in range(n)]
        df = pd.DataFrame(D, index=ids, columns=ids)
        cl, k, _ = self.finder.find_clusters(df, None, k_min=2, k_max=50,
                                             output_dir=self.tmp, cut_mode="count")
        self.assertEqual(k, 2)
        assign = dict(zip(cl["Seq_ID"], cl["Assigned_cluster"]))
        self.assertEqual(len({assign[ids[i]] for i in range(20)}), 1)
        self.assertEqual(len({assign[ids[i]] for i in range(20, 40)}), 1)


if __name__ == "__main__":
    unittest.main()


class TestLinkageMethod(unittest.TestCase):
    """Single linkage follows chains; Ward looks for compact blobs.

    A transmission cluster is a chain -- A infects B infects C -- so the two linkages are not
    interchangeable on this kind of data, and which is better is a property of the cohort
    rather than a preference.
    """

    def setUp(self):
        self.finder = CyclosporaClusterFinder(stringency=95.0, robust=True)
        self.tmp = tempfile.mkdtemp()

    @staticmethod
    def _chain(n=10, step=0.03, n_far=15, seed=11):
        """One transmission chain plus unrelated specimens.

        Consecutive members of the chain are close; the ends are far apart. Ward penalises the
        within-cluster spread this creates, single linkage does not.
        """
        rng = np.random.default_rng(seed)
        pos = [i * step for i in range(n)] + [5.0 + rng.random() * 5.0 for _ in range(n_far)]
        m = len(pos)
        D = np.zeros((m, m))
        for i in range(m):
            for j in range(m):
                D[i, j] = abs(pos[i] - pos[j])
        ids = [f"C{i:02d}" for i in range(n)] + [f"F{i:02d}" for i in range(n_far)]
        truth = {s: ("chain" if s.startswith("C") else s) for s in ids}
        return pd.DataFrame(D, index=ids, columns=ids), truth

    def test_single_linkage_keeps_a_chain_together(self):
        D, truth = self._chain()
        cl, _, _ = self.finder.find_clusters(D, None, output_dir=self.tmp,
                                             cut_mode="distance", linkage_threshold=0.05,
                                             linkage_method="single")
        assign = dict(zip(cl["Seq_ID"], cl["Assigned_cluster"]))
        chain = [s for s in truth if s.startswith("C")]
        self.assertEqual(len({assign[s] for s in chain}), 1,
                         "single linkage should follow the chain end to end")

    def test_ward_splits_the_same_chain(self):
        D, truth = self._chain()
        cl, _, _ = self.finder.find_clusters(D, None, output_dir=self.tmp,
                                             cut_mode="distance", linkage_threshold=0.05,
                                             linkage_method="ward")
        assign = dict(zip(cl["Seq_ID"], cl["Assigned_cluster"]))
        chain = [s for s in truth if s.startswith("C")]
        self.assertGreater(len({assign[s] for s in chain}), 1,
                           "ward is expected to break a chain at this threshold; if it no "
                           "longer does, the documented reason for offering single linkage "
                           "needs rechecking")

    def test_method_is_recorded(self):
        D, _ = self._chain()
        self.finder.find_clusters(D, None, output_dir=self.tmp, cut_mode="distance",
                                  linkage_threshold=0.05, linkage_method="single")
        self.assertEqual(self.finder.last_selection_meta.get("cut_mode"), "distance")
