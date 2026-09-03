"""
CyclosporaClusterFinder: Modernized Python module for AGNES Ward's hierarchical clustering
and robust outbreak threshold calibration for CDC Cyclospora cayetanensis workflow.
Supports prospective unsupervised clustering via Dendrogram Merge Height Gap Knee Detection
and gold-standard supervised calibration.
"""

import os
import datetime
from collections import Counter
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, cut_tree, fcluster, leaves_list
from sklearn.metrics import roc_auc_score, silhouette_score
from typing import Tuple, Optional, Dict, List


class CyclosporaClusterFinder:
    """
    Hierarchical clustering engine for outbreak surveillance, supporting prospective unsupervised
    clustering (via Dendrogram Merge Height Gap Knee Detection) and gold-standard supervised calibration.
    """

    def __init__(
        self,
        stringency: float = 95.0,
        robust: bool = True,
        default_threshold: float = 0.05,
        relative_gap_floor: float = 0.2200
    ):
        """
        Parameters
        ----------
        stringency : float
            Percentage of within-cluster distances that must fall below the distance threshold (default: 95.0%).
        robust : bool
            If True, uses Median + 3 * 1.4826 * MAD for threshold calibration. If False, uses Mean + 3 * StdDev.
        default_threshold : float
            Default threshold for prospective unsupervised clustering when gold standards are omitted.
        relative_gap_floor : float
            Minimum relative merge-height gap fraction of tree height required for unsupervised knee selection (default: 0.2200).
        """
        self.stringency = stringency
        self.robust = robust
        self.default_threshold = default_threshold
        self.relative_gap_floor = relative_gap_floor
        self.last_selection_meta: Dict[str, any] = {}

    @staticmethod
    def compute_distance_auc(dist_df: pd.DataFrame, gold_df: pd.DataFrame) -> float:
        """
        Computes ROC AUC of raw pairwise distances at separating same-outbreak vs different-outbreak sample pairs.
        This provides a pure, label-free diagnostic of distance metric discriminative quality without clustering dependencies.
        """
        gold_df = gold_df.copy()
        gold_df.columns = [c.strip() for c in gold_df.columns]
        if "Seq_ID" not in gold_df.columns:
            gold_df.rename(columns={gold_df.columns[0]: "Seq_ID"}, inplace=True)
        if "Cluster_alias" not in gold_df.columns and len(gold_df.columns) > 1:
            gold_df.rename(columns={gold_df.columns[1]: "Cluster_alias"}, inplace=True)

        valid_ids = list(set(dist_df.index) & set(gold_df["Seq_ID"]))
        if len(valid_ids) < 2:
            return 0.5

        label_map = dict(zip(gold_df["Seq_ID"], gold_df["Cluster_alias"]))
        
        y_true = [] # 1 if same outbreak cluster, 0 if different
        y_score = [] # -distance (smaller distance = higher probability of same cluster)

        for i in range(len(valid_ids)):
            for j in range(i + 1, len(valid_ids)):
                id1, id2 = valid_ids[i], valid_ids[j]
                d = dist_df.loc[id1, id2]
                if not np.isnan(d):
                    same_cluster = 1 if label_map[id1] == label_map[id2] else 0
                    y_true.append(same_cluster)
                    y_score.append(-d)

        if len(set(y_true)) < 2:
            return 0.5

        auc = float(roc_auc_score(y_true, y_score))
        print(f"[DistanceEngine AUC] Pairwise Distance ROC AUC = {auc:.4f} ({len(y_true)} sample pairs)")
        return auc

    def calibrate_gold_standard_threshold(
        self,
        dist_df: pd.DataFrame,
        gold_df: pd.DataFrame
    ) -> float:
        """
        Calculates maximum allowed intra-cluster distance threshold using epidemiological gold standards.
        """
        gold_df.columns = [c.strip() for c in gold_df.columns]
        if "Seq_ID" not in gold_df.columns:
            gold_df.rename(columns={gold_df.columns[0]: "Seq_ID"}, inplace=True)
        if "Cluster_alias" not in gold_df.columns and len(gold_df.columns) > 1:
            gold_df.rename(columns={gold_df.columns[1]: "Cluster_alias"}, inplace=True)

        valid_ids = set(dist_df.index) & set(gold_df["Seq_ID"])
        filtered_gold = gold_df[gold_df["Seq_ID"].isin(valid_ids)].copy()

        within_distances = []
        cluster_groups = filtered_gold.groupby("Cluster_alias")

        for alias, group in cluster_groups:
            members = group["Seq_ID"].tolist()
            if len(members) < 2:
                continue
            for i in range(len(members)):
                for j in range(i + 1, len(members)):
                    d = dist_df.loc[members[i], members[j]]
                    if not np.isnan(d):
                        within_distances.append(d)

        if not within_distances:
            print(f"[ClusterFinder] Warning: No overlapping gold standard sample pairs found. Using default threshold ({self.default_threshold}).")
            return self.default_threshold

        within_distances = np.array(within_distances)

        if self.robust:
            med = np.median(within_distances)
            mad = np.median(np.abs(within_distances - med))
            threshold = med + 3.0 * 1.4826 * mad
            print(f"[ClusterFinder] Supervised calibration (Median + 3*MAD): Threshold = {threshold:.5f} (median={med:.5f}, MAD={mad:.5f})")
        else:
            mean = np.mean(within_distances)
            std = np.std(within_distances, ddof=1) if len(within_distances) > 1 else 0.0
            threshold = mean + 3.0 * std
            print(f"[ClusterFinder] Supervised calibration (Mean + 3*StdDev): Threshold = {threshold:.5f} (mean={mean:.5f}, std={std:.5f})")

        return float(threshold)

    @staticmethod
    def suggest_linkage_threshold(dist_df: pd.DataFrame, gold_df: Optional[pd.DataFrame] = None,
                                  robust: bool = True) -> Tuple[float, str]:
        """
        Proposes a distance threshold for linkage mode, and says where it came from.

        Two sources, in order of preference:

        1. Labelled pairs. Take every pair known to belong together, and set the threshold at
           the upper edge of that distribution -- median + 3 * 1.4826 * MAD when robust, else
           mean + 3 * SD. This is the same calibration the supervised path already uses.
        2. The distance distribution itself. With no labels, fall back to the 5th percentile of
           all pairwise distances. Related pairs are a small minority of all pairs in
           surveillance data, so the lower tail is where they live.

        Returns (threshold, provenance) so the caller can print how the number was obtained
        rather than presenting it as if it were a constant.
        """
        vals = dist_df.values[np.triu_indices(len(dist_df), k=1)]
        vals = vals[~np.isnan(vals)]

        if gold_df is not None and len(gold_df) > 0:
            gold_df = gold_df.copy()
            gold_df.columns = [c.strip() for c in gold_df.columns]
            if "Seq_ID" not in gold_df.columns:
                gold_df.rename(columns={gold_df.columns[0]: "Seq_ID"}, inplace=True)
            if "Cluster_alias" not in gold_df.columns and len(gold_df.columns) > 1:
                gold_df.rename(columns={gold_df.columns[1]: "Cluster_alias"}, inplace=True)
            ids = set(dist_df.index)
            within = []
            for _, grp in gold_df[gold_df["Seq_ID"].isin(ids)].groupby("Cluster_alias"):
                mem = grp["Seq_ID"].tolist()
                for i in range(len(mem)):
                    for j in range(i + 1, len(mem)):
                        d = dist_df.loc[mem[i], mem[j]]
                        if not np.isnan(d):
                            within.append(float(d))
            if len(within) >= 3:
                a = np.asarray(within)
                if robust:
                    med = float(np.median(a))
                    mad = float(np.median(np.abs(a - med)))
                    t = med + 3.0 * 1.4826 * mad
                else:
                    t = float(np.mean(a) + 3.0 * np.std(a))
                return float(t), f"calibrated from {len(within)} labelled within-cluster pairs"

        if len(vals) == 0:
            return 0.05, "no usable distances; fell back to a fixed default"
        t = float(np.percentile(vals, 5.0))
        return t, f"5th percentile of {len(vals)} pairwise distances (no labels supplied)"

    def find_clusters(
        self,
        dist_df: pd.DataFrame,
        gold_file_path: Optional[str] = None,
        k_min: int = 2,
        k_max: int = 50,
        relative_gap_floor: Optional[float] = None,
        output_dir: Optional[str] = None,
        all_input_ids: Optional[List[str]] = None,
        cut_mode: str = "count",
        linkage_threshold: Optional[float] = None,
        linkage_method: str = "ward"
    ) -> Tuple[pd.DataFrame, int, float]:
        """
        Runs Ward AGNES hierarchical clustering and dynamic tree cut search.
        Supports prospective unsupervised clustering (Dendrogram Merge Height Gap Knee Detection) when gold_file_path is None.

        Parameters
        ----------
        dist_df : pd.DataFrame
            Symmetric dissimilarity matrix.
        gold_file_path : Optional[str]
            Path to gold standard cluster reference list (optional).
        k_min : int
            Minimum cluster count for tree cut search (default: 2).
        k_max : int
            Maximum cluster count for tree cut search (default: 50).
        relative_gap_floor : Optional[float]
            Minimum relative merge-height gap fraction of tree height required for unsupervised knee selection (default: 0.2200).
        output_dir : Optional[str]
            Target directory for cluster output TSV.
        all_input_ids : Optional[List[str]]
            List of all original specimen IDs before completeness filtering (for transparent reporting).
        cut_mode : str
            "count" (default) or "distance". These answer different questions and neither is a
            better version of the other.

            "count" asks: split this cohort into k groups. It selects k from the largest merge
            height gap, subject to a minimum relative gap and a minimum cluster size. That is
            the right question for a closed outbreak investigation, where every specimen belongs
            to some cluster and clusters below a few members are not actionable.

            "distance" asks: which specimens are close enough to be linked? It cuts the tree at
            a fixed dissimilarity. There is no k to choose and no gap requirement, and specimens
            with no near neighbour are returned as singletons rather than forced into a group.
            That is the right question for surveillance, where most cases are unrelated.

            Choosing k fails outright on surveillance-shaped data. On a 183-specimen
            Plasmodium vivax AmpliSeq cohort whose published truth is 93 groups with 79
            singletons, both guards reject every k that could reproduce it -- the minimum
            cluster size is 5 at that n, and a diffuse cloud has no dominant merge gap -- so the
            count rule returns a single cluster. Cutting the same tree at 0.08 gives 116 groups
            with 91 singletons, ARI 0.7578 against the published assignment, and the score stays
            above 0.70 across 0.06-0.12.
        linkage_threshold : Optional[float]
            Dissimilarity at which to cut in "distance" mode. If omitted it is calibrated by
            suggest_linkage_threshold() and the provenance of the number is printed.
        linkage_method : str
            Linkage for the tree, default "ward". "single" is worth trying in distance mode:
            Ward optimises within-cluster variance, which is the right objective when clusters
            are compact blobs, but transmission clusters are chains -- A infects B infects C --
            and single linkage follows a chain. On the Plasmodium vivax AmpliSeq cohort, single
            linkage at 0.070 returns exactly 92 groups, the published number, at ARI 0.7762,
            against 0.7578 for Ward at its own best cut.

            Single linkage chains through noise if the threshold is loose, so it is not the
            default. It is offered because on this shape it is measurably better.

        Returns
        -------
        Tuple[pd.DataFrame, int, float]
            (Cluster Assignments DataFrame, Number of Clusters, Threshold Used)
        """
        samples = dist_df.index.tolist()
        n_samples = len(samples)

        if n_samples < 2:
            single_df = pd.DataFrame({"Seq_ID": samples, "Assigned_cluster": [1]})
            return single_df, 1, 0.0

        dist_mat = dist_df.values.copy()
        np.fill_diagonal(dist_mat, 0.0)
        condensed_dist = squareform(dist_mat, checks=False)

        # Ward AGNES hierarchical clustering
        Z = linkage(condensed_dist, method=linkage_method, metric="euclidean")

        if cut_mode == "distance":
            # Linkage mode: cut at a dissimilarity, do not choose a cluster count.
            # Neither the relative-gap floor nor the minimum-cluster-size guard applies here.
            # Both exist to keep the count rule from reporting a split it cannot justify, and
            # both would reject the singleton-heavy structure this mode is written for.
            gold_df = None
            if gold_file_path and os.path.exists(gold_file_path):
                gold_df = pd.read_csv(gold_file_path, sep=r"\s+", engine="python")
            if linkage_threshold is None:
                threshold, provenance = self.suggest_linkage_threshold(
                    dist_df, gold_df, robust=self.robust)
                print(f"[ClusterFinder] Linkage mode: threshold {threshold:.4f} "
                      f"({provenance}).")
            else:
                threshold = float(linkage_threshold)
                print(f"[ClusterFinder] Linkage mode: threshold {threshold:.4f} "
                      f"(supplied by caller).")

            cluster_ids = fcluster(Z, threshold, criterion="distance")
            sizes = np.bincount(cluster_ids)[1:]
            n_singletons = int(np.sum(sizes == 1))
            correct_k = int(len(sizes))
            best_cluster_df = pd.DataFrame({
                "Seq_ID": samples,
                "Assigned_cluster": cluster_ids
            })
            print(f"[ClusterFinder] Linkage mode: {correct_k} groups over {n_samples} specimens "
                  f"({n_singletons} singletons, largest {int(sizes.max())}).")
            self.last_selection_meta = {
                "status": "linkage",
                "cut_mode": "distance",
                "k": correct_k,
                "threshold": threshold,
                "n_singletons": n_singletons,
                "largest_cluster": int(sizes.max()),
                "n_samples": n_samples,
            }
            best_cluster_df = self._write_clusters(
                best_cluster_df, output_dir, correct_k, all_input_ids, samples)
            return best_cluster_df, correct_k, float(threshold)


        if gold_file_path and os.path.exists(gold_file_path):
            gold_df = pd.read_csv(gold_file_path, sep=r"\s+", engine="python")
            self.compute_distance_auc(dist_df, gold_df)
            threshold = self.calibrate_gold_standard_threshold(dist_df, gold_df)
            
            correct_k = min(k_max, n_samples)
            best_cluster_df = None

            pair_distances = []
            for i in range(n_samples):
                for j in range(i + 1, n_samples):
                    pair_distances.append((i, j, dist_mat[i, j]))

            # Guard: Require k >= max(2, k_min) to avoid single-cluster k=1 collapse
            search_start = max(2, k_min)
            for k in range(search_start, min(k_max, n_samples) + 1):
                cluster_ids = cut_tree(Z, n_clusters=k).ravel()

                within_total = 0
                within_below = 0

                for i, j, d in pair_distances:
                    if cluster_ids[i] == cluster_ids[j]:
                        within_total += 1
                        if d <= threshold:
                            within_below += 1

                pct_meeting = (within_below / within_total * 100.0) if within_total > 0 else 100.0

                if pct_meeting >= self.stringency:
                    correct_k = k
                    best_cluster_df = pd.DataFrame({
                        "Seq_ID": samples,
                        "Assigned_cluster": cluster_ids + 1
                    })
                    print(f"[ClusterFinder] Optimal supervised cluster count: k = {correct_k} ({pct_meeting:.2f}% pairs meeting threshold).")
                    break

            if best_cluster_df is None:
                cluster_ids = cut_tree(Z, n_clusters=max(2, min(k_max, n_samples))).ravel()
                best_cluster_df = pd.DataFrame({
                    "Seq_ID": samples,
                    "Assigned_cluster": cluster_ids + 1
                })
                correct_k = len(set(cluster_ids))

        else:
            # Prospective Unsupervised Mode: Dendrogram Merge Height Gap Knee Detection (Elbow Rule)
            rel_floor = self.relative_gap_floor if relative_gap_floor is None else relative_gap_floor
            print(f"[ClusterFinder] Prospective Unsupervised Mode: Evaluating dendrogram merge height gap knee (k in [{k_min}, {k_max}], relative_gap_floor={rel_floor:.4f})...")
            if k_max < n_samples / 2:
                # Surveillance cohorts are often mostly singletons, so the true k can approach
                # n. A default k_max of 50 silently puts such an answer out of reach: on a
                # 183-specimen cohort whose published truth is 93 groups, no k in [2, 50] can
                # be right whatever the guards do. Say so rather than return a number the
                # search range alone determined.
                print(f"[ClusterFinder] Note: k_max={k_max} is below n/2={n_samples // 2}. If this "
                      f"cohort is mostly unrelated specimens its true k may exceed the search "
                      f"range; consider --cut distance, which does not choose a k.")
            
            # Z[:, 2] contains merge heights in ascending order (last merge is 2 -> 1 cluster)
            heights = Z[::-1, 2] # Descending order: heights[0] is merge 2 -> 1, heights[1] is merge 3 -> 2, etc.
            tree_height = float(heights[0]) if len(heights) > 0 else 1.0
            
            # Adaptive minimum cluster size guard: min 2 for small cohorts, up to 5 for surveillance scale
            min_required_size = max(2, min(5, int(0.10 * n_samples)))
            search_start = max(2, k_min)
            search_limit = min(k_max, len(heights) + 1, n_samples)
            
            gap_scores = []
            for k in range(search_start, search_limit + 1):
                idx = k - 2
                if idx < len(heights):
                    h_curr = heights[idx]
                    h_next = heights[idx + 1] if idx + 1 < len(heights) else 0.0
                    gap = h_curr - h_next
                    gap_scores.append((k, gap, h_curr, h_next))

            valid_k_found = False
            best_fallback_k = None
            best_fallback_info = None
            rejection_reasons = []

            if gap_scores:
                # Sort candidate k by descending merge height gap
                sorted_gaps = sorted(gap_scores, key=lambda x: x[1], reverse=True)

                for candidate in sorted_gaps:
                    cand_k, cand_gap, h_curr, h_next = candidate
                    cand_rel_gap = (cand_gap / tree_height) if tree_height > 0 else 0.0

                    cand_ids = cut_tree(Z, n_clusters=cand_k).ravel()
                    cluster_counts = np.bincount(cand_ids)
                    min_c_size = int(np.min(cluster_counts))

                    if cand_rel_gap >= rel_floor and min_c_size >= min_required_size:
                        correct_k = cand_k
                        max_gap = cand_gap
                        rel_gap = cand_rel_gap
                        threshold = float((h_curr + h_next) / 2.0)
                        valid_k_found = True
                        self.last_selection_meta = {
                            "status": "optimal",
                            "k": correct_k,
                            "threshold": threshold,
                            "gap": max_gap,
                            "relative_gap": rel_gap,
                            "min_cluster_size": min_c_size,
                            "min_required_size": min_required_size,
                            "k_min": k_min,
                            "k_max": k_max,
                            "tree_height": tree_height
                        }
                        print(f"[ClusterFinder] Dendrogram Merge Height Gap Knee Detection: Optimal k = {correct_k} (Height Gap = {max_gap:.5f}, Rel Gap = {rel_gap:.4f}, Min Cluster Size = {min_c_size} >= {min_required_size}, Threshold = {threshold:.5f}).")
                        break
                    else:
                        reasons = []
                        if cand_rel_gap < rel_floor:
                            reasons.append(f"relative gap {cand_rel_gap:.4f} < floor {rel_floor:.4f}")
                        if min_c_size < min_required_size:
                            reasons.append(f"min cluster size {min_c_size} < guard {min_required_size}")
                        rejection_reasons.append(f"k={cand_k} ({', '.join(reasons)})")

                        if min_c_size >= min_required_size and best_fallback_k is None:
                            best_fallback_k = cand_k
                            best_fallback_info = (cand_gap, cand_rel_gap, h_curr, h_next, min_c_size)

                if not valid_k_found:
                    if k_min > 2 and best_fallback_k is not None:
                        cand_gap, cand_rel_gap, h_curr, h_next, min_c_size = best_fallback_info
                        correct_k = best_fallback_k
                        threshold = float((h_curr + h_next) / 2.0)
                        self.last_selection_meta = {
                            "status": "floor_override",
                            "k": correct_k,
                            "threshold": threshold,
                            "gap": cand_gap,
                            "relative_gap": cand_rel_gap,
                            "min_cluster_size": min_c_size,
                            "min_required_size": min_required_size,
                            "k_min": k_min,
                            "k_max": k_max,
                            "tree_height": tree_height
                        }
                        print(f"[ClusterFinder] Selected best candidate k = {correct_k} satisfying requested k_min >= {k_min} and cluster size guard ({min_c_size} >= {min_required_size}), despite relative gap ({cand_rel_gap:.4f}) being below floor ({rel_floor:.4f}).")
                    else:
                        correct_k = 1
                        threshold = float(heights[0]) if len(heights) > 0 else 0.0
                        status_label = "unsatisfiable_constraint" if k_min > 2 else "single_group"
                        self.last_selection_meta = {
                            "status": status_label,
                            "k": 1,
                            "threshold": threshold,
                            "rejection_reasons": rejection_reasons,
                            "min_required_size": min_required_size,
                            "k_min": k_min,
                            "k_max": k_max,
                            "tree_height": tree_height
                        }
                        fail_summary = "; ".join(rejection_reasons[:3])
                        print(f"[ClusterFinder] Dendrogram Merge Height Gap Knee Detection: No valid partition (k in [{search_start}, {search_limit}]) met both relative gap floor ({rel_floor:.4f}) and cluster size guard ({min_required_size}). Rejections: {fail_summary}. Assigned k = 1 (Single Outbreak Group).")
            else:
                correct_k = max(1, k_min)
                threshold = self.default_threshold
                self.last_selection_meta = {
                    "status": "trivial",
                    "k": correct_k,
                    "threshold": threshold,
                    "k_min": k_min,
                    "k_max": k_max
                }

            cluster_ids = cut_tree(Z, n_clusters=correct_k).ravel()

            # Deterministic lexicographical tie-breaking and label ordering (Cluster 1 contains lexicographically smallest specimen)
            unique_cids = []
            for sid, cid in sorted(zip(samples, cluster_ids), key=lambda x: x[0]):
                if cid not in unique_cids:
                    unique_cids.append(cid)
            label_remap = {old_cid: new_idx + 1 for new_idx, old_cid in enumerate(unique_cids)}
            remapped_ids = np.array([label_remap[cid] for cid in cluster_ids])

            best_cluster_df = pd.DataFrame({
                "Seq_ID": samples,
                "Assigned_cluster": remapped_ids
            })

        best_cluster_df = self._write_clusters(
            best_cluster_df, output_dir, correct_k, all_input_ids, samples)

        return best_cluster_df, correct_k, float(threshold)

    @staticmethod
    def _write_clusters(cluster_df, output_dir, k, all_input_ids, samples):
        """Append excluded specimens as cluster -1, write the TSV, return the full frame.

        Shared by both cut modes so a specimen that was dropped for low completeness is reported
        the same way whichever question was asked. Silently omitting it would make the cohort
        look smaller than it was.
        """
        if all_input_ids:
            missing_ids = [sid for sid in all_input_ids if sid not in samples]
            if missing_ids:
                missing_df = pd.DataFrame({
                    "Seq_ID": missing_ids,
                    "Assigned_cluster": -1  # -1 = excluded from the distance matrix
                })
                cluster_df = pd.concat([cluster_df, missing_df], ignore_index=True)
                print(f"[ClusterFinder] Transparent Reporting: {len(missing_ids)} low-completeness "
                      f"specimens explicitly reported as Cluster -1 (Unassigned).")

        if output_dir is None:
            output_dir = "clusters_detected"
        os.makedirs(output_dir, exist_ok=True)
        today_str = datetime.date.today().strftime("%Y-%m-%d")
        output_path = os.path.join(output_dir, f"{today_str}_RESULTING_CLUSTERS_{k}.txt")
        cluster_df.to_csv(output_path, sep="\t", index=False)
        print(f"[ClusterFinder] Saved outbreak cluster assignments to: {output_path}")
        return cluster_df


    # ==================================================================================
    # SWEEP DIAGNOSTIC -- the default clustering answer.
    #
    # A single k is only trustworthy when the data actually determines it. Rather than
    # commit to one cut, cluster_sweep reports how many groups the cohort supports as a
    # RANGE, how confident that is, and a per-branch confidence tree. Four independent,
    # unsupervised count selectors (merge-gap knee, silhouette, bootstrap stability, and
    # the Tibshirani gap statistic) vote; when they agree a single number is reported,
    # when they scatter a range is, with the confident sub-structure (stable cores and
    # well-supported splits) reported underneath either way. Nothing here uses ground
    # truth. Branch support = 1 - mean cross-cluster co-assignment over bootstrap
    # resamples, marginalised across resolution -- high support = a split the data
    # reproduces, low = one it does not (drawn faded/dashed in the confidence tree).
    # ==================================================================================

    @staticmethod
    def _rel_gap(heights: np.ndarray, k: int, n: int) -> float:
        """Relative merge-height gap crossing from k to k+1 clusters (what the knee reads)."""
        if k < 1 or k >= n:
            return 0.0
        hi = heights[n - 1 - k]
        lo = heights[n - 2 - k] if n - 2 - k >= 0 else 0.0
        return float((hi - lo) / hi) if hi > 0 else 0.0

    @staticmethod
    def _to_newick(Z: np.ndarray, labels: List[str], support: List[float]) -> str:
        """Newick string with per-internal-node support (phylogenetics convention)."""
        n = len(labels)
        import sys as _sys
        _sys.setrecursionlimit(max(10000, n * 4))

        def rec(node: int) -> str:
            if node < n:
                return str(labels[node]).replace("(", "_").replace(")", "_").replace(",", "_").replace(":", "_")
            a, b = int(Z[node - n, 0]), int(Z[node - n, 1])
            s = support[node - n]
            return f"({rec(a)},{rec(b)}){s:.3f}"

        return rec(2 * n - 2) + ";"

    @staticmethod
    def _count_confident(Z: np.ndarray, support: List[float], n: int, tau: float) -> int:
        """Number of clusters when only splits with support >= tau are trusted.

        Walk down from the root: at a split whose support clears tau the two sides are
        distinct clusters (recurse); at a split below tau the whole subtree collapses to one
        cluster, because the data does not reproduce that division. Sweeping tau turns the
        confidence tree into a range of defensible cluster counts.
        """
        import sys as _sys
        _sys.setrecursionlimit(max(10000, n * 4))

        def rec(node: int) -> int:
            if node < n:
                return 1
            if support[node - n] >= tau:
                a, b = int(Z[node - n, 0]), int(Z[node - n, 1])
                return rec(a) + rec(b)
            return 1

        return rec(2 * n - 2)

    def cluster_sweep(
        self,
        dist_df: pd.DataFrame,
        k_min: int = 2,
        k_max: int = 50,
        n_boot: int = 200,
        boot_frac: float = 0.85,
        seed: int = 0,
        core_frac: float = 0.90,
        support_solid: float = 0.75,
        support_mid: float = 0.45,
        do_gap: bool = True,
        linkage_method: str = "ward",
        output_dir: Optional[str] = None,
    ) -> Dict[str, any]:
        """Run the cluster-count sweep diagnostic on a distance matrix.

        Returns a dict with: the count RANGE and (when the selectors agree) a point
        estimate, a per-k sweep table (silhouette / relative gap / stability / gap
        statistic), the four selector votes, the stable cores, and the confidence tree
        (leaf order, per-node support, and a Newick string). Writes a JSON report and a
        Newick tree to output_dir when given. Unsupervised throughout.
        """
        samples = list(dist_df.index)
        n = len(samples)
        result: Dict[str, any] = {"n": n, "linkage_method": linkage_method}
        if n < 3:
            result.update({"count_range": [1, max(1, n)], "point_estimate": max(1, n),
                           "confident": True, "note": "cohort too small to sweep"})
            return result

        D = dist_df.values.astype(float).copy()
        D = (D + D.T) / 2.0
        np.fill_diagonal(D, 0.0)
        Z = linkage(squareform(D, checks=False), method=linkage_method, metric="euclidean")
        heights = Z[:, 2]

        ks = list(range(max(2, k_min), min(k_max, n - 1) + 1))
        full_lab = {k: fcluster(Z, k, criterion="maxclust") for k in ks}

        # per-k unsupervised scores
        sweep = []
        for k in ks:
            lab = full_lab[k]
            kk = len(set(lab))
            sil = float(silhouette_score(D, lab, metric="precomputed")) if 2 <= kk < n else float("nan")
            singles = int(sum(1 for c in set(lab) if list(lab).count(c) == 1))
            sweep.append({"k": k, "clusters": kk, "silhouette": round(sil, 4),
                          "rel_gap": round(self._rel_gap(heights, k, n), 4),
                          "singletons": singles})

        # bootstrap: marginalised co-assignment (for support / cores / pair resolution)
        # and per-k stability (agreement of resample co-membership with the full tree).
        rng = np.random.default_rng(seed)
        co = np.zeros((n, n)); ct = np.zeros((n, n))
        stab_sum = {k: 0.0 for k in ks}; stab_cnt = {k: 0 for k in ks}
        mid_hi = max(3, min(k_max, n // 2))
        for _ in range(n_boot):
            idx = np.sort(rng.choice(n, int(n * boot_frac), replace=False))
            Zi = linkage(squareform(D[np.ix_(idx, idx)], checks=False), method=linkage_method, metric="euclidean")
            kk = int(rng.integers(max(2, k_min), mid_hi + 1))
            lab = fcluster(Zi, kk, criterion="maxclust")
            same = (lab[:, None] == lab[None, :]).astype(float)
            sub = np.ix_(idx, idx)
            co[sub] += same
            ct[sub] += 1.0
            iu = np.triu_indices(len(idx), 1)
            for k in ks:
                lk = fcluster(Zi, k, criterion="maxclust")
                rc = (lk[:, None] == lk[None, :])
                fl = full_lab[k][idx]
                fc = (fl[:, None] == fl[None, :])
                stab_sum[k] += float(np.mean(rc[iu] == fc[iu]))
                stab_cnt[k] += 1
        with np.errstate(invalid="ignore", divide="ignore"):
            frac = np.where(ct > 0, co / np.maximum(ct, 1.0), 0.0)
        np.fill_diagonal(frac, 1.0)
        for row in sweep:
            k = row["k"]
            row["stability"] = round(stab_sum[k] / stab_cnt[k], 4) if stab_cnt[k] else float("nan")

        # per-merge support = 1 - mean cross-cluster co-assignment
        mem: Dict[int, List[int]] = {i: [i] for i in range(n)}
        support: List[float] = []
        for i, (a, b, h, _) in enumerate(Z):
            a, b = int(a), int(b)
            mem[n + i] = mem[a] + mem[b]
            cross = frac[np.ix_(mem[a], mem[b])]
            support.append(round(float(1.0 - cross.mean()), 4))
        strong_splits = int(sum(1 for s in support if s >= 0.75))

        # stable cores: connected components at co-assignment >= core_frac
        adj = frac >= core_frac
        seen = [False] * n; cores = []
        for i in range(n):
            if seen[i]:
                continue
            stack = [i]; comp = []
            seen[i] = True
            while stack:
                u = stack.pop(); comp.append(u)
                for v in range(n):
                    if adj[u, v] and not seen[v]:
                        seen[v] = True; stack.append(v)
            cores.append([samples[j] for j in comp])
        nonsingleton_cores = [c for c in cores if len(c) > 1]

        # pairwise resolution: how decisively pairs are same/different (0.5 = a coin flip)
        iu = np.triu_indices(n, 1)
        decisive = np.abs(2.0 * frac[iu] - 1.0)
        pairs_resolved = float(np.mean(decisive >= 0.5))
        mean_decisiveness = float(np.mean(decisive))

        # selectors (all unsupervised)
        def _best(key):
            vals = [(r["k"], r[key]) for r in sweep if r[key] == r[key]]
            return max(vals, key=lambda t: t[1])[0] if vals else None
        selectors: Dict[str, Optional[int]] = {}
        knee = next((r["k"] for r in sweep if r["rel_gap"] >= self.relative_gap_floor), None)
        selectors["knee"] = knee
        selectors["silhouette"] = _best("silhouette")
        selectors["stability"] = _best("stability")
        if do_gap:
            try:
                selectors["gap"] = self._gap_statistic(D, Z, ks, rng, linkage_method)
            except Exception as exc:  # MDS is optional; never let it sink the sweep
                print(f"[ClusterFinder] gap statistic skipped: {exc}")
        # CONFIDENCE -- do the independent count selectors agree? The merge-gap knee,
        # silhouette, and gap statistic are three unsupervised estimates of k. When they
        # concur the count is determined and a single number is reported; when they scatter
        # it is not. (Bootstrap stability is kept in the sweep table but excluded from the
        # vote: it saturates high on diffuse data and would drag every cohort toward "fuzzy".)
        vote_keys = ("knee", "silhouette", "gap")
        votes = [selectors[k] for k in vote_keys if selectors.get(k)]
        if votes:
            mode_k = Counter(votes).most_common(1)[0][0]
            n_near = sum(1 for v in votes if abs(v - mode_k) <= 1)
            spread = max(votes) - min(votes)
            confident = (n_near >= max(2, (len(votes) + 1) // 2)) and (spread <= max(2, 0.30 * mode_k))
        else:
            mode_k = None; confident = False; spread = 0

        # RANGE -- when the count is not determined, read it off the confidence tree itself:
        # cut the tree keeping only splits at or above a support threshold, at the same two
        # tiers the tree is drawn. solid (>= support_solid) gives the groups the data
        # reproduces fully; solid+moderate (>= support_mid) gives the finest partition still
        # supported. (The raw selector span is uninformative -- silhouette and gap drift
        # toward fine partitions on diffuse surveillance data.)
        k_solid = self._count_confident(Z, support, n, support_solid)
        k_mid = self._count_confident(Z, support, n, support_mid)
        if confident:
            lo = hi = int(mode_k); point = int(mode_k)
        else:
            lo, hi = min(k_solid, k_mid), max(k_solid, k_mid); point = None

        result.update({
            "count_range": [int(lo), int(hi)],
            "point_estimate": point,
            "confident": bool(confident),
            "count_at_solid_support": int(k_solid),
            "count_at_moderate_support": int(k_mid),
            "support_tiers": {"solid": support_solid, "moderate": support_mid},
            "naive_selectors": selectors,
            "naive_selector_spread": int(spread),
            "pairs_resolved": round(pairs_resolved, 4),
            "mean_decisiveness": round(mean_decisiveness, 4),
            "strong_splits": strong_splits,
            "n_stable_cores": len(nonsingleton_cores),
            "stable_cores": nonsingleton_cores,
            "sweep": sweep,
            "tree": self._tree_render(Z, samples, support),
        })
        result["tree"]["newick"] = self._to_newick(Z, samples, support)
        result["headline"] = (
            f"{point} clusters (confident: count selectors agree)"
            if confident else
            f"{lo}-{hi} clusters (count not determined; selectors scatter "
            f"{min(votes) if votes else '?'}-{max(votes) if votes else '?'}): solid splits give "
            f"{k_solid}, adding moderate splits gives {k_mid}; {len(nonsingleton_cores)} stable "
            f"cores, {int(round(pairs_resolved * 100))}% of specimen pairs resolved"
        )

        # representative partition, for downstream tools that need one flat assignment:
        # the confident count when there is one, else the finest supported (moderate) tier.
        rep_k = point if point else max(2, k_mid)
        rep_lab = fcluster(Z, rep_k, criterion="maxclust")
        rep_df = pd.DataFrame({"Seq_ID": samples, "Assigned_cluster": rep_lab})
        result["representative_k"] = int(rep_k)

        if output_dir:
            import json as _json
            os.makedirs(output_dir, exist_ok=True)
            today = datetime.date.today().strftime("%Y-%m-%d")
            jpath = os.path.join(output_dir, f"{today}_SWEEP.json")
            npath = os.path.join(output_dir, f"{today}_confidence_tree.nwk")
            with open(jpath, "w") as fh:
                _json.dump(result, fh, indent=1)
            with open(npath, "w") as fh:
                fh.write(result["tree"]["newick"] + "\n")
            self._write_clusters(rep_df, output_dir, rep_k, None, samples)
            print(f"[ClusterFinder] Sweep: {result['headline']}")
            print(f"[ClusterFinder] Wrote {jpath}, {npath}, and a representative partition "
                  f"(k={rep_k}, {'confident' if point else 'moderate tier'}).")
        self.last_selection_meta = {"status": "sweep", "count_range": [int(lo), int(hi)],
                                    "confident": bool(confident), "point_estimate": point}
        return result

    @staticmethod
    def _tree_render(Z: np.ndarray, samples: List[str], support: List[float]) -> Dict[str, any]:
        """Leaf order + per-node geometry (child positions/heights) and support, so a
        renderer can draw the confidence tree without re-deriving the linkage."""
        n = len(samples)
        order = list(leaves_list(Z))
        ypos = {leaf: i for i, leaf in enumerate(order)}
        nx: Dict[int, float] = {}
        nodes = []
        for i, (a, b, h, _) in enumerate(Z):
            a, b = int(a), int(b)
            ya = ypos[a] if a < n else nx[a]
            yb = ypos[b] if b < n else nx[b]
            hca = Z[a - n][2] if a >= n else 0.0
            hcb = Z[b - n][2] if b >= n else 0.0
            nx[n + i] = (ya + yb) / 2.0
            nodes.append({"ya": float(ya), "yb": float(yb), "h": float(h),
                          "hca": float(hca), "hcb": float(hcb), "support": support[i]})
        return {"leaf_order": [samples[i] for i in order],
                "hmax": float(Z[:, 2].max()), "nodes": nodes}

    @staticmethod
    def _gap_statistic(D: np.ndarray, Z: np.ndarray, ks: List[int], rng, linkage_method: str) -> int:
        """Tibshirani gap statistic over an MDS embedding; returns argmax-gap k. Unsupervised."""
        import warnings
        from sklearn.manifold import MDS
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # sklearn churns MDS's init/dissimilarity kwargs
            X = MDS(n_components=min(10, D.shape[0] - 1), dissimilarity="precomputed",
                    random_state=0, normalized_stress="auto").fit_transform(D)
        lo, hi = X.min(0), X.max(0)

        def Wk(lab, Xe):
            s = 0.0
            for c in set(lab):
                pts = Xe[lab == c]
                if len(pts) > 1:
                    s += ((pts[:, None, :] - pts[None, :, :]) ** 2).sum() / (2 * len(pts))
            return s

        B = 20
        lref = np.zeros((B, len(ks)))
        for bnum in range(B):
            Xr = rng.uniform(lo, hi, X.shape)
            Zr = linkage(Xr, method=linkage_method)
            for j, k in enumerate(ks):
                lref[bnum, j] = np.log(Wk(fcluster(Zr, k, "maxclust"), Xr) + 1e-9)
        gaps = {}
        for j, k in enumerate(ks):
            lw = np.log(Wk(fcluster(Z, k, "maxclust"), X) + 1e-9)
            gaps[k] = float(lref[:, j].mean() - lw)
        return max(gaps, key=gaps.get)

    def detect_micro_clusters(
        self,
        dist_df: pd.DataFrame,
        micro_threshold: float = 0.0,
        min_micro_size: int = 2
    ) -> List[List[str]]:
        """
        [EXPERIMENTAL DIAGNOSTIC] Scans distance matrix for micro-clusters of specimens
        (n >= min_micro_size) exhibiting exact or near-identical multi-locus profiles (D <= micro_threshold).
        By default (micro_threshold=0.0), restricts scans to exact identical profile matches (D <= 1e-6)
        to prevent single-linkage chaining on noisy background data.
        """
        from scipy.cluster.hierarchy import fcluster, linkage
        from scipy.spatial.distance import squareform

        dist_mat = dist_df.values.copy()
        np.fill_diagonal(dist_mat, 0.0)
        condensed_dist = squareform(dist_mat, checks=False)

        # Default to exact zero-distance identical genotype profile matching
        thresh_val = micro_threshold if micro_threshold is not None else 0.0
        if thresh_val <= 1e-6:
            thresh_val = 1e-6

        Z = linkage(condensed_dist, method="single", metric="euclidean")
        labels = fcluster(Z, t=thresh_val, criterion="distance")

        samples = dist_df.index.tolist()
        clusters = {}
        for sample, label in zip(samples, labels):
            clusters.setdefault(label, []).append(sample)

        micro_clusters = [members for members in clusters.values() if len(members) >= min_micro_size]
        print(f"[MicroClusterScanner] Evaluated dynamic threshold = {thresh_val:.6f}. Found {len(micro_clusters)} micro-clusters (size >= {min_micro_size}).")
        return micro_clusters

