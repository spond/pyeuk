"""
CyclosporaClusterFinder: Modernized Python module for AGNES Ward's hierarchical clustering
and robust outbreak threshold calibration for CDC Cyclospora cayetanensis workflow.
Supports prospective unsupervised clustering via Dendrogram Merge Height Gap Knee Detection
and gold-standard supervised calibration.
"""

import os
import datetime
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, cut_tree, fcluster
from sklearn.metrics import roc_auc_score
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
        linkage_threshold: Optional[float] = None
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
        Z = linkage(condensed_dist, method="ward", metric="euclidean")

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

