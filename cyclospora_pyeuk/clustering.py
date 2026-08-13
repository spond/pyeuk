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
from scipy.cluster.hierarchy import linkage, cut_tree
from sklearn.metrics import roc_auc_score
from typing import Tuple, Optional, Dict, List


class CyclosporaClusterFinder:
    """
    Hierarchical clustering engine for outbreak surveillance, supporting prospective unsupervised
    clustering (via Dendrogram Merge Height Gap Knee Detection) and gold-standard supervised calibration.
    """

    def __init__(self, stringency: float = 95.0, robust: bool = True, default_threshold: float = 0.05):
        """
        Parameters
        ----------
        stringency : float
            Percentage of within-cluster distances that must fall below the distance threshold (default: 95.0%).
        robust : bool
            If True, uses Median + 3 * 1.4826 * MAD for threshold calibration. If False, uses Mean + 3 * StdDev.
        default_threshold : float
            Default threshold for prospective unsupervised clustering when gold standards are omitted.
        """
        self.stringency = stringency
        self.robust = robust
        self.default_threshold = default_threshold

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

    def find_clusters(
        self,
        dist_df: pd.DataFrame,
        gold_file_path: Optional[str] = None,
        k_min: int = 2,
        k_max: int = 50,
        output_dir: Optional[str] = None,
        all_input_ids: Optional[List[str]] = None
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
            Minimum cluster count for tree cut search (default: 2 to guard against k=1 single-cluster collapse).
        k_max : int
            Maximum cluster count for tree cut search (default: 50).
        output_dir : Optional[str]
            Target directory for cluster output TSV.
        all_input_ids : Optional[List[str]]
            List of all original specimen IDs before completeness filtering (for transparent reporting).

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
            print("[ClusterFinder] Prospective Unsupervised Mode: Evaluating dendrogram merge height gap knee...")
            
            # Z[:, 2] contains merge heights in ascending order (last merge is 2 -> 1 cluster)
            heights = Z[::-1, 2] # Descending order: heights[0] is merge 2 -> 1, heights[1] is merge 3 -> 2, etc.
            
            search_limit = min(k_max, len(heights))
            gap_scores = []

            for idx in range(search_limit - 1):
                k = idx + 2 # k clusters
                h_curr = heights[idx]
                h_next = heights[idx + 1]
                gap = h_curr - h_next
                gap_scores.append((k, gap, h_curr, h_next))

            if gap_scores:
                # Sort candidate k by descending merge height gap
                sorted_gaps = sorted(gap_scores, key=lambda x: x[1], reverse=True)
                tree_height = float(heights[0]) if len(heights) > 0 else 1.0
                
                # Minimum cluster size guard: at least 10% of samples or min 5 specimens
                min_required_size = max(5, int(0.10 * n_samples))
                
                valid_k_found = False
                for candidate in sorted_gaps:
                    cand_k, cand_gap, h_curr, h_next = candidate
                    cand_rel_gap = (cand_gap / tree_height) if tree_height > 0 else 0.0
                    
                    if cand_rel_gap < 0.2200:
                        break  # Remaining gap scores are below scale-free noise floor
                        
                    # Check cluster size distribution for candidate k
                    cand_ids = cut_tree(Z, n_clusters=cand_k).ravel()
                    cluster_counts = np.bincount(cand_ids)
                    min_c_size = int(np.min(cluster_counts))
                    
                    if min_c_size >= min_required_size:
                        correct_k = cand_k
                        max_gap = cand_gap
                        rel_gap = cand_rel_gap
                        threshold = float((h_curr + h_next) / 2.0)
                        valid_k_found = True
                        print(f"[ClusterFinder] Dendrogram Merge Height Gap Knee Detection: Optimal k = {correct_k} (Height Gap = {max_gap:.5f}, Rel Gap = {rel_gap:.4f}, Min Cluster Size = {min_c_size} >= {min_required_size}, Threshold = {threshold:.5f}).")
                        break
                    else:
                        print(f"[ClusterFinder] Outlier Split Guard: Rejecting candidate k = {cand_k} (Min Cluster Size = {min_c_size} < {min_required_size}). Searching for valid partition...")

                if not valid_k_found:
                    correct_k = 1
                    threshold = float(heights[0]) if len(heights) > 0 else 0.0
                    print(f"[ClusterFinder] Dendrogram Merge Height Gap Knee Detection: No valid outbreak partition found meeting minimum cluster size ({min_required_size}). Assigned k = 1 (Single Outbreak Group).")
            else:
                correct_k = 2
                threshold = self.default_threshold

            cluster_ids = cut_tree(Z, n_clusters=correct_k).ravel()
            best_cluster_df = pd.DataFrame({
                "Seq_ID": samples,
                "Assigned_cluster": cluster_ids + 1
            })

        # Transparently append low-completeness / excluded specimens if provided
        if all_input_ids:
            missing_ids = [sid for sid in all_input_ids if sid not in samples]
            if missing_ids:
                missing_df = pd.DataFrame({
                    "Seq_ID": missing_ids,
                    "Assigned_cluster": -1  # Cluster -1 indicates Low Completeness / Excluded from distance matrix
                })
                best_cluster_df = pd.concat([best_cluster_df, missing_df], ignore_index=True)
                print(f"[ClusterFinder] Transparent Reporting: {len(missing_ids)} low-completeness specimens explicitly reported as Cluster -1 (Unassigned).")

        if output_dir is None:
            output_dir = "clusters_detected"
        os.makedirs(output_dir, exist_ok=True)
        today_str = datetime.date.today().strftime("%Y-%m-%d")
        output_path = os.path.join(output_dir, f"{today_str}_RESULTING_CLUSTERS_{correct_k}.txt")

        best_cluster_df.to_csv(output_path, sep="\t", index=False)
        print(f"[ClusterFinder] Saved outbreak cluster assignments to: {output_path}")

        return best_cluster_df, correct_k, threshold
