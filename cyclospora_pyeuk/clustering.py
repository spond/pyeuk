"""
CyclosporaClusterFinder: Modernized Python module for AGNES Ward's hierarchical clustering
and robust outbreak threshold calibration for CDC Cyclospora cayetanensis workflow.
Supports prospective unsupervised clustering via Silhouette score optimization and gold-standard supervised calibration.
"""

import os
import datetime
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, cut_tree
from sklearn.metrics import silhouette_score
from typing import Tuple, Optional, Dict, List


class CyclosporaClusterFinder:
    """
    Hierarchical clustering engine for outbreak surveillance, supporting prospective unsupervised
    clustering (via Silhouette score maximization) and gold-standard supervised calibration.
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
        k_min: int = 1,
        k_max: int = 50,
        output_dir: Optional[str] = None,
        all_input_ids: Optional[List[str]] = None
    ) -> Tuple[pd.DataFrame, int, float]:
        """
        Runs Ward AGNES hierarchical clustering and dynamic tree cut search.
        Supports prospective unsupervised clustering (Silhouette score optimization) when gold_file_path is None.

        Parameters
        ----------
        dist_df : pd.DataFrame
            Symmetric dissimilarity matrix.
        gold_file_path : Optional[str]
            Path to gold standard cluster reference list (optional).
        k_min : int
            Minimum cluster count for tree cut search (default: 1).
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

        dist_mat = dist_df.values.copy()
        np.fill_diagonal(dist_mat, 0.0)
        condensed_dist = squareform(dist_mat, checks=False)

        # Ward AGNES hierarchical clustering
        Z = linkage(condensed_dist, method="ward", metric="euclidean")

        if gold_file_path and os.path.exists(gold_file_path):
            gold_df = pd.read_csv(gold_file_path, sep=r"\s+", engine="python")
            threshold = self.calibrate_gold_standard_threshold(dist_df, gold_df)
            
            correct_k = min(k_max, n_samples)
            best_cluster_df = None

            pair_distances = []
            for i in range(n_samples):
                for j in range(i + 1, n_samples):
                    pair_distances.append((i, j, dist_mat[i, j]))

            for k in range(k_min, min(k_max, n_samples) + 1):
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
                cluster_ids = cut_tree(Z, n_clusters=min(k_max, n_samples)).ravel()
                best_cluster_df = pd.DataFrame({
                    "Seq_ID": samples,
                    "Assigned_cluster": cluster_ids + 1
                })
                correct_k = min(k_max, n_samples)

        else:
            # Prospective Unsupervised Mode: Optimize Silhouette score over k
            print("[ClusterFinder] Prospective Unsupervised Mode: Optimizing Silhouette score over tree cuts...")
            best_k = 2
            best_sil = -1.0

            search_max = min(k_max, n_samples - 1)
            for k in range(2, max(3, search_max + 1)):
                labels = cut_tree(Z, n_clusters=k).ravel()
                sil = silhouette_score(dist_mat, labels, metric="precomputed")
                if sil > best_sil:
                    best_sil = sil
                    best_k = k

            correct_k = best_k
            cluster_ids = cut_tree(Z, n_clusters=correct_k).ravel()

            # Estimate threshold as merge height at k
            heights = Z[:, 2]
            threshold = float(heights[-(correct_k - 1)]) if len(heights) >= correct_k - 1 else self.default_threshold

            best_cluster_df = pd.DataFrame({
                "Seq_ID": samples,
                "Assigned_cluster": cluster_ids + 1
            })
            print(f"[ClusterFinder] Unsupervised Silhouette Maximization: Optimal k = {correct_k} (Silhouette Score = {best_sil:.4f}, Height Threshold = {threshold:.5f}).")

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
