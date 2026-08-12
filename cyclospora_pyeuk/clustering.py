"""
CyclosporaClusterFinder: Modernized Python module for AGNES Ward's hierarchical clustering
and robust outbreak threshold calibration for CDC Cyclospora cayetanensis workflow.
Replaces legacy R scripts cluster_counter.R and CLUSTER_FINDER.R.
"""

import os
import datetime
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, cut_tree
from typing import Tuple, Optional, Dict, List


class CyclosporaClusterFinder:
    """
    Hierarchical clustering and outbreak threshold calibration engine.
    """

    def __init__(self, stringency: float = 95.0, robust: bool = True):
        """
        Parameters
        ----------
        stringency : float
            Percentage of within-cluster distances that must fall below the gold standard threshold (default: 95.0%).
        robust : bool
            If True, uses Median + 3 * 1.4826 * MAD for threshold calibration. If False, uses Mean + 3 * StdDev.
        """
        self.stringency = stringency
        self.robust = robust

    def calibrate_gold_standard_threshold(
        self,
        dist_df: pd.DataFrame,
        gold_df: pd.DataFrame
    ) -> float:
        """
        Calculates maximum allowed intra-cluster distance threshold using epidemiological gold standards.
        """
        # Ensure gold_df has Seq_ID and Cluster_alias columns
        gold_df.columns = [c.strip() for c in gold_df.columns]
        if "Seq_ID" not in gold_df.columns:
            gold_df.rename(columns={gold_df.columns[0]: "Seq_ID"}, inplace=True)
        if "Cluster_alias" not in gold_df.columns and len(gold_df.columns) > 1:
            gold_df.rename(columns={gold_df.columns[1]: "Cluster_alias"}, inplace=True)

        # Intersect available samples
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
            print("[ClusterFinder] Warning: No overlapping gold standard sample pairs found. Falling back to default threshold (0.15).")
            return 0.15

        within_distances = np.array(within_distances)

        if self.robust:
            med = np.median(within_distances)
            mad = np.median(np.abs(within_distances - med))
            threshold = med + 3.0 * 1.4826 * mad
            print(f"[ClusterFinder] Robust calibration (Median + 3*MAD): Threshold = {threshold:.5f} (median={med:.5f}, MAD={mad:.5f})")
        else:
            mean = np.mean(within_distances)
            std = np.std(within_distances, ddof=1) if len(within_distances) > 1 else 0.0
            threshold = mean + 3.0 * std
            print(f"[ClusterFinder] Classical calibration (Mean + 3*StdDev): Threshold = {threshold:.5f} (mean={mean:.5f}, std={std:.5f})")

        return float(threshold)

    def find_clusters(
        self,
        dist_df: pd.DataFrame,
        gold_file_path: str,
        k_min: int = 1,
        k_max: int = 50,
        output_dir: Optional[str] = None
    ) -> Tuple[pd.DataFrame, int, float]:
        """
        Runs Ward AGNES hierarchical clustering and dynamic tree cut search.

        Returns
        -------
        Tuple[pd.DataFrame, int, float]
            (Cluster Assignments DataFrame, Correct Number of Clusters, Threshold Used)
        """
        # Read gold standard clusters file
        gold_df = pd.read_csv(gold_file_path, sep=r"\s+", engine="python")
        threshold = self.calibrate_gold_standard_threshold(dist_df, gold_df)

        samples = dist_df.index.tolist()
        n_samples = len(samples)

        # Convert distance matrix to condensed vector for linkage
        dist_mat = dist_df.values.copy()
        np.fill_diagonal(dist_mat, 0.0)
        condensed_dist = squareform(dist_mat, checks=False)

        # Ward AGNES hierarchical clustering (scipy.cluster.hierarchy.linkage method='ward')
        Z = linkage(condensed_dist, method="ward", metric="euclidean")

        # Evaluate cuttree for k in k_min..k_max
        correct_k = k_max
        best_cluster_df = None

        # Build pair distance lookup
        pair_distances = []
        for i in range(n_samples):
            for j in range(i + 1, n_samples):
                pair_distances.append((i, j, dist_mat[i, j]))

        for k in range(k_min, k_max + 1):
            cluster_ids = cut_tree(Z, n_clusters=k).ravel()
            
            # Check within-cluster distances meeting threshold
            within_total = 0
            within_below = 0

            for i, j, d in pair_distances:
                if cluster_ids[i] == cluster_ids[j]:
                    within_total += 1
                    if d < threshold:
                        within_below += 1

            pct_meeting = (within_below / within_total * 100.0) if within_total > 0 else 100.0

            if pct_meeting >= self.stringency:
                correct_k = k
                best_cluster_df = pd.DataFrame({
                    "Seq_ID": samples,
                    "Assigned_cluster": cluster_ids + 1  # 1-indexed cluster IDs
                })
                print(f"[ClusterFinder] Optimal cluster count found: {correct_k} clusters ({pct_meeting:.2f}% pairs meeting threshold).")
                break
            else:
                print(f"[ClusterFinder] k={k} clusters too small ({pct_meeting:.2f}% < {self.stringency}% meeting threshold).")

        if best_cluster_df is None:
            cluster_ids = cut_tree(Z, n_clusters=k_max).ravel()
            best_cluster_df = pd.DataFrame({
                "Seq_ID": samples,
                "Assigned_cluster": cluster_ids + 1
            })

        # Save output
        if output_dir is None:
            output_dir = "clusters_detected"
        os.makedirs(output_dir, exist_ok=True)
        today_str = datetime.date.today().strftime("%Y-%m-%d")
        output_path = os.path.join(output_dir, f"{today_str}_RESULTING_CLUSTERS_{correct_k}.txt")

        best_cluster_df.to_csv(output_path, sep="\t", index=False)
        print(f"[ClusterFinder] Saved outbreak cluster assignments to: {output_path}")

        return best_cluster_df, correct_k, threshold
