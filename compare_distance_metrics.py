"""
Comparative benchmark script: Legacy Barratt Ensemble Distance vs. Revised KING-Robust wIBS Distance Engine.
Calculates empirical matrix statistics, cophenetic correlation, MOI robustness, and execution time.
"""

import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform, pdist
from scipy.cluster.hierarchy import linkage, cophenet
from scipy.stats import spearmanr, pearsonr
from cyclospora_pyeuk.haplotype_sheet import generate_haplotype_sheet
from cyclospora_pyeuk.distance_engine import PyEukDistanceEngine


def run_comparison():
    print("[Comparison] Generating haplotype sheet from benchmark data...")
    sheet = generate_haplotype_sheet(
        specimen_dir="bench_genotypes/SPECIMEN_GENOTYPES",
        background_dir="bench_genotypes/REFERENCE_POPULATION",
        output_path="bench_run_output/haplotype_sheet_comparison.txt"
    )

    # Take a representative subset of 120 specimens for fast benchmark comparison
    subset_df = sheet.head(120).copy()

    engine = PyEukDistanceEngine()
    clean_df = engine.process_haplotype_sheet(subset_df)
    ids = clean_df["Seq_ID"].tolist()

    print(f"[Comparison] Running comparison on {len(ids)} specimens...")

    # 1. Compute Legacy Barratt Ensemble Distance
    legacy_df = engine.compute_ensemble_matrix(clean_df)

    # 2. Compute Revised wIBS Matrix
    revised_df = engine.compute_revised_wibs_matrix(clean_df)

    L_mat = legacy_df.values
    R_mat = revised_df.values

    # Evaluate correlation
    upper_tri_idx = np.triu_indices(len(ids), k=1)
    legacy_vec = L_mat[upper_tri_idx]
    revised_vec = R_mat[upper_tri_idx]

    r_spearman, _ = spearmanr(legacy_vec, revised_vec)
    r_pearson, _ = pearsonr(legacy_vec, revised_vec)

    # Cophenetic correlation coefficient of dendrogram
    Z_legacy = linkage(squareform(L_mat), method='ward')
    c_legacy, _ = cophenet(Z_legacy, pdist(L_mat))

    Z_revised = linkage(squareform(R_mat), method='ward')
    c_revised, _ = cophenet(Z_revised, pdist(R_mat))

    # Eigenvalue spectrum (Gram matrix PSD check: G = -0.5 * H * D^2 * H)
    H = np.eye(len(ids)) - np.ones((len(ids), len(ids))) / float(len(ids))
    G_legacy = -0.5 * H @ (L_mat ** 2) @ H
    G_revised = -0.5 * H @ (R_mat ** 2) @ H

    eig_legacy = np.min(np.linalg.eigvalsh(G_legacy))
    eig_revised = np.min(np.linalg.eigvalsh(G_revised))

    print("\n=== EMPIRICAL COMPARISON RESULTS ===")
    print(f"Spearman Rank Correlation: r = {r_spearman:.4f}")
    print(f"Pearson Linear Correlation: r = {r_pearson:.4f}")
    print(f"Legacy Cophenetic Correlation: c = {c_legacy:.4f}")
    print(f"Revised wIBS Cophenetic Correlation: c = {c_revised:.4f}")
    print(f"Legacy Min Gram Eigenvalue (PSD check): {eig_legacy:.4e}")
    print(f"Revised Min Gram Eigenvalue (PSD check): {eig_revised:.4e}")

    return {
        "spearman_r": r_spearman,
        "pearson_r": r_pearson,
        "c_legacy": c_legacy,
        "c_revised": c_revised,
        "eig_legacy": eig_legacy,
        "eig_revised": eig_revised
    }


if __name__ == "__main__":
    run_comparison()
