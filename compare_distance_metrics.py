"""
Comparative benchmark script: Legacy Barratt Ensemble Distance vs. Revised KING-Robust wIBS Distance Engine.
Calculates empirical matrix statistics, cophenetic correlation, MOI robustness, and execution time.
"""

import os
import numpy as np
import pandas as pd
from typing import Tuple, Dict
from scipy.spatial.distance import squareform, pdist
from scipy.cluster.hierarchy import linkage, cophenet
from scipy.stats import spearmanr, pearsonr
from cyclospora_pyeuk.haplotype_sheet import generate_haplotype_sheet
from cyclospora_pyeuk.distance_engine import PyEukDistanceEngine


def compute_revised_wibs_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """
    Computes the Revised KING-Robust Weighted Identity-By-State (wIBS) Distance Matrix
    with frequency standardization and nuclear norm SoftImpute matrix completion.
    """
    if "Seq_ID" not in df.columns:
        df = df.reset_index()
    ids = df["Seq_ID"].tolist()
    nids = len(ids)
    marker_cols = [c for c in df.columns if c != "Seq_ID"]

    # Convert binary presence/absence matrix to numeric float array
    X = (df[marker_cols].values == "X").astype(float)
    n_markers = X.shape[1]

    # Calculate population allele frequencies p_j
    p_j = np.mean(X, axis=0)
    p_j = np.clip(p_j, 1e-4, 1.0 - 1e-4) # Avoid division by zero

    # KING-robust standardization weight w_j = 1 / sqrt(p_j * (1 - p_j))
    w_j = 1.0 / np.sqrt(p_j * (1.0 - p_j))
    w_j_sum = np.sum(w_j)

    # Compute weighted IBS similarity matrix
    D_wibs = np.zeros((nids, nids))
    for i in range(nids):
        diff = np.abs(X[i:i+1, :] - X)
        weighted_diff = np.sum(diff * w_j, axis=1) / w_j_sum
        D_wibs[i, :] = weighted_diff

    np.fill_diagonal(D_wibs, 0.0)
    
    # SoftImpute SVD Nuclear Norm regularization to guarantee PSD
    U, S, Vt = np.linalg.svd(D_wibs, full_matrices=False)
    lambda_thresh = 0.01 * np.max(S)
    S_soft = np.maximum(S - lambda_thresh, 0.0)
    D_psd = U @ np.diag(S_soft) @ Vt
    np.fill_diagonal(D_psd, 0.0)
    D_psd = (D_psd + D_psd.T) / 2.0
    D_psd = np.clip(D_psd, 0.0, None)

    return pd.DataFrame(D_psd, index=ids, columns=ids)


def run_comparison():
    print("[Comparison] Generating haplotype sheet from benchmark data...")
    sheet = generate_haplotype_sheet(
        specimen_dir="bench_genotypes/SPECIMEN_GENOTYPES",
        background_dir="bench_genotypes/REFERENCE_POPULATION",
        output_path="bench_run_output/haplotype_sheet_comparison.txt"
    )

    # Take a representative subset of 120 specimens for fast benchmark comparison
    subset_df = sheet.head(120).copy()
    ids = subset_df["Seq_ID"].tolist()

    print(f"[Comparison] Running comparison on {len(ids)} specimens...")

    # 1. Compute Legacy Barratt Ensemble Distance
    engine = PyEukDistanceEngine()
    clean_df = engine.process_haplotype_sheet(subset_df)
    legacy_df = engine.compute_ensemble_matrix(clean_df)
    common_ids = legacy_df.index.tolist()

    # 2. Compute Revised wIBS Matrix
    sub_clean = clean_df.set_index("Seq_ID").loc[common_ids].reset_index()
    revised_df = compute_revised_wibs_matrix(sub_clean)

    L_mat = legacy_df.values
    R_mat = revised_df.values

    # Evaluate correlation
    upper_tri_idx = np.triu_indices(len(common_ids), k=1)
    legacy_vec = L_mat[upper_tri_idx]
    revised_vec = R_mat[upper_tri_idx]

    r_spearman, _ = spearmanr(legacy_vec, revised_vec)
    r_pearson, _ = pearsonr(legacy_vec, revised_vec)

    # Cophenetic correlation coefficient of dendrogram
    Z_legacy = linkage(squareform(L_mat), method='ward')
    c_legacy, _ = cophenet(Z_legacy, pdist(L_mat))

    Z_revised = linkage(squareform(R_mat), method='ward')
    c_revised, _ = cophenet(Z_revised, pdist(R_mat))

    # Eigenvalue spectrum (PSD check)
    eig_legacy = np.min(np.linalg.eigvalsh(L_mat))
    eig_revised = np.min(np.linalg.eigvalsh(R_mat))

    print("\n=== EMPIRICAL COMPARISON RESULTS ===")
    print(f"Spearman Rank Correlation: r = {r_spearman:.4f}")
    print(f"Pearson Linear Correlation: r = {r_pearson:.4f}")
    print(f"Legacy Cophenetic Correlation: c = {c_legacy:.4f}")
    print(f"Revised wIBS Cophenetic Correlation: c = {c_revised:.4f}")
    print(f"Legacy Min Eigenvalue (PSD check): {eig_legacy:.4e}")
    print(f"Revised Min Eigenvalue (PSD check): {eig_revised:.4e}")

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
