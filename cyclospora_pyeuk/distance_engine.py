"""
Modernized CDC Cyclospora cayetanensis MLST typing, Eukaryotyping ensemble distance engine, and revised wIBS distance matrix.
Re-engineered from original CDC High Sierra ALPHA release.
"""

import os
import numpy as np
import pandas as pd
from typing import Tuple, List, Dict, Optional
from numba import njit, prange


@njit(parallel=True, fastmath=True)
def _fast_numba_wibs(X: np.ndarray, w_j: np.ndarray, locus_mask: np.ndarray) -> np.ndarray:
    """
    Numba JIT C-kernel for KING-robust Weighted Identity-By-State (wIBS) pairwise distance computation
    with pairwise-complete locus dropout handling.

    Parameters
    ----------
    X : np.ndarray (n_samples, n_markers)
        Binary presence/absence matrix (1.0 for present, 0.0 for absent/uncalled).
    w_j : np.ndarray (n_markers,)
        KING standardization weight 1 / sqrt(p_j * (1 - p_j)).
    locus_mask : np.ndarray (n_samples, n_markers)
        Boolean mask indicating if the marker belongs to a locus window called in specimen i.
    """
    nids, n_markers = X.shape
    D_wibs = np.zeros((nids, nids))

    for i in prange(nids):
        for j in range(i + 1, nids):
            # Pairwise complete called markers mask
            valid_markers = locus_mask[i] & locus_mask[j]
            weight_sum = 0.0
            weighted_diff = 0.0

            for k in range(n_markers):
                if valid_markers[k]:
                    wk = w_j[k]
                    weight_sum += wk
                    diff = abs(X[i, k] - X[j, k])
                    weighted_diff += diff * wk

            if weight_sum > 0.0:
                dist = weighted_diff / weight_sum
            else:
                dist = 0.0

            D_wibs[i, j] = dist
            D_wibs[j, i] = dist

    return D_wibs


class PyEukDistanceEngine:
    """
    Modernized Distance Engine implementing Barratt Heuristic, Plucinski Bayesian LLR,
    and KING-Robust Weighted Identity-By-State (wIBS) with exact Gram matrix PSD projection.
    """

    def __init__(self, epsilon: float = 0.3072, min_completeness: float = 0.10):
        self.epsilon = epsilon
        self.min_completeness = min_completeness

    def process_haplotype_sheet(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Cleans and filters haplotype data sheet, enforcing minimum locus completeness.
        """
        if "Seq_ID" not in df.columns:
            df = df.reset_index()

        marker_cols = [c for c in df.columns if c != "Seq_ID"]

        # Calculate specimen completeness (fraction of non-empty calls)
        specimen_calls = (df[marker_cols] == "X").sum(axis=1)
        max_possible = len(marker_cols)
        completeness = specimen_calls / max_possible

        eligible_mask = completeness >= self.min_completeness
        cleandata = df[eligible_mask].copy()

        print(f"[PyEuk] Filtered dataset: {len(cleandata)} / {len(df)} specimens passed completeness criteria.")
        return cleandata

    def compute_revised_wibs_matrix(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Computes Revised KING-Robust Weighted Identity-By-State (wIBS) Distance Matrix
        with pairwise-complete locus dropout handling and exact Gram matrix PSD projection.
        """
        clean_df = self.process_haplotype_sheet(df)
        ids = clean_df["Seq_ID"].tolist()
        nids = len(ids)
        marker_cols = [c for c in clean_df.columns if c != "Seq_ID"]

        X = (clean_df[marker_cols].values == "X").astype(np.float64)

        # Map marker columns to locus windows
        def get_locus_name(col):
            sub = col.split("_Hap_")[0]
            if "Junction" in sub or ("Mt_" in sub and "Cmt" in sub):
                return "Mt_Cmt"
            return sub

        locus_names = [get_locus_name(c) for c in marker_cols]

        # Determine locus call status per specimen (locus_mask: True if specimen has >=1 allele called in locus window)
        locus_mask = np.zeros_like(X, dtype=bool)
        unique_loci = list(set(locus_names))

        for loc in unique_loci:
            loc_indices = [k for k, name in enumerate(locus_names) if name == loc]
            loc_sub = X[:, loc_indices]
            called_specimens = loc_sub.sum(axis=1) > 0
            for idx in loc_indices:
                locus_mask[:, idx] = called_specimens

        # Calculate population allele frequencies p_j over called instances
        p_j = np.zeros(X.shape[1])
        for k in range(X.shape[1]):
            called_count = locus_mask[:, k].sum()
            if called_count > 0:
                p_j[k] = X[locus_mask[:, k], k].mean()
            else:
                p_j[k] = 1e-4

        p_j = np.clip(p_j, 1e-4, 1.0 - 1e-4)

        # KING-robust weight w_j = 1 / sqrt(p_j * (1 - p_j))
        w_j = 1.0 / np.sqrt(p_j * (1.0 - p_j))

        # Parallel Numba execution
        print(f"[PyEuk-wIBS] Executing Numba JIT C-kernel on {nids} specimens across {len(unique_loci)} locus windows...")
        D_wibs = _fast_numba_wibs(X, w_j, locus_mask)

        # Exact Positive Semi-Definite (PSD) projection via Gram Matrix double centering
        H = np.eye(nids) - np.ones((nids, nids)) / float(nids)
        G = -0.5 * H @ (D_wibs ** 2) @ H

        evals, evecs = np.linalg.eigh(G)
        evals_clipped = np.maximum(evals, 0.0)
        G_psd = evecs @ np.diag(evals_clipped) @ evecs.T

        diag_g = np.diag(G_psd)
        D_sq = np.clip(diag_g[:, None] + diag_g[None, :] - 2.0 * G_psd, 0.0, None)
        D_psd = np.sqrt(D_sq)
        D_psd = (D_psd + D_psd.T) / 2.0
        np.fill_diagonal(D_psd, 0.0)

        print("[PyEuk-wIBS] Revised wIBS Distance Matrix successfully computed (PSD Guaranteed: lambda_min >= 0.0).")
        return pd.DataFrame(D_psd, index=ids, columns=ids)

    def compute_ensemble_matrix(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Computes the dual-ensemble matrix (Barratt Heuristic + Plucinski Bayesian LLR Rank Integration).
        """
        clean_df = self.process_haplotype_sheet(df)
        ids, loci_data, ploidy = self._extract_locus_data(clean_df)

        print(f"[PyEuk] Computing Barratt's Heuristic Distance...")
        D_heur = self.compute_heuristic_distance(ids, loci_data)

        print(f"[PyEuk] Computing Plucinski's Bayesian Distance...")
        D_bayes = self.compute_bayesian_distance(ids, loci_data, ploidy)

        print(f"[PyEuk] Performing Uniform Fractional Rank Integration...")
        D_ensemble = self.rank_integrate(D_heur, D_bayes)

        return pd.DataFrame(D_ensemble, index=ids, columns=ids)

    def _extract_locus_data(self, df: pd.DataFrame) -> Tuple[List[str], List[Dict], np.ndarray]:
        ids = df["Seq_ID"].tolist()
        nids = len(ids)
        marker_cols = [c for c in df.columns if c != "Seq_ID"]

        def get_locus_name(col):
            sub = col.split("_Hap_")[0]
            if "Junction" in sub or ("Mt_" in sub and "Cmt" in sub):
                return "Mt_Cmt"
            return sub

        locus_names = []
        for col in marker_cols:
            loc = get_locus_name(col)
            if loc not in locus_names:
                locus_names.append(loc)

        nloci = len(locus_names)
        ploidy = np.array([1 if loc.startswith("Mt") else 2 for loc in locus_names])

        loci_data = []
        for j, loc in enumerate(locus_names):
            loc_cols = [c for c in marker_cols if get_locus_name(c) == loc]
            sub = df[loc_cols].values == "X"
            alleles = loc_cols
            nalleles = len(alleles)

            counts = sub.sum(axis=0)
            total_obs = counts.sum()
            freqs = counts / total_obs if total_obs > 0 else np.zeros(nalleles)

            specimen_alleles = []
            for i in range(nids):
                present = [alleles[k] for k in range(nalleles) if sub[i, k]]
                specimen_alleles.append(present)

            loci_data.append({
                "locus": loc,
                "ploidy": ploidy[j],
                "cols": loc_cols,
                "alleles": alleles,
                "freqs": freqs,
                "freq_map": dict(zip(alleles, freqs)),
                "specimen_alleles": specimen_alleles,
                "presence_matrix": sub
            })

        return ids, loci_data, ploidy

    def compute_heuristic_distance(self, ids: List[str], loci_data: List[Dict]) -> np.ndarray:
        nids = len(ids)
        nloci = len(loci_data)

        H_nu = np.zeros(nloci)
        for j in range(nloci):
            p = loci_data[j]["freqs"]
            p_pos = p[p > 0]
            H_nu[j] = -np.sum(p_pos * np.log2(p_pos)) if len(p_pos) > 0 else 0.0

        locus_dists = np.full((nloci, nids, nids), np.nan)

        for j in range(nloci):
            loc_info = loci_data[j]
            ploidy = loc_info["ploidy"]
            spec_alleles = loc_info["specimen_alleles"]
            h_val = H_nu[j]

            for i1 in range(nids):
                v1 = spec_alleles[i1]
                if len(v1) == 0:
                    continue
                for i2 in range(i1, nids):
                    v2 = spec_alleles[i2]
                    if len(v2) == 0:
                        continue

                    x = len(set(v1) | set(v2))
                    n_min = min(len(v1), len(v2))
                    shared = set(v1) & set(v2)
                    y = len(shared)

                    if ploidy > 1:
                        w = x * (n_min > 1) + 4 * (n_min == 1 and x == 2) + (1 + x) * (n_min == 1 and x > 2)
                        z = 3 * (y == 1 and n_min == 1 and x > 2) + 2 * (y > 0 and n_min > 1)
                        delta_raw = w if y == 0 else (4.0 - z)
                    else:
                        w = x * (n_min > 1) + (1 + x) * (n_min == 1 and x > 1)
                        z = 3 * (y == 1 and n_min == 1 and x > 1) + 2 * (y > 0 and n_min > 1)
                        delta_raw = w if y == 0 else (4.0 - z)

                    d_final = h_val * delta_raw
                    locus_dists[j, i1, i2] = d_final
                    locus_dists[j, i2, i1] = d_final

        with np.errstate(all="ignore"):
            mean_dist = np.nanmean(locus_dists, axis=0)
        np.fill_diagonal(mean_dist, 0.0)
        return mean_dist

    def compute_bayesian_distance(self, ids: List[str], loci_data: List[Dict], ploidy: np.ndarray) -> np.ndarray:
        nids = len(ids)
        nloci = len(loci_data)
        D_bayes = np.zeros((nids, nids))

        for i1 in range(nids):
            for i2 in range(i1 + 1, nids):
                L0, L1, L2 = 0.0, 0.0, 0.0
                valid_loci = 0

                for j in range(nloci):
                    v1 = loci_data[j]["specimen_alleles"][i1]
                    v2 = loci_data[j]["specimen_alleles"][i2]
                    fmap = loci_data[j]["freq_map"]

                    if len(v1) == 0 or len(v2) == 0:
                        continue

                    valid_loci += 1
                    p1 = [max(fmap.get(a, 1e-4), 1e-4) for a in v1]
                    p2 = [max(fmap.get(a, 1e-4), 1e-4) for a in v2]

                    l0 = np.sum(np.log(p1)) + np.sum(np.log(p2))

                    shared = list(set(v1) & set(v2))
                    if len(shared) > 0:
                        l1_vals = []
                        for a in shared:
                            idx1 = v1.index(a)
                            p1_sub = [p1[k] for k in range(len(p1)) if k != idx1]
                            l1_vals.append(np.sum(np.log(p1_sub)) + np.sum(np.log(p2)))
                        l1 = max(l1_vals)
                    else:
                        l1 = np.log(self.epsilon * min(p1 + p2))

                    if ploidy[j] > 1 and len(v1) > 1 and len(v2) > 1:
                        import itertools
                        pairs1 = list(itertools.combinations(v1, 2))
                        pairs2 = list(itertools.combinations(v2, 2))
                        l2_vals = []
                        for pair1 in pairs1:
                            for pair2 in pairs2:
                                if sorted(pair1) == sorted(pair2):
                                    idx_sub = [k for k, a in enumerate(v1) if a not in pair1]
                                    p1_sub = [p1[k] for k in idx_sub]
                                    l2_vals.append(np.sum(np.log(p1_sub)) + np.sum(np.log(p2)))
                        l2 = max(l2_vals) if len(l2_vals) > 0 else np.log(self.epsilon * min(p1 + p2))
                    else:
                        l2 = l1

                    L0 += l0
                    L1 += l1
                    L2 += l2

                if valid_loci > 0:
                    L_arr = np.array([L0, L1, L2])
                    L_max = np.max(L_arr)
                    exp_L = np.exp(L_arr - L_max)
                    P = exp_L / np.sum(exp_L)
                    d_b = 0.0 * P[2] + 1.0 * P[1] + 2.0 * P[0]
                else:
                    d_b = 2.0

                D_bayes[i1, i2] = d_b
                D_bayes[i2, i1] = d_b

        return D_bayes

    def rank_integrate(self, D1: np.ndarray, D2: np.ndarray) -> np.ndarray:
        from scipy.stats import rankdata
        nids = D1.shape[0]

        triu_idx = np.triu_indices(nids, k=1)
        v1 = D1[triu_idx]
        v2 = D2[triu_idx]

        r1 = rankdata(v1, method="average") / len(v1)
        r2 = rankdata(v2, method="average") / len(v2)

        r_ens = 0.5 * r1 + 0.5 * r2

        D_ens = np.zeros((nids, nids))
        D_ens[triu_idx] = r_ens
        D_ens = D_ens + D_ens.T
        return D_ens
