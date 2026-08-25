"""
Modernized CDC Cyclospora cayetanensis MLST typing, Eukaryotyping ensemble distance engine, and revised wIBS distance matrix.
Re-engineered from original CDC High Sierra ALPHA release.
"""

import os
import re
import numpy as np
import pandas as pd
from typing import Tuple, List, Dict, Optional, Union
from numba import njit, prange


def parse_locus_name(col: str) -> str:
    r"""
    Extracts the homologous locus window name from a marker or haplotype identifier.
    Correctly handles:
    - CDC format: Nu_378_PART_A_Hap_4 -> Nu_378_PART_A
    - De novo format: gp60_L752bp.H01_9180 -> gp60
    - Haplotype delimiters: _Hap_\d+, .H\d+, _H\d+, _NOVEL_\d+, .X_\d+
    - Preserves locus names with underscores like Cp_HSP70, Cp_HSP90
    - Strips amplicon length tags like _L245bp
    """
    sub = str(col).strip()
    patterns = [
        r"_Hap_\d+",
        r"\.H\d+(_[A-F0-9]+)?",
        r"_H\d+(_[A-F0-9]+)?",
        r"_(NOVEL|novel)_\d+",
        r"\.X_\d+",
    ]
    for pat in patterns:
        sub = re.sub(pat, "", sub)
    sub = re.sub(r"_L\d+bp$", "", sub)
    sub = sub.rstrip("_")
    if "Junction" in sub or ("Mt_" in sub and "Cmt" in sub):
        return "Mt_Cmt"
    return sub


@njit(parallel=True, fastmath=True)
def _fast_numba_wibs(
    X: np.ndarray,
    w_j: np.ndarray,
    locus_mask: np.ndarray,
    locus_ploidy: np.ndarray
) -> np.ndarray:
    """
    Numba JIT C-kernel for KING-robust Weighted Identity-By-State (wIBS) pairwise distance computation
    with pairwise-complete locus dropout handling and locus-specific ploidy normalization.

    Parameters
    ----------
    X : np.ndarray (n_samples, n_markers)
        Binary presence/absence matrix (1.0 for present, 0.0 for absent/uncalled).
    w_j : np.ndarray (n_markers,)
        KING standardization weight 1 / sqrt(p_j * (1 - p_j)).
    locus_mask : np.ndarray (n_samples, n_markers)
        Boolean mask indicating if the marker belongs to a locus window called in specimen i.
    locus_ploidy : np.ndarray (n_markers,)
        Base ploidy per marker locus (e.g. 1.0 for haploid, 2.0 for diploid).
    """
    nids, n_markers = X.shape
    D_wibs = np.zeros((nids, nids))

    for i in prange(nids):
        for j in range(i + 1, nids):
            # Pairwise complete called markers mask
            valid_markers = locus_mask[i] & locus_mask[j]
            weight_sum = 0.0
            weighted_diff = 0.0
            has_valid_locus = False

            for k in range(n_markers):
                if valid_markers[k]:
                    has_valid_locus = True
                    wk = w_j[k]
                    weight_sum += wk
                    diff = abs(X[i, k] - X[j, k]) / locus_ploidy[k]
                    weighted_diff += diff * wk

            if not has_valid_locus:
                dist = 1.0  # Disjoint locus coverage (0 shared loci) receives maximum distance
            elif weight_sum > 0.0:
                dist = weighted_diff / weight_sum
            else:
                dist = 0.0  # Shared loci are all monomorphic in the cohort

            D_wibs[i, j] = dist
            D_wibs[j, i] = dist

    return D_wibs


class PyEukDistanceEngine:
    """
    Modernized Distance Engine implementing Barratt Heuristic, Plucinski Bayesian LLR,
    and Heterozygosity-Weighted Identity-By-State (wIBS) with exact Gram matrix PSD projection.
    """

    def __init__(
        self,
        epsilon: float = 0.3072,
        min_completeness: float = 0.10,
        ploidy: Optional[Union[int, Dict[str, int]]] = None,
        weight_mode: str = "heterozygosity",
        min_maf: float = 0.0,
        project_psd: bool = True
    ):
        self.epsilon = epsilon
        self.min_completeness = min_completeness
        self.ploidy = ploidy
        self.weight_mode = weight_mode
        self.min_maf = min_maf
        self.project_psd = project_psd

    def process_haplotype_sheet(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Cleans and filters haplotype data sheet, enforcing minimum locus completeness.
        """
        if "Seq_ID" not in df.columns:
            df = df.reset_index()

        marker_cols = [c for c in df.columns if c != "Seq_ID"]

        # Group marker columns by locus window
        locus_to_cols = {}
        for col in marker_cols:
            loc = parse_locus_name(col)
            locus_to_cols.setdefault(loc, []).append(col)

        nloci = len(locus_to_cols)
        if nloci > 0:
            # Specimen has a locus called if at least one allele column in that locus has 'X'
            locus_called_counts = np.zeros(len(df), dtype=int)
            for loc, cols in locus_to_cols.items():
                loc_has_call = (df[cols] == "X").any(axis=1).values
                locus_called_counts += loc_has_call.astype(int)
            completeness = locus_called_counts / float(nloci)
        else:
            completeness = np.ones(len(df))

        eligible_mask = completeness >= self.min_completeness
        cleandata = df[eligible_mask].copy()

        print(f"[PyEuk] Filtered dataset: {len(cleandata)} / {len(df)} specimens passed completeness criteria (min_locus_completeness >= {self.min_completeness}).")
        return cleandata

    def compute_revised_wibs_matrix(
        self,
        df: pd.DataFrame,
        weight_mode: Optional[str] = None,
        min_maf: Optional[float] = None,
        project_psd: Optional[bool] = None
    ) -> pd.DataFrame:
        """
        Computes Revised Weighted Identity-By-State (wIBS) Distance Matrix
        with pairwise-complete locus dropout handling, configurable allele weighting,
        and exact Gram matrix PSD projection.

        Parameters
        ----------
        df : pd.DataFrame
            Haplotype binary presence/absence indicator sheet.
        weight_mode : str, optional
            Allele weighting scheme:
            - 'heterozygosity' / 'het' (default): w = 2*p*(1-p) (expected heterozygosity / binomial variance for presence/absence indicators).
            - 'inverted_king' / 'inv_king': w = sqrt(p*(1-p)) (binomial standard deviation).
            - 'king': w = 1 / sqrt(p*(1-p)) (KING standardization for centered genotype dosage).
            - 'uniform' / 'unweighted': w = 1.0 (unweighted polymorphic matching).
        min_maf : float, optional
            Minor allele frequency filter threshold (default: 0.0). Columns with p < min_maf or p > 1 - min_maf receive weight 0.
        project_psd : bool, optional
            Whether to project the distance matrix onto the Positive Semi-Definite (PSD) cone via Gram double centering (default: True).
        """
        w_mode = (weight_mode or self.weight_mode).lower().replace("-", "_")
        maf_thresh = self.min_maf if min_maf is None else min_maf
        do_psd = self.project_psd if project_psd is None else project_psd

        clean_df = self.process_haplotype_sheet(df)
        ids = clean_df["Seq_ID"].tolist()
        nids = len(ids)
        marker_cols = [c for c in clean_df.columns if c != "Seq_ID"]

        X = (clean_df[marker_cols].values == "X").astype(np.float64)

        # Map marker columns to locus windows
        locus_names = [parse_locus_name(c) for c in marker_cols]

        def get_loc_ploidy(loc_name: str) -> float:
            if self.ploidy is None:
                return 1.0 if loc_name.startswith("Mt") else 2.0
            elif isinstance(self.ploidy, (int, float)):
                return float(self.ploidy)
            elif isinstance(self.ploidy, dict):
                return float(self.ploidy.get(loc_name, 2.0))
            return 2.0

        locus_ploidy = np.array([get_loc_ploidy(name) for name in locus_names], dtype=np.float64)

        # Determine locus call status per specimen (locus_mask: True if specimen has >=1 allele called in locus window)
        locus_mask = np.zeros_like(X, dtype=bool)
        unique_loci = list(set(locus_names))

        for loc in unique_loci:
            loc_indices = [k for k, name in enumerate(locus_names) if name == loc]
            loc_sub = X[:, loc_indices]
            called_specimens = loc_sub.sum(axis=1) > 0
            for idx in loc_indices:
                locus_mask[:, idx] = called_specimens

        # Calculate population allele frequencies p_j over called instances and assign weights
        p_j = np.zeros(X.shape[1], dtype=np.float64)
        w_j = np.zeros(X.shape[1], dtype=np.float64)

        for k in range(X.shape[1]):
            called_count = locus_mask[:, k].sum()
            if called_count > 0:
                p_val = float(X[locus_mask[:, k], k].mean())
                p_j[k] = p_val
                # Informative polymorphic alleles receive weights according to weight_mode
                # Unobserved columns (p == 0.0) and fixed columns (p == 1.0) receive weight 0.0
                if 0.0 < p_val < 1.0:
                    if maf_thresh > 0.0 and (p_val < maf_thresh or p_val > 1.0 - maf_thresh):
                        w_j[k] = 0.0
                    elif w_mode in ["heterozygosity", "het"]:
                        # Expected heterozygosity / binomial variance: 2 * p * (1 - p)
                        w_j[k] = 2.0 * p_val * (1.0 - p_val)
                    elif w_mode in ["inverted_king", "inv_king"]:
                        # Binomial standard deviation: sqrt(p * (1 - p))
                        w_j[k] = np.sqrt(p_val * (1.0 - p_val))
                    elif w_mode in ["uniform", "unweighted"]:
                        w_j[k] = 1.0
                    elif w_mode == "king":
                        # KING standardization for centered genotype dosage: 1 / sqrt(p * (1 - p))
                        p_clip = np.clip(p_val, 1.0 / (2.0 * called_count), 1.0 - 1.0 / (2.0 * called_count))
                        w_j[k] = 1.0 / np.sqrt(p_clip * (1.0 - p_clip))
                    else:
                        raise ValueError(f"Unknown weight_mode: '{w_mode}'. Choose from 'heterozygosity', 'inverted_king', 'king', 'uniform'.")

        # Parallel Numba execution
        print(f"[PyEuk-wIBS] Executing Numba JIT C-kernel on {nids} specimens across {len(unique_loci)} locus windows (weight_mode='{w_mode}', min_maf={maf_thresh})...")
        D_wibs = _fast_numba_wibs(X, w_j, locus_mask, locus_ploidy)

        if not do_psd:
            print("[PyEuk-wIBS] Returning raw wIBS Distance Matrix (PSD projection disabled).")
            return pd.DataFrame(D_wibs, index=ids, columns=ids)

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

    def compute_snp_weighted_wibs_matrix(
        self,
        df: pd.DataFrame,
        sequence_map: Optional[Dict[str, str]] = None,
        fasta_path: Optional[str] = None,
        weight_mode: Optional[str] = None,
        min_maf: Optional[float] = None,
        project_psd: Optional[bool] = None
    ) -> pd.DataFrame:
        """
        Computes Reference-Free SNP-Weighted wIBS Distance Matrix (PyEuk v0.3.0).
        Calculates de novo pairwise sequence alignment Hamming distances directly between called haplotype sequences,
        eliminating reliance on external reference database files.
        """
        import glob
        w_mode = (weight_mode or self.weight_mode).lower().replace("-", "_")
        maf_thresh = self.min_maf if min_maf is None else min_maf
        do_psd = self.project_psd if project_psd is None else project_psd

        clean_df = self.process_haplotype_sheet(df)
        ids = clean_df["Seq_ID"].tolist()
        nids = len(ids)
        marker_cols = [c for c in clean_df.columns if c != "Seq_ID"]

        # Load de novo sequence map or parse input files if provided
        records = {}
        fpaths = []
        if sequence_map:
            records.update(sequence_map)
        elif fasta_path and os.path.exists(fasta_path):
            fpaths = [fasta_path]
        else:
            fpaths = glob.glob("cdc_reference_data/**/*.fasta", recursive=True) + glob.glob("cdc_reference_data/**/*.fa", recursive=True)

        for fp in fpaths:
            if "Illumina" in fp or "MAPPING" in fp:
                continue
            with open(fp) as f:
                cur_header = None
                seq_lines = []
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        if cur_header:
                            records[cur_header] = "".join(seq_lines)
                        cur_header = line[1:].strip()
                        seq_lines = []
                    else:
                        seq_lines.append(line)
                if cur_header:
                    records[cur_header] = "".join(seq_lines)

        def norm_h(name):
            if "." in name:
                name = name.split(".")[-1]
            return name.replace("-", "_").strip()

        fasta_map = {norm_h(k): seq for k, seq in records.items()}

        def seq_dist(seq1, seq2):
            min_len = min(len(seq1), len(seq2))
            mismatches = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a != b)
            mismatches += abs(len(seq1) - len(seq2))
            max_len = max(len(seq1), len(seq2))
            return float(mismatches) / float(max_len) if max_len > 0 else 0.0

        # Group marker columns by locus window
        locus_map = {}
        for col in marker_cols:
            loc = parse_locus_name(col)
            locus_map.setdefault(loc, []).append(col)

        # Calculate locus weights based on population marker allele frequencies
        locus_weights = {}
        for loc, cols in locus_map.items():
            sub = (clean_df[cols].values == "X").astype(np.float64)
            p_j = sub.mean(axis=0)
            w_cols = np.zeros(len(p_j), dtype=np.float64)
            for k, p_val in enumerate(p_j):
                if 0.0 < p_val < 1.0:
                    if maf_thresh > 0.0 and (p_val < maf_thresh or p_val > 1.0 - maf_thresh):
                        w_cols[k] = 0.0
                    elif w_mode in ["heterozygosity", "het"]:
                        w_cols[k] = 2.0 * p_val * (1.0 - p_val)
                    elif w_mode in ["inverted_king", "inv_king"]:
                        w_cols[k] = np.sqrt(p_val * (1.0 - p_val))
                    elif w_mode in ["uniform", "unweighted"]:
                        w_cols[k] = 1.0
                    elif w_mode == "king":
                        p_clip = np.clip(p_val, 1.0 / (2.0 * len(clean_df)), 1.0 - 1.0 / (2.0 * len(clean_df)))
                        w_cols[k] = 1.0 / np.sqrt(p_clip * (1.0 - p_clip))
            locus_weights[loc] = float(np.mean(w_cols)) if len(w_cols) > 0 and np.mean(w_cols) > 0 else 1.0

        dist_mat = np.zeros((nids, nids), dtype=np.float64)
        sample_active = []
        for idx, row in clean_df.iterrows():
            locus_active = {}
            for locus, cols in locus_map.items():
                active_cols = [c for c in cols if row[c] == "X"]
                if active_cols:
                    locus_active[locus] = active_cols
            sample_active.append(locus_active)

        for i in range(nids):
            for j in range(i + 1, nids):
                common_loci = set(sample_active[i].keys()) & set(sample_active[j].keys())
                if not common_loci:
                    dist_mat[i, j] = 1.0  # Disjoint locus coverage receives maximum distance (1.0)
                    dist_mat[j, i] = 1.0
                    continue

                weighted_scores = []
                weight_sum = 0.0
                for loc in common_loci:
                    w_loc = locus_weights[loc]
                    haps_i = sample_active[i][loc]
                    haps_j = sample_active[j][loc]
                    if set(haps_i) & set(haps_j):
                        loc_d = 0.0
                    else:
                        min_d = 1.0
                        for hi in haps_i:
                            seq_i = fasta_map.get(norm_h(hi), "")
                            for hj in haps_j:
                                seq_j = fasta_map.get(norm_h(hj), "")
                                if seq_i and seq_j:
                                     d = seq_dist(seq_i, seq_j)
                                else:
                                    d = 1.0
                                if d < min_d:
                                    min_d = d
                        loc_d = min_d
                    weighted_scores.append(loc_d * w_loc)
                    weight_sum += w_loc

                mean_dist = float(sum(weighted_scores) / weight_sum) if weight_sum > 0.0 else 1.0
                dist_mat[i, j] = mean_dist
                dist_mat[j, i] = mean_dist

        if not do_psd:
            print("[PyEuk-SNP-wIBS] Returning raw SNP-wIBS Distance Matrix (PSD projection disabled).")
            return pd.DataFrame(dist_mat, index=ids, columns=ids)

        # Gram matrix PSD Projection
        H = np.eye(nids) - np.ones((nids, nids)) / float(nids)
        G = -0.5 * H @ (dist_mat ** 2) @ H
        evals, evecs = np.linalg.eigh(G)
        evals_clipped = np.maximum(evals, 0.0)
        G_psd = evecs @ np.diag(evals_clipped) @ evecs.T

        diag_g = np.diag(G_psd)
        D_sq = np.clip(diag_g[:, None] + diag_g[None, :] - 2.0 * G_psd, 0.0, None)
        D_psd = np.sqrt(D_sq)
        D_psd = (D_psd + D_psd.T) / 2.0
        np.fill_diagonal(D_psd, 0.0)

        print("[PyEuk-SNP-wIBS] SNP-Weighted wIBS Distance Matrix successfully computed (PSD Guaranteed: lambda_min >= 0.0).")
        return pd.DataFrame(D_psd, index=ids, columns=ids)

    def compute_ensemble_matrix(
        self,
        df: pd.DataFrame,
        project_psd: Optional[bool] = None
    ) -> pd.DataFrame:
        """
        Computes the dual-ensemble matrix (Barratt Heuristic + Plucinski Bayesian LLR Rank Integration)
        with exact Gram matrix PSD projection.
        """
        do_psd = self.project_psd if project_psd is None else project_psd
        clean_df = self.process_haplotype_sheet(df)
        ids, loci_data, ploidy = self._extract_locus_data(clean_df)
        nids = len(ids)

        print(f"[PyEuk] Computing Barratt's Heuristic Distance...")
        D_heur = self.compute_heuristic_distance(ids, loci_data)

        print(f"[PyEuk] Computing Plucinski's Bayesian Distance...")
        D_bayes = self.compute_bayesian_distance(ids, loci_data, ploidy)

        print(f"[PyEuk] Performing Uniform Fractional Rank Integration...")
        D_ensemble = self.rank_integrate(D_heur, D_bayes)

        # Exact Positive Semi-Definite (PSD) projection via Gram Matrix double centering
        if do_psd and nids > 2:
            H = np.eye(nids) - np.ones((nids, nids)) / float(nids)
            G = -0.5 * H @ (D_ensemble ** 2) @ H
            evals, evecs = np.linalg.eigh(G)
            evals_clipped = np.maximum(evals, 0.0)
            G_psd = evecs @ np.diag(evals_clipped) @ evecs.T
            diag_g = np.diag(G_psd)
            D_sq = np.clip(diag_g[:, None] + diag_g[None, :] - 2.0 * G_psd, 0.0, None)
            D_psd = np.sqrt(D_sq)
            D_psd = (D_psd + D_psd.T) / 2.0
            np.fill_diagonal(D_psd, 0.0)
            D_ensemble = D_psd

        return pd.DataFrame(D_ensemble, index=ids, columns=ids)

    def _extract_locus_data(self, df: pd.DataFrame) -> Tuple[List[str], List[Dict], np.ndarray]:
        ids = df["Seq_ID"].tolist()
        nids = len(ids)
        marker_cols = [c for c in df.columns if c != "Seq_ID"]

        locus_names = []
        for col in marker_cols:
            loc = parse_locus_name(col)
            if loc not in locus_names:
                locus_names.append(loc)

        nloci = len(locus_names)
        if self.ploidy is None:
            ploidy = np.array([1 if loc.startswith("Mt") else 2 for loc in locus_names])
        elif isinstance(self.ploidy, (int, float)):
            ploidy = np.full(nloci, int(self.ploidy))
        elif isinstance(self.ploidy, dict):
            ploidy = np.array([self.ploidy.get(loc, 2) for loc in locus_names])
        else:
            ploidy = np.array([1 if loc.startswith("Mt") else 2 for loc in locus_names])

        loci_data = []
        for j, loc in enumerate(locus_names):
            loc_cols = [c for c in marker_cols if parse_locus_name(c) == loc]
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
        max_d = np.nanmax(mean_dist) if not np.isnan(np.nanmax(mean_dist)) and np.nanmax(mean_dist) > 0 else 1.0
        mean_dist = np.nan_to_num(mean_dist, nan=max_d)
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
