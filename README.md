# PyEuk: Modern *Cyclospora* Genotyping & Outbreak Detection

[![Version](https://img.shields.io/badge/version-0.3.0-blue.svg)](https://github.com/spond/pyeuk)
[![Python](https://img.shields.io/badge/python-3.8%2B-green.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-Apache%202.0-orange.svg)](LICENSE)
[![Speedup](https://img.shields.io/badge/speedup-99.2x-brightgreen.svg)]()

**PyEuk** is a high-performance Python framework for molecular typing, genetic distance estimation, and foodborne outbreak cluster detection in the human apicomplexan parasite ***Cyclospora cayetanensis***. 

It replaces legacy, brittle heuristics with a fast, mathematically rigorous distance engine and automated, label-free hierarchical clustering.

---

## 🚀 Core Driver Features

### 1. Flexible Haplotype Ingestion & De Novo Discovery
* **External Assembly Ingestion (`-a / --assembled-fasta`)**: Directly ingest assembled FASTA contigs from SPAdes, Flye, MEGAHIT, or Galaxy pipelines without manual BLAST parsing.
* **Reference-Free De Novo Discovery (`--de-novo`)**: Discover homologous loci and phased haplotypes directly from sequence contigs without requiring pre-existing reference databases.
* **Deterministic Naming Scheme**: Mints content-addressable identifiers (`<Locus>_L<Length>bp.H<Rank>_<Hash4>`, e.g., `Nu_378_PART_A_L245bp.H01_508B`) embedding locus anchor, amplicon length, cohort frequency rank, and an MD5 sequence hash for global cross-lab reproducibility.

### 2. Dropout-Robust Genetic Distance Engine
* **KING-Weighted Identity-by-State (wIBS)**: Evaluates pairwise genetic dissimilarity across multi-locus marker panels, properly weighting population allele frequencies and handling co-infections.
* **Pairwise-Complete Dropout Tolerance**: Dissimilarity is computed only over mutually amplified loci, preventing PCR sequencing dropouts from triggering artificial distance spikes.
* **Gram Matrix PSD Projection**: Guarantees positive semi-definite Euclidean metric geometry ($\lambda_{\min} \ge 0.0$) for valid, mathematically sound Ward hierarchical clustering.
* **Vectorized Acceleration**: Accelerated via vectorized NumPy and Numba JIT kernels (**99.2× faster** than legacy R scripts, processing 1,000+ specimens in seconds).

### 3. Automated Label-Free Outbreak Clustering
* **Unsupervised Knee Detection**: Automatically determines the optimal number of outbreak clusters ($k$) via scale-free relative merge-height gap analysis (`rel_gap ≥ 0.2200`), eliminating the need for labeled training data or manual cutoffs.
* **100% Deterministic Tree Cuts**: Uses lexicographical tie-breaking to eliminate non-deterministic clustering artifacts across runs.
* **Outlier & Noise Guards**: Incorporates cohort size guards to prevent false cluster splitting on background surveillance samples.

---

## 📦 Installation

```bash
git clone https://github.com/spond/pyeuk.git
cd pyeuk
pip install -e .
```

**Requirements:** Python ≥ 3.8, `numpy`, `scipy`, `pandas`, `scikit-learn`, `numba`.

---

## ⚡ Quickstart

### 1. Run Complete Pipeline on Genotype Calls
```bash
# Process specimen call files or .zip archives and detect outbreak clusters
cyclospora-typing run-all \
    -s ./specimens_dir_or_zip \
    -o ./outbreak_results
```

### 2. Ingest Assembled FASTA Contigs (Reference-Free De Novo)
```bash
# Ingest external FASTA assemblies and discover haplotypes de novo
cyclospora-typing run-all \
    -a ./assembled_contigs.fasta \
    --de-novo \
    -o ./de_novo_results
```

### 3. Modular Step-by-Step CLI Commands
```bash
# Step 1: Generate binary presence/absence haplotype sheet
cyclospora-typing generate-sheet -s ./specimens -o haplotype_sheet.txt

# Step 2: Compute pairwise KING-wIBS distance matrix
cyclospora-typing eukaryotyping -i haplotype_sheet.txt -o distance_matrix.csv --wibs

# Step 3: Run prospective outbreak clustering
cyclospora-typing cluster -m distance_matrix.csv -o clusters_detected
```

---

## 💻 Python API

```python
import pandas as pd
from cyclospora_pyeuk import (
    generate_haplotype_sheet,
    learn_de_novo_haplotypes,
    PyEukDistanceEngine,
    CyclosporaClusterFinder
)

# 1. Option A: Ingest standard genotype calls
sheet_df = generate_haplotype_sheet("specimens_dir_or_zip.zip")

# 1. Option B: Discover haplotypes reference-free from assembled FASTA contigs
# sheet_df, learned_refs = learn_de_novo_haplotypes("cohort_assemblies.fasta")

# 2. Compute KING-wIBS distance matrix (robust to dropouts)
engine = PyEukDistanceEngine()
clean_df = engine.process_haplotype_sheet(sheet_df)
dist_df = engine.compute_revised_wibs_matrix(clean_df)

# 3. Detect outbreak clusters automatically (unsupervised)
finder = CyclosporaClusterFinder()
clusters_df, k, thresh = finder.find_clusters(dist_df, output_dir="results")

print(f"Detected {k} outbreak clusters across {len(clusters_df)} specimens.")
```

---

## 📊 Performance & Validation Highlights

| Benchmark Dataset | Specimen Count ($N$) | Selected $k$ (Label-Free) | Adjusted Rand Index (ARI) | Notes |
| :--- | :---: | :---: | :---: | :--- |
| **CDC Outbreak Benchmark** | 153 | **$k = 2$** | **0.9721** | 99.1% sensitivity, 98.1% specificity against gold standard |
| **Expanded Surveillance Cohort** | 203 | **$k = 2$** | **1.0000** | Perfect 1-to-1 recovery of multi-state outbreaks |
| **De Novo Reference-Free Run** | 12 | **$k = 2$** | **1.0000** | 100% concordance with 0 reference database guidance |

* **Speedup**: Distance matrix computation on $N = 1,078$ national surveillance specimens drops from **24.6 minutes to 14.9 seconds** (99.2× faster).
* **Metric Validity**: Gram matrix PSD projection ensures $\lambda_{\min} \ge 0.0$, eliminating distorted hierarchical tree geometries.

---

## 🌐 Genomic Resources & Documentation

* 🧬 **UCSC BRC-Analytics Pathogen Portal**: [brc-analytics.org](https://brc-analytics.org/) — Reference genome tracks, gene models, and visual browser hubs for *Cyclospora cayetanensis* assemblies ([`GCA_002893315.1`](https://genome.ucsc.edu/cgi-bin/hgTracks?db=GCA_002893315.1) and [`GCA_002893485.1`](https://genome.ucsc.edu/cgi-bin/hgTracks?db=GCA_002893485.1)).
* 📄 **Technical Reports**: Detailed mathematical audits and validation documents are available in the [`docs/`](docs/) directory.

---

## 📄 License & Public Domain

This repository is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/). Code contributions are licensed under the Apache Software License v2.
