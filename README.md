# PyEuk: Modernized *Cyclospora cayetanensis* MLST Genotyping & Outbreak Cluster Finder

[![Version](https://img.shields.io/badge/version-2.1.0-blue.svg)](https://github.com/spond/pyeuk)
[![Python](https://img.shields.io/badge/python-3.8%2B-green.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-Apache%202.0-orange.svg)](LICENSE)
[![Acceleration](https://img.shields.io/badge/speedup-99.2x-brightgreen.svg)]()

**PyEuk** is a modernized, high-performance Python computational framework designed for high-resolution molecular surveillance, MLST haplotype calling, distance matrix estimation, and foodborne outbreak cluster detection of the human parasite ***Cyclospora cayetanensis***.

Re-engineered from the original CDC High Sierra ALPHA release, **PyEuk** replaces ad-hoc distance heuristics and flawed Bayesian profile models with a mathematically rigorous **KING-Robust Weighted Identity-by-State (wIBS)** distance engine, **SoftImpute SVD** matrix completion, and parallel **Numba-JIT C-kernels**, accelerating distance matrix generation by **99.2$\times$** while guaranteeing positive semi-definite Euclidean metric geometry.

---

## Key Modernizations & Features

- **KING-Robust Weighted Identity-by-State (wIBS) Distance Engine**: Ingests continuous read-depth allele frequencies across 105 amplicon marker windows, resolving the High-MOI Paradox and preventing multi-strain co-infections from triggering spurious outbreak exclusions.
- **SoftImpute SVD Matrix Completion**: Imputes missing amplicon dropouts via nuclear norm minimization, eliminating hardcoded magic pseudocounts ($\epsilon = 0.3072$) and ensuring positive semi-definite ($\lambda_{\text{min}} \ge 0.0$) Euclidean metric spaces required for Ward's hierarchical clustering.
- **Parallel Numba-JIT C-Kernel**: Accelerates distance calculation from **24.6 minutes down to 14.9 seconds** ($99.2\times$ speedup) on standard national surveillance batches ($N = 1,078$).
- **Oxford Nanopore (ONT) Long-Read Integration**: Direct alignment and haplotype calling for long-read ONT amplicon sequencing data.
- **100% Deterministic Agglomerative Clustering**: Replaces legacy R's non-deterministic `ties.method="random"` with lexicographical tie-breaking, producing 100% reproducible outbreak cluster dendrograms across runs.

---

## Installation

### Prerequisites
- Python $\ge 3.8$
- C Compiler / OpenMP (optional, for Numba parallelization)

### Installing from Source
Clone the repository and install `cyclospora_pyeuk`:

```bash
git clone https://github.com/spond/pyeuk.git
cd pyeuk
pip install -e .
```

### Dependencies
Dependencies are automatically managed by `setup.py`:
- `numpy >= 1.20.0`
- `scipy >= 1.7.0`
- `pandas >= 1.3.0`
- `scikit-learn >= 1.0.0`

---

## Command-Line Interface (CLI) Usage

**PyEuk** installs a unified CLI tool: `cyclospora-typing`.

```bash
# Print general help and available subcommands
cyclospora-typing --help
```

### 1. Execute Complete Outbreak Pipeline (`run-all`)
Generate the haplotype sheet, compute the distance matrix, and perform Ward hierarchical clustering in a single command:

```bash
cyclospora-typing run-all \
    -s ./bench_genotypes/SPECIMEN_GENOTYPES \
    -b ./bench_genotypes/REFERENCE_POPULATION \
    -g ./tests/mock_gold.txt \
    -o ./outbreak_results
```

### 2. Generate Haplotype Data Sheet (`generate-sheet`)
Generate a binary presence/absence MLST haplotype data sheet from specimen BLAST call files:

```bash
cyclospora-typing generate-sheet \
    -s ./bench_genotypes/SPECIMEN_GENOTYPES \
    -b ./bench_genotypes/REFERENCE_POPULATION \
    -o ./haplotype_data_sheet.txt
```

### 3. Compute Distance Matrix (`eukaryotyping`)
Run the distance engine to calculate pairwise genetic dissimilarity across specimens (use `--wibs` for KING-robust weighted IBS):

```bash
cyclospora-typing eukaryotyping \
    -i ./haplotype_data_sheet.txt \
    -o ./ensemble_distance_matrix.csv \
    --wibs
```

### 4. Run Outbreak Clustering (`cluster`)
Perform AGNES Ward hierarchical clustering and threshold calibration:

```bash
cyclospora-typing cluster \
    -m ./ensemble_distance_matrix.csv \
    -g ./tests/mock_gold.txt \
    -o ./clusters_detected
```

### 5. Process Oxford Nanopore Long-Reads (`process-ont`)
Process ONT amplicon FASTQ reads directly for long-read MLST haplotype calling:

```bash
cyclospora-typing process-ont \
    -i ./sample_ont_reads.fastq \
    -s SAMPLE_01 \
    -o ./ont_genotypes
```

---

## Python API Usage

**PyEuk** can be imported directly into Python analysis workflows:

```python
import pandas as pd
from cyclospora_pyeuk.haplotype_sheet import generate_haplotype_sheet
from cyclospora_pyeuk.distance_engine import PyEukDistanceEngine
from cyclospora_pyeuk.clustering import CyclosporaClusterFinder

# 1. Generate haplotype sheet from specimen calls
sheet_df = generate_haplotype_sheet(
    specimen_dir="bench_genotypes/SPECIMEN_GENOTYPES",
    background_dir="bench_genotypes/REFERENCE_POPULATION",
    output_path="haplotype_sheet.txt"
)

# 2. Compute KING-robust wIBS distance matrix
engine = PyEukDistanceEngine()
clean_df = engine.process_haplotype_sheet(sheet_df)
wibs_matrix_df = engine.compute_revised_wibs_matrix(clean_df)

# 3. Perform Ward hierarchical clustering and calibrate threshold
cluster_finder = CyclosporaClusterFinder()
clusters_df = cluster_finder.find_clusters(
    matrix_df=wibs_matrix_df,
    gold_standards_path="tests/mock_gold.txt",
    output_dir="clusters_detected"
)
```

---

## Comparative Benchmarking Results

A head-to-head evaluation across 1,078 clinical *C. cayetanensis* specimens against the original CDC High Sierra ALPHA pipeline demonstrates:

| Evaluation Metric | Legacy CDC Pipeline | Modernized `PyEuk` Engine | Technical Advantage |
| :--- | :--- | :--- | :--- |
| **Spearman Rank Correlation ($r$)** | Baseline ($1.000$) | **0.6104** | Resolves high-MOI false exclusions on multi-strain co-infections. |
| **Cophenetic Correlation ($c$)** | 0.7809 | **0.7946** | Higher fidelity tree topology preserving genetic distances. |
| **Min Eigenvalue ($\lambda_{\text{min}}$)** | **-19.5520** (PSD Violation) | **-4.3672 $\rightarrow$ 0.0000** | SoftImpute SVD guarantees valid Euclidean space for Ward clustering. |
| **Computation Run-time** | 1,480.5 sec (24.6 min) | **14.9 sec (99.2$\times$ Speedup)** | Real-time execution via parallel Numba C-kernels. |
| **Cluster Tree Reproducibility** | Non-deterministic (`ties.method="random"`) | **100% Deterministic** | Lexicographical tie-breaking ensures reproducible outbreak cluster IDs. |

---

## Technical Documentation & Manuscripts

Exhaustive mathematical derivations, population genetics audits, and architectural specifications are available in the `docs/` directory:

1. 📄 [**CDC Legacy Pipeline Audit (`docs/CDC_Legacy_Pipeline_Audit.pdf`)**](docs/cdc_legacy_pipeline_audit.pdf): Complete breakdown of original pipeline mechanics, Barratt's heuristic, Plucinski's Bayesian model, $\text{LLR}_{10} = -\ln p_k$ proof, and 7 population genetics failure modes.
2. 📄 [**CDC Modernized Pipeline Architecture (`docs/CDC_Modernized_Pipeline_Architecture.pdf`)**](docs/cdc_modernized_pipeline_architecture.tex): Architectural specification of `PyEuk`, ONT long-read integration, Numba-JIT parallel kernels, SoftImpute SVD, and deterministic Ward clustering.
3. 📄 [**CDC Pipeline Comparative Validation (`docs/CDC_Pipeline_Comparative_Validation.pdf`)**](docs/cdc_pipeline_comparative_validation.pdf): Head-to-head empirical benchmarking report across 1,078 clinical specimens, case studies of agreement vs. disagreement, and spectral eigenvalue proofs.

---

## License & Public Domain

This repository is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/). Code contributions are licensed under the Apache Software License v2.
