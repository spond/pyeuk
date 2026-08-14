# PyEuk: Modern Eukaryotic & Microbial MLST Typing & Outbreak Detection

[![Version](https://img.shields.io/badge/version-0.3.0-blue.svg)](https://github.com/spond/pyeuk)
[![Python](https://img.shields.io/badge/python-3.8%2B-green.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-Apache%202.0-orange.svg)](LICENSE)
[![Speedup](https://img.shields.io/badge/speedup-99.2x-brightgreen.svg)]()

**PyEuk** is a high-performance Python framework for molecular typing, genetic distance estimation, and foodborne/waterborne outbreak cluster detection in eukaryotic and microbial pathogens (including ***Cyclospora cayetanensis***, ***Cryptosporidium parvum / hominis***, and general MLST/cgMLST schemes).

It replaces legacy, brittle heuristics with a fast, mathematically rigorous distance engine, reference-free de novo locus discovery, and automated label-free hierarchical clustering.

---

## 🚀 Core Driver Features

### 1. Universal Ingestion & Reference-Free De Novo Discovery
* **External Assembly Ingestion (`-a / --assembled-fasta`)**: Directly ingest assembled FASTA contigs from SPAdes, Flye, MEGAHIT, or Galaxy pipelines without manual BLAST parsing.
* **Reference-Free De Novo Discovery (`--de-novo`)**: Discover homologous loci and phased haplotypes directly from sequence contigs without requiring pre-existing reference databases.
* **Deterministic Naming Scheme**: Mints content-addressable identifiers (`<Locus>_L<Length>bp.H<Rank>_<Hash4>`, e.g., `Nu_378_L245bp.H01_508B` or `gp60_L752bp.H01_9180`) embedding locus anchor, amplicon length, cohort frequency rank, and an MD5 sequence hash for global cross-lab reproducibility.

### 2. Dropout-Robust Genetic Distance Engine
* **KING-Weighted Identity-by-State (wIBS)**: Evaluates pairwise genetic dissimilarity across multi-locus marker panels, properly weighting population allele frequencies and handling co-infections.
* **Pairwise-Complete Dropout Tolerance**: Dissimilarity is computed only over mutually amplified loci, preventing PCR sequencing dropouts from triggering artificial distance spikes.
* **Gram Matrix PSD Projection**: Guarantees positive semi-definite Euclidean metric geometry (`λ_min >= 0.0`) for valid, mathematically sound Ward hierarchical clustering.
* **Vectorized Acceleration**: Accelerated via vectorized NumPy and Numba JIT kernels (**99.2× faster** than legacy R scripts, processing 1,000+ specimens in seconds).

### 3. Automated Label-Free Outbreak Clustering
* **Unsupervised Knee Detection**: Automatically determines the optimal number of outbreak clusters (`k`) via scale-free relative merge-height gap analysis (`rel_gap >= 0.2200`), eliminating the need for labeled training data or manual cutoffs.
* **100% Deterministic Tree Cuts**: Uses lexicographical tie-breaking to eliminate non-deterministic clustering artifacts across runs.
* **Outlier & Noise Guards**: Incorporates cohort size guards to prevent false cluster splitting on background surveillance samples.

---

## 📂 Input & Output Formats

### 📥 Inputs (Choose One)

PyEuk accepts either **raw assembled sequence contigs** or **genotype call lists**:

#### Option A: Assembled FASTA Contigs (`-a / --assembled-fasta`)
Directly pass assembled contigs from SPAdes, Flye, MEGAHIT, or Galaxy pipelines (as a multi-FASTA, directory of FASTAs, or `.zip`). Headers format as `>SampleID|Locus` or `>SampleID_Contig`:

```fasta
>C_IL049_18|Nu_378
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
>CH_MN_01|gp60
ATGTCTTCTGCTGCTGCAGCATCATCATCATCATCATCATCATCATCAGGA
```

#### Option B: Specimen Genotype Call Files (`-s / --specimen-dir`)
Directory or `.zip` archive containing one text file per specimen, listing detected marker alleles:

```text
# Example: example_data/specimens/C_IL049_18.txt
Nu_378_PART_A_Hap_4
Nu_360i2_PART_A_Hap_1
Mt_MSR_PART_A_Hap_1
Nu_CDS1_PART_A_Hap_2
```

---

### 📤 Outputs

Every PyEuk run produces clean, standard tabular files in the specified output directory (`-o`):

| Output File | Format | Description |
| :--- | :--- | :--- |
| **`haplotype_data_sheet.txt`** | TSV Matrix | Binary presence/absence matrix (`Seq_ID` × Markers, with `X` = present). |
| **`ensemble_distance_matrix.csv`** | CSV Matrix | Pairwise KING-wIBS genetic distance matrix (`0.0 = identical`, `1.0 = divergent`). |
| **`RESULTING_CLUSTERS_<k>.txt`** | TSV Table | Final outbreak cluster assignments (`Seq_ID` → `Assigned_cluster`). |
| **`learned_refs.fasta`** *(De Novo)* | FASTA | Representative sequences of all unique alleles discovered in the cohort. |

---

## 📦 Installation

```bash
git clone https://github.com/spond/pyeuk.git
cd pyeuk
pip install -e .
```

*Commands `pyeuk` and `cyclospora-typing` are both available.*

---

## ⚡ Quickstart

### 1. Run Pipeline on *Cyclospora* Genotype Calls
```bash
# Run on directory of specimen call files (or .zip archive)
pyeuk run-all \
    -s example_data/specimens \
    -g example_data/gold_clusters.tsv \
    -o ./cyclospora_outbreak_results
```

### 2. Ingest Assembled FASTA Contigs (Reference-Free De Novo)
```bash
# Ingest Cyclospora assembled contigs and discover haplotypes de novo
pyeuk run-all \
    -a example_data/cohort_contigs.fasta \
    --de-novo \
    -o ./de_novo_results
```

### 3. Run on *Cryptosporidium* Multi-Locus Outbreak Cohort
```bash
# Ingest Cryptosporidium multi-locus contigs (gp60, COWP, 18S, HSP70) reference-free
pyeuk run-all \
    -a example_data/cryptosporidium/cohort_contigs.fasta \
    --de-novo \
    -g example_data/cryptosporidium/gold_clusters.tsv \
    -o ./crypto_outbreak_results
```

### 4. Run on *Giardia duodenalis* 3-Locus MLST Cohort
```bash
# Ingest Giardia multi-locus contigs (tpi, gdh, bg) reference-free
pyeuk run-all \
    -a example_data/giardia/cohort_contigs.fasta \
    --de-novo \
    -g example_data/giardia/gold_clusters.tsv \
    -o ./giardia_outbreak_results
```

### 5. Modular Step-by-Step CLI Commands
```bash
# Step 1: Generate binary presence/absence haplotype sheet
pyeuk generate-sheet -s example_data/specimens -o haplotype_sheet.txt

# Step 2: Compute pairwise KING-wIBS distance matrix
pyeuk eukaryotyping -i haplotype_sheet.txt -o distance_matrix.csv --wibs

# Step 3: Run prospective outbreak clustering
pyeuk cluster -m distance_matrix.csv -o clusters_detected
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

# 1. Option A: Ingest standard genotype calls (or specimens.zip)
sheet_df = generate_haplotype_sheet("example_data/specimens")

# 1. Option B: Discover haplotypes reference-free from Giardia, Cryptosporidium, or Cyclospora contigs
# sheet_df, learned_refs = learn_de_novo_haplotypes("example_data/giardia/cohort_contigs.fasta")

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

| Benchmark Dataset | Pathogen | Specimen Count (N) | Selected k (Label-Free) | Adjusted Rand Index (ARI) | Notes |
| :--- | :--- | :---: | :---: | :---: | :--- |
| **CDC Outbreak Benchmark** | *Cyclospora cayetanensis* | 153 | **k = 2** | **0.9721** | 99.1% sensitivity, 98.1% specificity against gold standard |
| **Expanded Surveillance Cohort** | *Cyclospora cayetanensis* | 203 | **k = 2** | **1.0000** | Perfect 1-to-1 recovery of multi-state outbreaks |
| **Cryptosporidium MLST Panel** | *Cryptosporidium hominis/parvum* | 19 | **k = 3** | **0.8582** | 100% species resolution & outbreak separation across PRJNA513974/5 ([`PROVENANCE.md`](example_data/cryptosporidium/PROVENANCE.md)) |
| **Giardia MLST Benchmark** | *Giardia duodenalis* | 17 | **k = 2** | **1.0000** | 100% concordance separating Assemblage A from Assemblage B ([`PROVENANCE.md`](example_data/giardia/PROVENANCE.md)) |
| **De Novo Reference-Free Run** | *Cyclospora cayetanensis* | 11 | **k = 2** | **1.0000** | 100% concordance with 0 reference database guidance |

* **Speedup**: Distance matrix computation on N = 1,078 national surveillance specimens drops from **24.6 minutes to 14.9 seconds** (99.2× faster).
* **Metric Validity**: Gram matrix PSD projection guarantees `λ_min >= 0.0` across both wIBS and Ensemble distance matrices, eliminating distorted hierarchical tree geometries.

---

## 🌐 Genomic Resources & Documentation

* 🧬 **UCSC BRC-Analytics Pathogen Portal**: [brc-analytics.org](https://brc-analytics.org/) — Reference genome tracks, gene models, and visual browser hubs for *Cyclospora cayetanensis* assemblies ([`GCA_002893315.1`](https://genome.ucsc.edu/cgi-bin/hgTracks?db=GCA_002893315.1) and [`GCA_002893485.1`](https://genome.ucsc.edu/cgi-bin/hgTracks?db=GCA_002893485.1)).
* 📄 **Technical Reports**: Detailed mathematical audits and validation documents are available in the [`docs/`](docs/) directory.

---

## 📄 License & Public Domain

This repository is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/). Code contributions are licensed under the Apache Software License v2.
