# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.6.0] - 2026-08-27

### Changed
* **Top-Level Package Renamed to `pyeuk`**:
  - Renamed top-level Python import package from `cyclospora_pyeuk` to `pyeuk` to reflect panel- and species-agnostic eukaryotic MLST genotyping capabilities across *Cyclospora*, *Cryptosporidium*, *Giardia*, *Plasmodium vivax*, and nematodes (resolves #18).
  - Maintained complete backward compatibility via a `cyclospora_pyeuk` deprecation shim package emitting `DeprecationWarning` and re-exporting all API symbols and submodules.
  - Updated `pyproject.toml` distribution configuration (`packages = ["pyeuk", "pyeuk.amplicon", "cyclospora_pyeuk", "cyclospora_pyeuk.amplicon"]`).

---

## [0.5.0] - 2026-08-27

### Added
* **Amplicon Front End (`cyclospora_pyeuk.amplicon`)**:
  - `pyeuk define-windows`: Chooses cohort analysis windows and spannable cores from read data without requiring curated reference BEDs (resolves #14, #15).
  - `pyeuk call-haplotypes`: Directly calls window haplotypes on single spanning reads with left-aligned indel normalization and HGVS-like content-derived nomenclature (resolves #14, #15).
  - `pyeuk build-sheet`: Assembles specimen-by-haplotype presence/absence sheets, haplotype mapping catalogs, and long-format observation tables with per-window frequency and minor-allele filtering (resolves #14, #15).
  - Configured optional dependency `cyclospora_pyeuk[amplicon]` (`pysam>=0.22`) with lazy import guards ensuring core functionality installs and runs cleanly without pysam.
* **Unified Naming Contract (`cyclospora_pyeuk.naming`)**:
  - Centralized bidirectional parsing and formatting via `name_haplotype()` and `parse_locus_name()`, enforcing invariant round-tripping across CDC legacy, de novo, windowed amplicon, and compact formatting styles (resolves #14).
* **Distance-Threshold Linkage Clustering**:
  - Added `--cut distance` mode to `CyclosporaClusterFinder` and CLI `cluster` / `run-all` subcommands to cut hierarchical dendrograms at fixed dissimilarity thresholds for surveillance cohorts with abundant singletons (resolves #15).
  - Added `suggest_linkage_threshold()` supporting robust MAD-based calibration from labelled within-cluster pairs or 5th percentile distance distribution heuristics with transparent provenance reporting (resolves #15).
  - Added `--linkage-method` (`ward`, `single`, `average`, `complete`) to follow transmission chains natively (resolves #15).
  - Added informative warning when `k_max < n/2` in unsupervised count mode.

---

## [0.4.0] - 2026-08-25

### Changed
* **Default Allele Frequency Weighting Scheme**:
  - Changed default `weight_mode` in `PyEukDistanceEngine` to `"heterozygosity"` (`w = 2 * p * (1 - p)`), scaling with binomial variance for uncentered binary presence/absence indicator matrices. This concentrates weight on balanced, outbreak-discriminating loci while attenuating rare singletons (resolves #11).
  - Supported weight modes: `heterozygosity` (default), `inverted-king` (`w = sqrt(p * (1 - p))`), `king` (`w = 1 / sqrt(p * (1 - p))` for centered dosage comparison), and `uniform` (`w = 1.0`).
* **Honored `k_min` in Label-Free Clustering**:
  - In `CyclosporaClusterFinder.find_clusters`, candidate search now begins at `search_start = max(2, k_min)`. When `k_min > 2` is specified and no candidate clears the relative gap floor, `find_clusters` selects the highest-scoring candidate satisfying `k >= k_min` and cluster size guards rather than collapsing to `k = 1` (resolves #12).
* **Locus-Level Completeness and Distance Matching**:
  - Fixed denominator calculation to match shared locus presence at the locus level, eliminating dead-weight denominator dilution from unobserved allele column placeholders (resolves #8).
  - Locus completeness calculation made strictly invariant to unused/unobserved allele column presence (resolves #7).

### Added
* **Minor Allele Frequency (MAF) Filter**:
  - Added `--min-maf` parameter (default `0.0`) to `PyEukDistanceEngine`, `pyeuk eukaryotyping`, and `pyeuk run-all` to filter near-invariant singleton columns (`p < min_maf` or `p > 1 - min_maf`) (resolves #11).
* **Raw Pairwise Distance Export**:
  - Added `--no-psd` flag (`project_psd=False`) across `compute_revised_wibs_matrix`, `compute_snp_weighted_wibs_matrix`, and `compute_ensemble_matrix` to allow exporting exact raw distance matrices without Gram PSD projection loss (resolves #11).
* **Configurable Relative Gap Floor**:
  - Exposed `relative_gap_floor: float = 0.2200` parameter on `CyclosporaClusterFinder` and via `--relative-gap-floor` on CLI `cluster` and `run-all` commands (resolves #12).
  - Added `finder.last_selection_meta` dictionary recording selection status (`optimal`, `floor_override`, `unsatisfiable_constraint`, `single_group`, or `trivial`), gap metrics, and failure diagnostics (resolves #12).
* **Transparent Quality Reporting**:
  - Specimens failing completeness criteria are explicitly included in cluster outputs as Cluster -1 (Unassigned) (resolves #7).
* **Dataset Documentation & Provenance**:
  - Clarified in `example_data/cryptosporidium/PROVENANCE.md` and `example_data/giardia/PROVENANCE.md` that *Cryptosporidium* and *Giardia* evaluation panels are synthetic mosaic/assemblage benchmark fixtures constructed from GenBank references (resolves #10).

---

## [0.3.0] - 2026-08-14

### Added
* **Multi-Pathogen Support**:
  - Extended genotyping pipeline from *Cyclospora cayetanensis* to *Cryptosporidium* and *Giardia duodenalis*.
  - Added `--ploidy` argument across distance engines to support haploid, diploid, and polyploid eukaryotic pathogens.
* **De Novo Reference-Free Discovery**:
  - Added `learn_de_novo_haplotypes` and `--de-novo` flag to ingest assembled FASTA contigs and discover loci and alleles reference-free.
* **Nanopore & PacBio Long-Read Ingestion**:
  - Added `process-ont` CLI subcommand and `NanoporeAmpliconProcessor` with quality filtering and dynamic alignment against reference databases.
* **Gram Matrix PSD Projection**:
  - Integrated eigenvalue clipping on double-centered Gram matrices to guarantee positive semi-definiteness.
* **Packaging & Distribution**:
  - Added `pyproject.toml` and standard CLI entrypoints (`pyeuk`, `cyclospora-typing`).
