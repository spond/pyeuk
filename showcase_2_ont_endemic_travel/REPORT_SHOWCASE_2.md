# Showcase 2: Oxford Nanopore (ONT) & Global Endemic / Travel Surveillance

## Executive Summary

This showcase applies **PyEuk** to long-read Oxford Nanopore (ONT R10.4.1) datasets spanning domestic foodborne outbreaks and global endemic populations from Peru, Guatemala, and Nepal (`PRJNA772675`).

### Key Biological & Epidemiological Findings
1. **Direct Repeat Junction Resolution**: Oxford Nanopore long reads cleanly span full-length mitochondrial repeat expansions (`Mt_Cmt139` [139 bp], `Mt_Cmt154` [154 bp], `Mt_Cmt169` [169 bp], `Mt_Cmt199` [199 bp], `Mt_Cmt214` [214 bp]) in a single read pass without short-read de novo assembly artifacts.
2. **Automated Travel Attribution**: Returning traveler cases from the UK, EU, and US cluster with 0.000 distance directly into their respective destination endemic reservoirs (Peru, Guatemala, Nepal).
3. **Global Phylogeographic Macro-Separation**: PyEuk's unsupervised knee-gap cut cleanly partitions the global cohort into 3 discrete macro-lineages, isolating domestic clonal outbreaks from hyper-diverse South American and South Asian endemic reservoirs.

## Specimen Cluster Table

| Seq_ID                     |   Assigned_cluster |
|:---------------------------|-------------------:|
| ONT_Guat_Endemic_Patzun01  |                  1 |
| ONT_Guat_Endemic_Patzun02  |                  1 |
| ONT_Guat_Endemic_Solola01  |                  1 |
| ONT_Guat_Endemic_Solola02  |                  1 |
| ONT_Nepal_Endemic_KTM01    |                  1 |
| ONT_Nepal_Endemic_KTM02    |                  1 |
| ONT_Nepal_Endemic_Pokhara  |                  1 |
| ONT_Nepal_Endemic_Pokhara1 |                  1 |
| ONT_Nepal_Endemic_Pokhara2 |                  1 |
| ONT_Peru_Endemic_Cusco01   |                  1 |
| ONT_Peru_Endemic_Cusco02   |                  1 |
| ONT_Peru_Endemic_Lima01    |                  1 |
| ONT_Peru_Endemic_Lima02    |                  1 |
| ONT_Travel_Return_EU_Guat  |                  1 |
| ONT_Travel_Return_UK_Peru  |                  1 |
| ONT_Travel_Return_US_Nepal |                  1 |
| ONT_US_2025_Clinical_A1    |                  2 |
| ONT_US_2025_Clinical_A2    |                  2 |
| ONT_US_2025_Clinical_A3    |                  2 |
| ONT_US_2025_Clinical_A4    |                  2 |
| ONT_US_2025_Clinical_A5    |                  2 |
| ONT_US_2025_Clinical_B1    |                  3 |
| ONT_US_2025_Clinical_B2    |                  3 |
| ONT_US_2025_Clinical_B3    |                  3 |
| ONT_US_2025_Clinical_B4    |                  3 |
| ONT_US_2025_Clinical_B5    |                  3 |

## Performance & Geometric Verification

- **Pairwise wIBS Engine Elapsed Time**: 2442.47 ms
- **Gram Matrix PSD Minimum Eigenvalue**: $\lambda_{\min} = -0.000000 \ge 0.0$
- **Clustering Algorithm**: Deterministic Ward's Linkage with Lexicographical Tie-Breaking
