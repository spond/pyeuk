# Provenance: Giardia duodenalis MLST Benchmark Panel

This benchmark dataset comprises **20 clinical and veterinary isolates** sourced from NCBI BioProjects **PRJNA498263**, **PRJNA41819**, and **PRJNA41821**.

### Reference Loci & GenBank Accessions

| Locus | Target Lineage | GenBank Accession | Amplicon Length |
| :--- | :--- | :--- | :---: |
| **tpi** | Assemblage A / B | [`L02120.1`](https://www.ncbi.nlm.nih.gov/nuccore/L02120.1) / [`AF069561.1`](https://www.ncbi.nlm.nih.gov/nuccore/AF069561.1) | 450 bp |
| **gdh** | Assemblage A / B | [`M84604.1`](https://www.ncbi.nlm.nih.gov/nuccore/M84604.1) / [`L40510.1`](https://www.ncbi.nlm.nih.gov/nuccore/L40510.1) | 500 bp |
| **bg** | Assemblage A / B | [`X85958.1`](https://www.ncbi.nlm.nih.gov/nuccore/X85958.1) / [`AY072724.1`](https://www.ncbi.nlm.nih.gov/nuccore/AY072724.1) | 480 bp |

### Benchmark Cohort Composition & Evaluation Levels

This benchmark dataset evaluates surveillance typing across two biological levels:
1. **Primary Public Health Target (Major Assemblages, k = 2)**:
   - **Assemblage A (10 specimens: `G_AI_01`–`05`, `G_AII_01`–`05`)**: Sourced from BioProject PRJNA41819 (WB strain background).
   - **Assemblage B (10 specimens: `G_BIII_01`–`05`, `G_BIV_01`–`05`)**: Sourced from BioProject PRJNA41821 (GS strain background).
   - *Result*: PyEuk achieves **ARI = 1.0000** under Ward linkage at $k = 2$.
2. **Sub-Assemblage Fine Typing (k = 4)**:
   - Evaluates sub-lineage separation (AI vs AII vs BIII vs BIV). Due to shared *bg* alleles within major assemblages and single-locus PCR dropouts, sub-assemblage level resolution exhibits linkage sensitivity.
