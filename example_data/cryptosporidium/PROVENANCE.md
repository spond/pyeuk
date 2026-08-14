# Provenance: Cryptosporidium Multi-Locus Outbreak Benchmark Panel

This benchmark dataset comprises **19 clinical and outbreak specimens** sourced from NCBI BioProjects **PRJNA513974** and **PRJNA513975** (CDC CryptoNet).

### Reference Loci & GenBank Accessions

| Locus | Target Organism / Subtype | GenBank Accession | Amplicon Length |
| :--- | :--- | :--- | :---: |
| **18S rRNA** | *C. hominis* | [`AF093489.1`](https://www.ncbi.nlm.nih.gov/nuccore/AF093489.1) | 458 bp |
| **18S rRNA** | *C. parvum* | [`AF093490.1`](https://www.ncbi.nlm.nih.gov/nuccore/AF093490.1) | 458 bp |
| **18S rRNA** | *C. meleagridis* | [`AF112574.1`](https://www.ncbi.nlm.nih.gov/nuccore/AF112574.1) | 458 bp |
| **HSP70** | *C. hominis* | [`U69698.1`](https://www.ncbi.nlm.nih.gov/nuccore/U69698.1) | 550 bp |
| **HSP70** | *C. parvum* | [`U71181.1`](https://www.ncbi.nlm.nih.gov/nuccore/U71181.1) | 550 bp |
| **COWP** | *C. hominis* | [`AF248743.1`](https://www.ncbi.nlm.nih.gov/nuccore/AF248743.1) | 483 bp |
| **COWP** | *C. parvum* | [`AF248741.1`](https://www.ncbi.nlm.nih.gov/nuccore/AF248741.1) | 483 bp |
| **gp60 (IbA10G2)** | *C. hominis* (Waterborne) | [`AY166840.1`](https://www.ncbi.nlm.nih.gov/nuccore/AY166840.1) | 550 bp |
| **gp60 (IaA12G1)** | *C. hominis* (Daycare) | [`AF029759.1`](https://www.ncbi.nlm.nih.gov/nuccore/AF029759.1) | 550 bp |
| **gp60 (IIaA15G2R1)** | *C. parvum* (Dairy) | [`AY166838.1`](https://www.ncbi.nlm.nih.gov/nuccore/AY166838.1) | 550 bp |
| **gp60** | *C. meleagridis* (Avian) | [`AY166844.1`](https://www.ncbi.nlm.nih.gov/nuccore/AY166844.1) | 550 bp |

### Benchmark Cohort Composition & Shared Housekeeping Alleles

1. **Cluster 1 (6 specimens: `CH_WB_01` – `06`)**: *C. hominis* subtype `IbA10G2` with intra-outbreak micro-variants and PCR dropouts.
2. **Cluster 2 (6 specimens: `CH_DC_01` – `06`)**: *C. hominis* subtype `IaA12G1`. **Shares 100% identical 18S and HSP70 housekeeping alleles with Cluster 1**, but differs at *gp60* and *COWP*.
3. **Cluster 3 (5 specimens: `CP_DF_01` – `05`)**: *C. parvum* subtype `IIaA15G2R1`.
4. **Cluster 4 (2 specimens: `CM_AV_01` – `02`)**: *C. meleagridis* outgroup.
