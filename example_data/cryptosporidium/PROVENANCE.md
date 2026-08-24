# Provenance: Cryptosporidium Multi-Locus Outbreak Benchmark Panel (Synthetic)

This synthetic benchmark dataset comprises **24 designed test specimens across 4 clusters** constructed in code from authentic GenBank reference alleles (NCBI Nuccore accessions below), with engineered single-nucleotide polymorphisms and simulated locus dropout to validate multi-locus frequency-weighted clustering against mosaic allele-sharing structures.

### Reference Loci & GenBank Accessions

| Locus | Lineage / Subtype | GenBank Accession | Amplicon Length | Single-Locus ARI |
| :--- | :--- | :--- | :---: | :---: |
| **18S rRNA** | *C. hominis / parvum* | [`AF093489.1`](https://www.ncbi.nlm.nih.gov/nuccore/AF093489.1) / [`AF093490.1`](https://www.ncbi.nlm.nih.gov/nuccore/AF093490.1) | 458 bp | **0.4651** |
| **HSP70** | *C. hominis / parvum* | [`U69698.1`](https://www.ncbi.nlm.nih.gov/nuccore/U69698.1) / [`U71181.1`](https://www.ncbi.nlm.nih.gov/nuccore/U71181.1) | 550 bp | **0.2133** |
| **COWP** | *C. hominis / parvum* | [`AF248743.1`](https://www.ncbi.nlm.nih.gov/nuccore/AF248743.1) / [`AF248741.1`](https://www.ncbi.nlm.nih.gov/nuccore/AF248741.1) | 483 bp | **0.3349** |
| **gp60** | Subtype IbA10G2 / IaA12G1 | [`AY166840.1`](https://www.ncbi.nlm.nih.gov/nuccore/AY166840.1) / [`AF029759.1`](https://www.ncbi.nlm.nih.gov/nuccore/AF029759.1) | 550 bp | **0.3571** |

### Benchmark Design: Shared Alleles across Clusters

To evaluate frequency-weighted distance estimation and verify that the distance engine cannot take single-locus shortcuts, the panel is constructed so every locus is shared across multiple clusters (mosaic sharing):
- **18S_A** is shared by Cluster 1 and Cluster 2; **18S_B** is shared by Cluster 3 and Cluster 4.
- **HSP70_A** is shared by Clusters 1, 2, and 3; **HSP70_B** is specific to Cluster 4.
- **COWP_A** is shared by Cluster 1 and Cluster 3; **COWP_B** is shared by Cluster 2 and Cluster 4.
- **gp60_A** is shared by Cluster 1 and Cluster 3; **gp60_B** is shared by Cluster 2 and Cluster 4.

Because no single locus separates the cohorts (all single-locus ARIs < 0.50), resolving the true outbreaks requires combining multi-locus information. In this synthetic test fixture, naive unweighted Hamming/Jaccard drops to **ARI = 0.6503**, while PyEuk KING-wIBS achieves **ARI = 0.8836–1.0000**.
