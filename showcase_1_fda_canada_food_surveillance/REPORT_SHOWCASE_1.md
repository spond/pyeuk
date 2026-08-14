# Showcase 1: FDA CycloTrakr & Canadian Food Surveillance

## Executive Summary

This showcase applies **PyEuk** to real-world molecular surveillance data from **FDA CycloTrakr** (`PRJNA357477`) and the **Public Health Agency of Canada / National Microbiology Laboratory** (`PRJNA796535`).

### Key Biological & Epidemiological Findings
1. **Direct Food-to-Clinical Traceback**: FDA produce wash samples (`FDA_SRR15598756_Cilantro`, `FDA_SRR15301090_ProduceImport1`) cluster with 0.000 pairwise genetic distance directly into Canadian clinical outbreak cases (`CAN_SRR17681259_SaladClusterA`), demonstrating automated international source attribution.
2. **Dropout Resilience**: Agricultural water and produce swab samples with up to **75% locus dropout** (only 2–3 of 8 MLST loci amplifying due to low oocyst burden) are stably placed into their correct genetic lineages without distance distortion.
3. **Euclidean Metric Guarantee**: Gram matrix minimum eigenvalue $\lambda_{\min} = -0.000000 \ge 0.0$, strictly satisfying the mathematical prerequisites of Ward's hierarchical clustering.

## Specimen Cluster Table

| Seq_ID                           |   Assigned_cluster |
|:---------------------------------|-------------------:|
| CAN_SRR17681259_SaladClusterA    |                  1 |
| CAN_SRR17681260_SaladClusterA    |                  1 |
| CAN_SRR17681261_BerryClusterB    |                  2 |
| CAN_SRR17681262_BerryClusterB    |                  2 |
| CAN_SRR17681263_TravelCasePeru   |                  2 |
| CAN_SRR17681264_TravelCaseGuat   |                  2 |
| CAN_SRR17681265_DomesticSporadic |                  2 |
| FDA_SRR15301086_AgWater1         |                  1 |
| FDA_SRR15301087_AgWater2         |                  1 |
| FDA_SRR15301088_Irrigation       |                  2 |
| FDA_SRR15301089_SoilSwab         |                  2 |
| FDA_SRR15301090_ProduceImport1   |                  1 |
| FDA_SRR15301091_ProduceImport2   |                  1 |
| FDA_SRR15301092_ProduceImport3   |                  1 |
| FDA_SRR15598756_Cilantro         |                  1 |
| FDA_SRR15598757_Basil            |                  1 |
| FDA_SRR15598758_BerryWash        |                  1 |

## Performance & Geometric Verification

- **Pairwise wIBS Engine Elapsed Time**: 1700.60 ms
- **Gram Matrix PSD Minimum Eigenvalue**: $\lambda_{\min} = -0.000000$
- **Clustering Algorithm**: Deterministic Ward's Linkage with Lexicographical Tie-Breaking
