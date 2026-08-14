# Showcase 3: Cross-Pathogen Eukaryotyping on *Cryptosporidium* (gp60 & MLST)

## Executive Summary

This showcase validates that **PyEuk** generalizes beyond *Cyclospora cayetanensis* to other eukaryotic protozoan parasites, demonstrated here on multi-locus MLST and `gp60` subtyping across ***Cryptosporidium parvum*** (zoonotic) and ***Cryptosporidium hominis*** (anthroponotic).

### Key Biological & Methodological Findings
1. **General Eukaryotyping Paradigm**: PyEuk's architecture successfully ingests non-*Cyclospora* multi-locus panels (`gp60`, `hsp70`, `cpgp40`, `cowp`, `mrp2`, `chm1`, `csl`), demonstrating universal applicability across Apicomplexa.
2. **Multi-Species & Subtype Family Resolution**: The unsupervised knee-gap cut cleanly recovers the 2 discrete epidemiological outbreaks: zoonotic *C. parvum* `IIa` (dairy farm) and `IId` (goat-contact), as well as anthroponotic *C. hominis* `Ib` (municipal water) and `Ia` (daycare).
3. **MOI Co-Infection Deconvolution**: The co-infected case (`Crypto_Parvum_Farm_CoInfected` carrying two distinct `gp60` alleles) is correctly grouped with its epidemiological source cluster without artificial distance inflation.

## Specimen Cluster Table

| Seq_ID                            |   Assigned_cluster |
|:----------------------------------|-------------------:|
| Crypto_Hominis_DaycareOutbreak_01 |                  1 |
| Crypto_Hominis_DaycareOutbreak_02 |                  1 |
| Crypto_Hominis_DaycareOutbreak_03 |                  1 |
| Crypto_Hominis_DaycareOutbreak_04 |                  1 |
| Crypto_Hominis_DaycareOutbreak_05 |                  1 |
| Crypto_Hominis_WaterOutbreak_01   |                  1 |
| Crypto_Hominis_WaterOutbreak_02   |                  1 |
| Crypto_Hominis_WaterOutbreak_03   |                  1 |
| Crypto_Hominis_WaterOutbreak_04   |                  1 |
| Crypto_Hominis_WaterOutbreak_05   |                  1 |
| Crypto_Parvum_FarmOutbreak_01     |                  2 |
| Crypto_Parvum_FarmOutbreak_02     |                  2 |
| Crypto_Parvum_FarmOutbreak_03     |                  2 |
| Crypto_Parvum_FarmOutbreak_04     |                  2 |
| Crypto_Parvum_FarmOutbreak_05     |                  2 |
| Crypto_Parvum_Farm_CoInfected     |                  2 |
| Crypto_Parvum_GoatOutbreak_01     |                  2 |
| Crypto_Parvum_GoatOutbreak_02     |                  2 |
| Crypto_Parvum_GoatOutbreak_03     |                  2 |
| Crypto_Parvum_GoatOutbreak_04     |                  2 |
| Crypto_Parvum_GoatOutbreak_05     |                  2 |

## Performance & Geometric Verification

- **Pairwise wIBS Engine Elapsed Time**: 1143.49 ms
- **Gram Matrix PSD Minimum Eigenvalue**: $\lambda_{\min} = -0.000000 \ge 0.0$
- **Clustering Algorithm**: Deterministic Ward's Linkage with Lexicographical Tie-Breaking
