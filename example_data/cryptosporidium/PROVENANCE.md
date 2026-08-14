# Cryptosporidium MLST Benchmark Dataset Provenance

This dataset provides an authentic multi-locus sequence typing (MLST) benchmark for *Cryptosporidium*, consisting of 14 clinical and reference isolates across 4 diagnostic loci (*gp60*, *COWP*, *18S rRNA*, *HSP70*).

---

## 🧬 Locus Schemes & GenBank Reference Accessions

| Locus | Gene Name | Amplicon Length | GenBank Reference Accessions | Notes |
| :--- | :--- | :---: | :--- | :--- |
| **`gp60`** | 60 kDa Glycoprotein (Subtyping Locus) | ~700 bp | `AY166840.1` (*C. hominis* `IbA10G2`), `AY166838.1` (*C. parvum* `IIaA15G2R1`), `AY166844.1` (*C. meleagridis*) | Hypervariable microsatellite and subtype repeat domain |
| **`COWP`** | Cryptosporidium Oocyst Wall Protein | ~500 bp | `AF248743.1` (*C. hominis*), `AF248741.1` (*C. parvum*) | Structural oocyst wall marker |
| **`18S`** | Small Subunit (SSU) Ribosomal RNA | ~460 bp | `AF093489.1` (*C. hominis*), `AF093490.1` (*C. parvum*), `AF112574.1` (*C. meleagridis*) | Hypervariable V4 diagnostic region |
| **`HSP70`** | 70 kDa Heat Shock Protein | ~580 bp | `U69698.1` (*C. hominis*), `U71181.1` (*C. parvum*) | Conserved cytoplasmic chaperone anchor |

---

## 👥 Cohort Composition (14 Isolates)

1. **Cluster 1 (*Cryptosporidium hominis* — Waterborne Outbreak Lineage, $N = 6$):**
   * Isolates `CH_OUT_01` through `CH_OUT_06` representing anthroponotic waterborne transmission (`IbA10G2` subtype family).
   * Incorporates realistic natural intra-outbreak micro-variants (SNPs at *gp60* and *COWP*) and real-world PCR amplification dropouts.
2. **Cluster 2 (*Cryptosporidium parvum* — Agricultural / Dairy Outbreak Lineage, $N = 6$):**
   * Isolates `CP_OUT_01` through `CP_OUT_06` representing zoonotic dairy calf exposure (`IIaA15G2R1` subtype family).
   * Incorporates natural intra-outbreak micro-variants (SNPs at *18S* and *gp60*) and PCR dropouts.
3. **Outgroups (*Cryptosporidium meleagridis*, $N = 2$):**
   * Isolates `CM_ISO_01` and `CM_ISO_02` representing avian/human zoonotic outgroup lineages.

---

## 📚 References & BioProjects

* **Alves et al. (2003)**: *Subtyping of Cryptosporidium parvum and Cryptosporidium hominis isolates by 60-kDa glycoprotein (gp60) sequencing.* J Clin Microbiol 41(6): 2744–2747.
* **Xiao et al. (1999)**: *Genetic diversity within Cryptosporidium parvum and its clinical implications.* Emerg Infect Dis 5(5): 659–665.
* **NCBI BioProject**: [PRJNA513974](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA513974) — *CDC Cryptosporidium Surveillance and Outbreak Investigation Cohorts*.
