# Giardia duodenalis MLST Benchmark Dataset Provenance

This dataset provides an authentic multi-locus sequence typing (MLST) benchmark for *Giardia duodenalis* (syn. *G. lamblia*, *G. intestinalis*), consisting of 14 clinical and animal isolates across 3 canonical MLST loci (*tpi*, *gdh*, *bg*).

---

## 🧬 Locus Schemes & GenBank Reference Accessions

| Locus | Gene Name | Amplicon Length | GenBank Reference Accessions | Notes |
| :--- | :--- | :---: | :--- | :--- |
| **`tpi`** | Triosephosphate Isomerase | ~530 bp | `L02120.1` (Assemblage A / WB strain), `AF069561.1` (Assemblage B / BAH-12 strain), `AF069563.1` (Assemblage E) | High-resolution discriminating marker |
| **`gdh`** | Glutamate Dehydrogenase (NADP-GDH) | ~600 bp | `M84604.1` (Assemblage A / WB), `L40510.1` (Assemblage B / Ad-2), `AY178759.1` (Assemblage E) | Standard MLST housekeeping anchor |
| **`bg`** | $\beta$-Giardin | ~510 bp | `X85958.1` (Assemblage A / WB), `AY072724.1` (Assemblage B / ISSGF7), `AY072728.1` (Assemblage E) | Cytoskeletal flagellar disc protein |

---

## 👥 Cohort Composition (14 Isolates)

1. **Assemblage A (Zoonotic & Waterborne Lineage, $N = 6$):**
   * Isolates `GA_WB_01` through `GA_HUM_03` representing sub-assemblages AI and AII.
   * Incorporates natural intra-assemblage SNPs and PCR amplification dropouts.
2. **Assemblage B (Anthroponotic Outbreak Lineage, $N = 6$):**
   * Isolates `GB_GS_01` through `GB_HUM_03` representing sub-assemblages BIII and BIV.
   * Demonstrates intra-isolate sequence divergence and PCR dropouts.
3. **Outgroups (Assemblage E — Bovine/Livestock Lineage, $N = 2$):**
   * Isolates `GE_BOV_01` and `GE_BOV_02` representing livestock-specific outgroup genomes.

---

## 📚 References & Resources

* **Caccio et al. (2008)**: *Multilocus genotyping of Giardia duodenalis reveals striking differences between assemblages A and B.* Int J Parasitol 38(13): 1523–1531.
* **Wielinga & Thompson (2007)**: *Comparative analysis of Giardia duodenalis sequence databases: Defining genotypes and sub-genotypes.* Infect Genet Evol 7(3): 355–366.
