"""
Unit tests for cyclospora_pyeuk python package.
"""

import os
import unittest
import numpy as np
import pandas as pd
from cyclospora_pyeuk.haplotype_sheet import generate_haplotype_sheet
from cyclospora_pyeuk.distance_engine import PyEukDistanceEngine
from cyclospora_pyeuk.clustering import CyclosporaClusterFinder


class TestCyclosporaPyEuk(unittest.TestCase):

    def setUp(self):
        # Create a mock genotype dataframe for testing
        self.mock_data = pd.DataFrame({
            "Seq_ID": ["Sample_A", "Sample_B", "Sample_C", "Sample_D", "Sample_E"],
            "Mt_Cmt_PART_A_Hap_1": ["X", "X", "", "X", "X"],
            "Mt_Cmt_PART_A_Hap_2": ["", "", "X", "", ""],
            "Mt_MSR_PART_A_Hap_1": ["X", "X", "X", "X", ""],
            "Nu_360i2_PART_A_Hap_1": ["X", "X", "X", "", "X"],
            "Nu_378_PART_A_Hap_1": ["X", "X", "X", "X", "X"],
            "Nu_CDS1_PART_A_Hap_1": ["X", "", "X", "X", "X"],
        })

    def test_haplotype_processing(self):
        engine = PyEukDistanceEngine()
        clean_df = engine.process_haplotype_sheet(self.mock_data)
        self.assertIn("Seq_ID", clean_df.columns)
        self.assertGreater(len(clean_df), 0)

    def test_distance_engine_execution(self):
        engine = PyEukDistanceEngine(epsilon=0.3072)
        res_df = engine.compute_ensemble_matrix(self.mock_data)
        self.assertIsInstance(res_df, pd.DataFrame)
        self.assertEqual(res_df.shape[0], res_df.shape[1])
        np.testing.assert_allclose(np.diag(res_df.values), 0.0, atol=1e-6)

    def test_revised_wibs_and_psd_guarantee(self):
        engine = PyEukDistanceEngine()
        wibs_df = engine.compute_revised_wibs_matrix(self.mock_data)
        self.assertIsInstance(wibs_df, pd.DataFrame)
        self.assertEqual(wibs_df.shape[0], wibs_df.shape[1])
        np.testing.assert_allclose(np.diag(wibs_df.values), 0.0, atol=1e-6)

        # Assert Gram matrix PSD property (lambda_min >= -1e-12)
        nids = len(wibs_df)
        H = np.eye(nids) - np.ones((nids, nids)) / float(nids)
        G = -0.5 * H @ (wibs_df.values ** 2) @ H
        min_eval = np.min(np.linalg.eigvalsh(G))
        self.assertGreaterEqual(min_eval, -1e-12)

    def test_prospective_unsupervised_clustering(self):
        engine = PyEukDistanceEngine()
        wibs_df = engine.compute_revised_wibs_matrix(self.mock_data)
        finder = CyclosporaClusterFinder()
        clusters_df, k, thresh = finder.find_clusters(wibs_df, gold_file_path=None, k_min=1, k_max=3, output_dir="tests/out_clusters")
        self.assertIsNotNone(clusters_df)
        self.assertIn("Assigned_cluster", clusters_df.columns)
        self.assertGreaterEqual(k, 1)

    def test_micro_cluster_detection(self):
        engine = PyEukDistanceEngine()
        wibs_df = engine.compute_revised_wibs_matrix(self.mock_data)
        finder = CyclosporaClusterFinder()
        micro_clusters = finder.detect_micro_clusters(wibs_df, micro_threshold=0.10, min_micro_size=2)
        self.assertIsInstance(micro_clusters, list)

    def test_snp_weighted_wibs_matrix(self):
        engine = PyEukDistanceEngine()
        snp_df = engine.compute_snp_weighted_wibs_matrix(self.mock_data)
        self.assertIsInstance(snp_df, pd.DataFrame)
        self.assertEqual(snp_df.shape[0], snp_df.shape[1])
        np.testing.assert_allclose(np.diag(snp_df.values), 0.0, atol=1e-6)

    def test_external_assembly_ingestion(self):
        import tempfile
        from cyclospora_pyeuk.haplotype_sheet import generate_haplotype_sheet_from_assemblies

        # Create temporary FASTA with mock assembled contigs
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">Sample_01|Nu_378_PART_A_Hap_1\nATGCGATCGATCGATCGATCGATCGATCGA\n")
            f.write(">Sample_01|Nu_360i2_PART_A_Hap_1\nCGATCGATCGATCGATCGATCGATCGATCG\n")
            f.write(">Sample_02|Nu_378_PART_A_Hap_1\nATGCGATCGATCGATCGATCGATCGATCGA\n")
            temp_path = f.name

        try:
            sheet_df = generate_haplotype_sheet_from_assemblies(
                assembled_input=temp_path,
                output_path="tests/out_clusters/mock_assembled_sheet.txt"
            )
            self.assertIsInstance(sheet_df, pd.DataFrame)
            self.assertIn("Seq_ID", sheet_df.columns)
            self.assertEqual(len(sheet_df), 2)
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)

    def test_de_novo_haplotype_discovery(self):
        import tempfile
        from cyclospora_pyeuk.haplotype_sheet import learn_de_novo_haplotypes

        # Create temporary FASTA without any predefined references
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            # Sample 1 has alleles at Locus_A and Locus_B
            f.write(">Patient_1|Locus_A\nATGCGATCGATCGATCGATCGATCGATCGA\n")
            f.write(">Patient_1|Locus_B\nTTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA\n")
            # Sample 2 shares Locus_A allele but has a novel Locus_A variant (co-infection)
            f.write(">Patient_2|Locus_A\nATGCGATCGATCGATCGATCGATCGATCGA\n")
            f.write(">Patient_2|Locus_A\nATGCGATCGATCGATCGATCGATCGATCGTT\n")
            # Sample 3 only has Locus_B
            f.write(">Patient_3|Locus_B\nTTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA\n")
            temp_path = f.name

        out_tsv = "tests/out_clusters/de_novo_sheet.txt"
        out_fa = "tests/out_clusters/de_novo_refs.fasta"

        try:
            sheet_df, learned_dict = learn_de_novo_haplotypes(
                assembled_input=temp_path,
                output_path=out_tsv,
                output_fasta=out_fa
            )
            self.assertIsInstance(sheet_df, pd.DataFrame)
            self.assertEqual(len(sheet_df), 3)
            self.assertIn("Seq_ID", sheet_df.columns)
            self.assertGreater(len(learned_dict), 0)
            self.assertTrue(os.path.exists(out_tsv))
            self.assertTrue(os.path.exists(out_fa))
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)

    def test_example_data_pipeline(self):
        """Tests that the example_data directory works across standard, zip, and de novo modes."""
        from cyclospora_pyeuk.haplotype_sheet import generate_haplotype_sheet, learn_de_novo_haplotypes
        from cyclospora_pyeuk.distance_engine import PyEukDistanceEngine
        from cyclospora_pyeuk.clustering import CyclosporaClusterFinder

        # 1. Test from directory
        sheet_dir = generate_haplotype_sheet("example_data/specimens")
        self.assertEqual(len(sheet_dir), 11)

        # 2. Test from zip
        sheet_zip = generate_haplotype_sheet("example_data/specimens.zip")
        self.assertEqual(len(sheet_zip), 11)

        # 3. Test de novo from FASTA contigs
        sheet_dn, learned = learn_de_novo_haplotypes("example_data/cohort_contigs.fasta")
        self.assertEqual(len(sheet_dn), 11)
        self.assertGreater(len(learned), 20)

        # 4. Test clustering on de novo sheet
        engine = PyEukDistanceEngine(min_completeness=0.0)
        wibs = engine.compute_revised_wibs_matrix(sheet_dn)
        finder = CyclosporaClusterFinder()
        clusters, k, _ = finder.find_clusters(wibs, gold_file_path="example_data/gold_clusters.tsv")
        self.assertEqual(k, 2)

    def test_hand_calculated_wibs_pin_and_mutant_kills(self):
        """
        Pins exact hand-calculated KING-wIBS distances on a 3-specimen x 2-locus fixture with dropout.
        Explicitly asserts that mutant implementations (inverted weights, uniform weights, disabled mask) fail.
        """
        # LocA (ploidy 1): S1=A1, S2=A1, S3=A2
        # LocB (ploidy 2): S1=B1, S2=B2, S3=uncalled
        df = pd.DataFrame([
            {"Seq_ID": "S1", "LocA_Hap_1": "X", "LocA_Hap_2": "",  "LocB_Hap_1": "X", "LocB_Hap_2": ""},
            {"Seq_ID": "S2", "LocA_Hap_1": "X", "LocA_Hap_2": "",  "LocB_Hap_1": "",  "LocB_Hap_2": "X"},
            {"Seq_ID": "S3", "LocA_Hap_1": "",  "LocA_Hap_2": "X", "LocB_Hap_1": "",  "LocB_Hap_2": ""},
        ])
        engine = PyEukDistanceEngine(min_completeness=0.0, ploidy={"LocA": 1, "LocB": 2})
        mat = engine.compute_revised_wibs_matrix(df)

        expected_d12 = 0.24264068711928537
        expected_d13 = 1.0
        expected_d23 = 1.0

        # Exact distance pin
        np.testing.assert_allclose(mat.loc["S1", "S2"], expected_d12, rtol=1e-5)
        np.testing.assert_allclose(mat.loc["S1", "S3"], expected_d13, rtol=1e-5)
        np.testing.assert_allclose(mat.loc["S2", "S3"], expected_d23, rtol=1e-5)

        # Assert mutant outputs diverge from true value
        # Mutant 1: Inverted weights sqrt(p(1-p))
        w_inv = np.array([np.sqrt(2/9), np.sqrt(2/9), np.sqrt(1/4), np.sqrt(1/4)])
        d12_mut1 = (0.5*w_inv[2] + 0.5*w_inv[3]) / (2*w_inv[0] + 2*w_inv[2])
        self.assertNotAlmostEqual(float(mat.loc["S1", "S2"]), d12_mut1, places=3)

        # Mutant 2: Uniform weights (w=1)
        d12_mut2 = 1.0 / 4.0
        self.assertNotAlmostEqual(float(mat.loc["S1", "S2"]), d12_mut2, places=3)

        # Mutant 3: Mask disabled (uncalled counted as mismatches)
        d13_mut3 = (2*2.12132034 + 1.0) / (2*2.12132034 + 4.0)
        self.assertNotAlmostEqual(float(mat.loc["S1", "S3"]), d13_mut3, places=2)

    def test_bayesian_plucinski_distance_pin_and_mutant_kill(self):
        """
        Pins exact Plucinski Bayesian distance on a 3-specimen x 2-locus fixture and explicitly asserts
        that reversing the posterior mapping (d_b = 2*P[2] + 1*P[1] + 0*P[0]) fails the test suite.
        """
        df = pd.DataFrame([
            {"Seq_ID": "S1", "Nu_1_Hap_1": "X", "Nu_1_Hap_2": "", "Nu_2_Hap_1": "X", "Nu_2_Hap_2": ""},
            {"Seq_ID": "S2", "Nu_1_Hap_1": "X", "Nu_1_Hap_2": "", "Nu_2_Hap_1": "X", "Nu_2_Hap_2": ""},
            {"Seq_ID": "S3", "Nu_1_Hap_1": "",  "Nu_1_Hap_2": "X", "Nu_2_Hap_1": "",  "Nu_2_Hap_2": "X"},
        ])
        engine = PyEukDistanceEngine(min_completeness=0.0)
        ids, loci_data, ploidy = engine._extract_locus_data(df)
        d_bayes = engine.compute_bayesian_distance(ids, loci_data, ploidy)

        # Exact pinned Bayesian distance values kill reversed mapping (which produces 1.227273)
        np.testing.assert_allclose(d_bayes[0, 1], 0.772727, rtol=1e-4)
        np.testing.assert_allclose(d_bayes[0, 2], 1.552873, rtol=1e-4)
        self.assertLess(d_bayes[0, 1], d_bayes[0, 2])

    def test_ploidy_differentiation(self):
        """Tests that ploidy=1, ploidy=2, and locus-specific ploidy produce distinct, mathematically valid matrices."""
        df = pd.DataFrame([
            {"Seq_ID": "S1", "Nu_L1_Hap_1": "X", "Nu_L1_Hap_2": "",  "Mt_L2_Hap_1": "X", "Mt_L2_Hap_2": ""},
            {"Seq_ID": "S2", "Nu_L1_Hap_1": "",  "Nu_L1_Hap_2": "X", "Mt_L2_Hap_1": "X", "Mt_L2_Hap_2": ""},
            {"Seq_ID": "S3", "Nu_L1_Hap_1": "X", "Nu_L1_Hap_2": "X", "Mt_L2_Hap_1": "",  "Mt_L2_Hap_2": "X"},
        ])
        eng_hap = PyEukDistanceEngine(min_completeness=0.0, ploidy=1)
        eng_dip = PyEukDistanceEngine(min_completeness=0.0, ploidy=2)
        eng_def = PyEukDistanceEngine(min_completeness=0.0, ploidy=None)

        m_hap = eng_hap.compute_revised_wibs_matrix(df)
        m_dip = eng_dip.compute_revised_wibs_matrix(df)
        m_def = eng_def.compute_revised_wibs_matrix(df)

        self.assertFalse(np.allclose(m_hap.values, m_dip.values))
        self.assertFalse(np.allclose(m_hap.values, m_def.values))
        self.assertFalse(np.allclose(m_dip.values, m_def.values))

    def test_gram_matrix_psd_spectrum(self):
        """Verifies PSD property (lambda_min >= -1e-12) and positive spectrum on both wIBS and Ensemble matrices."""
        df = pd.DataFrame([
            {"Seq_ID": f"S{i}", "L1_Hap_1": "X" if i%2==0 else "", "L1_Hap_2": "X" if i%2!=0 else "",
             "L2_Hap_1": "X" if i%3==0 else "", "L2_Hap_2": "X" if i%3!=0 else ""}
            for i in range(10)
        ])
        engine = PyEukDistanceEngine(min_completeness=0.0)
        wibs_mat = engine.compute_revised_wibs_matrix(df)
        ens_mat = engine.compute_ensemble_matrix(df)

        for name, mat in [("wIBS", wibs_mat), ("Ensemble", ens_mat)]:
            D = mat.values
            n = len(D)
            H = np.eye(n) - np.ones((n, n))/n
            G = -0.5 * H @ (D**2) @ H
            evals = np.linalg.eigvalsh(G)
            self.assertGreaterEqual(float(np.min(evals)), -1e-12, f"{name} failed PSD lower bound")
            self.assertGreater(float(np.max(evals)), 0.05, f"{name} failed positive spectrum test")

    def test_parse_locus_name_edge_cases(self):
        """Tests that locus name parsing handles diverse delimiter formats without false splits."""
        from cyclospora_pyeuk.distance_engine import parse_locus_name
        cases = {
            "Cp_HSP70_a1": "Cp_HSP70_a1",
            "Cp_HSP90_H01_508B": "Cp_HSP90",
            "gp60_L752bp.H01_9180": "gp60",
            "Nu_378_PART_A_Hap_4": "Nu_378_PART_A",
            "Mt_MSR_PART_A_Hap_1": "Mt_MSR_PART_A",
            "Mt_Cmt_Junction_PART_A_Hap_1": "Mt_Cmt",
            "tpi_L456bp.H02_58FE": "tpi",
            "18S_L458bp.H01_8F96": "18S",
            "COWP_L483bp.H01_462E": "COWP"
        }
        for col, expected in cases.items():
            self.assertEqual(parse_locus_name(col), expected)

    def test_mad_threshold_calibration_math(self):
        """Verifies exact Median + 3*MAD threshold formula calibration."""
        finder = CyclosporaClusterFinder(robust=True)
        gold_df = pd.DataFrame({
            "Seq_ID": ["A", "B", "C"],
            "True_Cluster": [1, 1, 1]
        })
        dist_df = pd.DataFrame([
            [0.0, 0.1, 0.3],
            [0.1, 0.0, 0.2],
            [0.3, 0.2, 0.0]
        ], index=["A", "B", "C"], columns=["A", "B", "C"])
        
        thresh = finder.calibrate_gold_standard_threshold(dist_df, gold_df)
        np.testing.assert_allclose(thresh, 0.64478, rtol=1e-4)

    def test_snp_weighted_wibs_with_real_sequences(self):
        """Verifies that SNP-weighted wIBS computes exact alignment Hamming distances."""
        seq1 = "A" * 100
        seq2 = "A" * 99 + "C"
        seq3 = "A" * 90 + "C" * 10

        seq_map = {
            "L1_H01": seq1,
            "L1_H02": seq2,
            "L1_H03": seq3
        }
        df = pd.DataFrame([
            {"Seq_ID": "S1", "L1_H01": "X", "L1_H02": "",  "L1_H03": ""},
            {"Seq_ID": "S2", "L1_H01": "",  "L1_H02": "X", "L1_H03": ""},
            {"Seq_ID": "S3", "L1_H01": "",  "L1_H02": "",  "L1_H03": "X"},
        ])
        engine = PyEukDistanceEngine(min_completeness=0.0)
        snp_mat = engine.compute_snp_weighted_wibs_matrix(df, sequence_map=seq_map)
        self.assertLess(snp_mat.loc["S1", "S2"], snp_mat.loc["S1", "S3"])

    def test_cryptosporidium_pipeline(self):
        """Tests that PyEuk runs reference-free de novo typing and outbreak clustering on Cryptosporidium."""
        from cyclospora_pyeuk.haplotype_sheet import learn_de_novo_haplotypes
        from cyclospora_pyeuk.distance_engine import PyEukDistanceEngine
        from cyclospora_pyeuk.clustering import CyclosporaClusterFinder
        from sklearn.metrics import adjusted_rand_score

        sheet_dn, learned = learn_de_novo_haplotypes("example_data/cryptosporidium/cohort_contigs.fasta")
        self.assertEqual(len(sheet_dn), 24)
        self.assertGreater(len(learned), 0)

        engine = PyEukDistanceEngine(min_completeness=0.0, ploidy=1)
        wibs = engine.compute_revised_wibs_matrix(sheet_dn)

        finder = CyclosporaClusterFinder()
        clusters_df, k, _ = finder.find_clusters(wibs, gold_file_path="example_data/cryptosporidium/gold_clusters.tsv")
        
        gold_df = pd.read_csv("example_data/cryptosporidium/gold_clusters.tsv", sep="\t")
        truth_map = dict(zip(gold_df["Seq_ID"], gold_df["True_Cluster"]))
        eval_df = clusters_df[clusters_df["Assigned_cluster"] != -1].copy()
        eval_df["True_Cluster"] = eval_df["Seq_ID"].map(truth_map)
        
        ari = adjusted_rand_score(eval_df["True_Cluster"], eval_df["Assigned_cluster"])
        self.assertGreaterEqual(ari, 0.40)

    def test_giardia_pipeline(self):
        """Tests that PyEuk runs reference-free de novo typing and outbreak clustering on Giardia."""
        from cyclospora_pyeuk.haplotype_sheet import learn_de_novo_haplotypes
        from cyclospora_pyeuk.distance_engine import PyEukDistanceEngine
        from cyclospora_pyeuk.clustering import CyclosporaClusterFinder
        from sklearn.metrics import adjusted_rand_score

        sheet_dn, learned = learn_de_novo_haplotypes("example_data/giardia/cohort_contigs.fasta")
        self.assertEqual(len(sheet_dn), 20)
        self.assertGreater(len(learned), 0)

        engine = PyEukDistanceEngine(min_completeness=0.0, ploidy=2)
        wibs = engine.compute_revised_wibs_matrix(sheet_dn)

        finder = CyclosporaClusterFinder()
        clusters_df, k, _ = finder.find_clusters(wibs, gold_file_path="example_data/giardia/gold_clusters.tsv")
        
        gold_df = pd.read_csv("example_data/giardia/gold_clusters.tsv", sep="\t")
        truth_map = dict(zip(gold_df["Seq_ID"], gold_df["True_Cluster"]))
        eval_df = clusters_df[clusters_df["Assigned_cluster"] != -1].copy()
        eval_df["True_Cluster"] = eval_df["Seq_ID"].map(truth_map)
        
        # Verify separation of Assemblage A and Assemblage B lineages
        ari = adjusted_rand_score(eval_df["True_Cluster"], eval_df["Assigned_cluster"])
        self.assertGreaterEqual(ari, 0.90)

    def test_unweighted_hamming_baseline_falls_below_pyeuk(self):
        """
        Verifies that on the non-trivial Cryptosporidium panel with shared background alleles:
        1. No single locus alone achieves ARI == 1.0 (all single-locus ARIs < 0.50).
        2. Under Complete linkage, Plain unweighted Hamming drops to ARI = 0.6503 while PyEuk KING-wIBS achieves ARI = 1.0000.
        3. Under Ward linkage, PyEuk KING-wIBS achieves ARI = 0.8836.
        """
        from cyclospora_pyeuk.haplotype_sheet import learn_de_novo_haplotypes
        from cyclospora_pyeuk.distance_engine import PyEukDistanceEngine
        from scipy.spatial.distance import pdist, squareform
        from scipy.cluster.hierarchy import linkage, cut_tree
        from sklearn.metrics import adjusted_rand_score

        sheet_df, _ = learn_de_novo_haplotypes("example_data/cryptosporidium/cohort_contigs.fasta")
        gold_df = pd.read_csv("example_data/cryptosporidium/gold_clusters.tsv", sep="\t")
        truth_map = dict(zip(gold_df["Seq_ID"], gold_df["True_Cluster"]))
        truth = [truth_map[s] for s in sheet_df["Seq_ID"]]

        marker_cols = [c for c in sheet_df.columns if c != "Seq_ID"]

        # 1. Assert no single locus separates the cohort
        for locus in ["18S", "HSP70", "COWP", "gp60"]:
            l_cols = [c for c in marker_cols if c.startswith(locus)]
            X_l = (sheet_df[l_cols].values == "X").astype(float)
            _, inv = np.unique(X_l, axis=0, return_inverse=True)
            ari_l = adjusted_rand_score(truth, inv)
            self.assertLess(ari_l, 0.50, f"Locus {locus} has ARI {ari_l} >= 0.50")

        # 2. Complete linkage: Plain unweighted Hamming = 0.6503 vs PyEuk KING-wIBS = 1.0000
        X = (sheet_df[marker_cols].values == "X").astype(float)
        d_ham = squareform(pdist(X, metric="hamming"))
        c_ham_comp = cut_tree(linkage(squareform(d_ham), method="complete"), n_clusters=4).ravel()
        ari_ham_comp = adjusted_rand_score(truth, c_ham_comp)
        self.assertAlmostEqual(ari_ham_comp, 0.6503, places=2)

        eng = PyEukDistanceEngine(min_completeness=0.0, ploidy=1)
        wibs_mat = eng.compute_revised_wibs_matrix(sheet_df)
        c_wibs_comp = cut_tree(linkage(squareform(wibs_mat.values), method="complete"), n_clusters=4).ravel()
        ari_wibs_comp = adjusted_rand_score(truth, c_wibs_comp)
        self.assertEqual(ari_wibs_comp, 1.0000)

        # 3. Ward linkage: PyEuk KING-wIBS achieves ARI = 0.8836
        c_wibs_ward = cut_tree(linkage(squareform(wibs_mat.values), method="ward"), n_clusters=4).ravel()
        ari_wibs_ward = adjusted_rand_score(truth, c_wibs_ward)
        self.assertAlmostEqual(ari_wibs_ward, 0.8836, places=2)

    def test_disjoint_locus_coverage_receives_max_distance_and_never_zero(self):
        """
        Adversarial test for Finding 1:
        Specimens with complementary locus dropouts (0 mutually called loci) must receive
        maximum distance (1.0), NEVER distance 0.0 (genetic identity).
        """
        from cyclospora_pyeuk.distance_engine import PyEukDistanceEngine
        
        # Specimen 1: only Locus 1 called
        # Specimen 2: only Locus 2 called (0 overlap with Specimen 1)
        # Specimen 3: both Locus 1 and Locus 2 called
        df_disjoint = pd.DataFrame([
            {"Seq_ID": "SPEC_01", "L1_H01": "X", "L1_H02": "", "L2_H01": "", "L2_H02": ""},
            {"Seq_ID": "SPEC_02", "L1_H01": "", "L1_H02": "", "L2_H01": "X", "L2_H02": ""},
            {"Seq_ID": "SPEC_03", "L1_H01": "X", "L1_H02": "", "L2_H01": "X", "L2_H02": ""}
        ]).fillna("")

        eng = PyEukDistanceEngine(min_completeness=0.0, ploidy=1)
        
        # wIBS distance
        wibs_mat = eng.compute_revised_wibs_matrix(df_disjoint)
        self.assertAlmostEqual(wibs_mat.loc["SPEC_01", "SPEC_02"], 1.0, places=5)
        self.assertNotEqual(wibs_mat.loc["SPEC_01", "SPEC_02"], 0.0)

        # SNP-wIBS distance
        snp_mat = eng.compute_snp_weighted_wibs_matrix(df_disjoint)
        self.assertAlmostEqual(snp_mat.loc["SPEC_01", "SPEC_02"], 1.0, places=5)
        self.assertNotEqual(snp_mat.loc["SPEC_01", "SPEC_02"], 0.0)

    def test_background_specimen_separation_and_collision_guard(self):
        """
        Adversarial test for Finding 2:
        Haplotype sheet generation must prevent sample ID collisions between specimen and background dirs
        and must not silently emit background controls as specimen case rows.
        """
        from cyclospora_pyeuk.haplotype_sheet import generate_haplotype_sheet
        import tempfile
        
        with tempfile.TemporaryDirectory() as spec_dir, tempfile.TemporaryDirectory() as bg_dir:
            # Create a specimen file
            with open(os.path.join(spec_dir, "SAMPLE_001.txt"), "w") as f:
                f.write("L1_H01\t100.0\n")
            
            # Create duplicate ID in background directory
            with open(os.path.join(bg_dir, "SAMPLE_001.txt"), "w") as f:
                f.write("L1_H02\t100.0\n")
                
            # Must raise ValueError on ID collision
            with self.assertRaises(ValueError):
                generate_haplotype_sheet(specimen_dir=spec_dir, background_dir=bg_dir)

    def test_locus_completeness_invariant_to_unused_allele_columns(self):
        """
        Adversarial test for gpt-5.6-sol finding:
        Specimen locus completeness must measure the fraction of called loci, NOT allele-column density.
        Adding 100 unused allele columns to the sheet must not alter locus completeness or filter out valid specimens.
        """
        from cyclospora_pyeuk.distance_engine import PyEukDistanceEngine
        
        # Specimen called at all 4 loci (1 allele each out of 4 columns)
        df_base = pd.DataFrame([
            {"Seq_ID": "SPEC_01", "18S_H01": "X", "HSP70_H01": "X", "COWP_H01": "X", "gp60_H01": "X"}
        ])
        
        # Expanded schema with 100 unused population allele columns across the 4 loci
        df_expanded = df_base.copy()
        for loc in ["18S", "HSP70", "COWP", "gp60"]:
            for a in range(2, 27):
                df_expanded[f"{loc}_H{a:02d}"] = ""
                
        # Total columns = 100 + 1 (Seq_ID). Allele column density would be 4 / 100 = 0.04 (< 0.50)
        # But true locus completeness is 4 / 4 = 1.0 (>= 0.50)
        engine = PyEukDistanceEngine(min_completeness=0.50)
        clean_df = engine.process_haplotype_sheet(df_expanded)
        self.assertEqual(len(clean_df), 1)
        self.assertEqual(clean_df.iloc[0]["Seq_ID"], "SPEC_01")

    def test_fetch_bioproject_cohort_strictly_validates_accessions(self):
        """
        Adversarial test for Finding 4:
        efetch_fasta_batch must fail loudly with RuntimeError on missing accessions rather than substituting synthetic homopolymers.
        """
        from scripts.fetch_bioproject_cohort import efetch_fasta_batch
        
        # Test with a non-existent bogus accession
        with self.assertRaises(RuntimeError):
            efetch_fasta_batch(["NON_EXISTENT_ACCESSION_999999"])

    def test_wibs_distance_invariant_to_unobserved_allele_columns(self):
        """
        Adversarial test for Issue #8:
        Unobserved allele columns (p == 0.0) and fixed columns (p == 1.0) must not inflate the wIBS denominator
        with dead weight. Adding unobserved reference panel columns must leave pairwise distances bit-identical.
        """
        from cyclospora_pyeuk.distance_engine import PyEukDistanceEngine
        
        df_base = pd.DataFrame([
            {"Seq_ID": "SPEC_01", "L1_H01": "X", "L1_H02": "", "L2_H01": "X", "L2_H02": ""},
            {"Seq_ID": "SPEC_02", "L1_H01": "", "L1_H02": "X", "L2_H01": "X", "L2_H02": ""},
            {"Seq_ID": "SPEC_03", "L1_H01": "X", "L1_H02": "", "L2_H01": "", "L2_H02": "X"},
            {"Seq_ID": "SPEC_04", "L1_H01": "", "L1_H02": "X", "L2_H01": "", "L2_H02": "X"},
        ]).fillna("")
        
        eng = PyEukDistanceEngine(min_completeness=0.10, ploidy=1)
        D_base = eng.compute_revised_wibs_matrix(df_base)
        
        # Add 100 unobserved columns (simulating a reference catalog with many unobserved alleles)
        df_expanded = df_base.copy()
        for i in range(3, 53):
            df_expanded[f"L1_H{i:02d}"] = ""
            df_expanded[f"L2_H{i:02d}"] = ""
            
        D_expanded = eng.compute_revised_wibs_matrix(df_expanded)
        
        # Pairwise distance matrix must be identical
        np.testing.assert_allclose(D_base.values, D_expanded.values, atol=1e-7)
        self.assertGreater(D_base.values.max(), 0.10, "Genetic distance should not be compressed to near-zero.")


if __name__ == "__main__":
    unittest.main()




