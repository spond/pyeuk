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


if __name__ == "__main__":
    unittest.main()




