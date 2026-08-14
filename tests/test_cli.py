"""
Functional tests for PyEuk CLI entrypoints and arguments.
"""

import os
import sys
import shutil
import unittest
from cyclospora_pyeuk.cli import main

class TestPyEukCLI(unittest.TestCase):

    def setUp(self):
        self.test_out = "tests/test_cli_output"
        os.makedirs(self.test_out, exist_ok=True)

    def tearDown(self):
        if os.path.exists(self.test_out):
            shutil.rmtree(self.test_out)

    def test_cli_eukaryotyping_wibs_and_ploidy(self):
        sheet_path = "example_data/specimens/NC_001_18_BLAST_results.txt" # We will test on generated sheet
        # Generate sheet first
        sys.argv = [
            "pyeuk", "generate-sheet",
            "-s", "example_data/specimens",
            "-o", os.path.join(self.test_out, "sheet.txt")
        ]
        main()
        self.assertTrue(os.path.exists(os.path.join(self.test_out, "sheet.txt")))

        # Run eukaryotyping with wibs and ploidy
        out_matrix = os.path.join(self.test_out, "dist_matrix.csv")
        sys.argv = [
            "pyeuk", "eukaryotyping",
            "-i", os.path.join(self.test_out, "sheet.txt"),
            "-o", out_matrix,
            "--wibs",
            "--ploidy", "2"
        ]
        main()
        self.assertTrue(os.path.exists(out_matrix))

        # Run cluster
        sys.argv = [
            "pyeuk", "cluster",
            "-m", out_matrix,
            "-g", "example_data/gold_clusters.tsv",
            "-o", os.path.join(self.test_out, "clusters")
        ]
        main()
        self.assertTrue(os.path.exists(os.path.join(self.test_out, "clusters")))

    def test_cli_run_all_de_novo(self):
        out_dir = os.path.join(self.test_out, "run_all_crypto")
        sys.argv = [
            "pyeuk", "run-all",
            "-a", "example_data/cryptosporidium/cohort_contigs.fasta",
            "--de-novo",
            "-g", "example_data/cryptosporidium/gold_clusters.tsv",
            "--ploidy", "1",
            "--min-completeness", "0.0",
            "-o", out_dir
        ]
        main()
        self.assertTrue(os.path.exists(os.path.join(out_dir, "haplotype_data_sheet.txt")))
        self.assertTrue(os.path.exists(os.path.join(out_dir, "learned_refs.fasta")))
        self.assertTrue(os.path.exists(os.path.join(out_dir, "ensemble_distance_matrix.csv")))

if __name__ == "__main__":
    unittest.main()
