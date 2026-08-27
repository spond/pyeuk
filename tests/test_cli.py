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

    def test_cli_process_ont_with_reference(self):
        # Create test FASTQ and test reference FASTA
        test_ref_fa = os.path.join(self.test_out, "ref.fasta")
        with open(test_ref_fa, "w") as f:
            f.write(">L1_H01_TEST\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n")

        test_fastq = os.path.join(self.test_out, "reads.fastq")
        read_seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" * 5
        read_qual = "I" * len(read_seq)
        with open(test_fastq, "w") as f:
            for i in range(10):
                f.write(f"@read_{i}\n{read_seq}\n+\n{read_qual}\n")

        ont_out = os.path.join(self.test_out, "ont_out")
        sys.argv = [
            "pyeuk", "process-ont",
            "-i", test_fastq,
            "-s", "SAMPLE_001",
            "-r", test_ref_fa,
            "-o", ont_out,
            "--qscore", "10.0"
        ]
        main()
        call_file = os.path.join(ont_out, "SAMPLE_001")
        self.assertTrue(os.path.exists(call_file))
        import pandas as pd
        df = pd.read_csv(call_file, sep="\t", header=None)
        self.assertFalse(df.empty)
        self.assertEqual(df.iloc[0, 0], "L1_H01_TEST")

    def test_cli_weight_mode_and_no_psd_options(self):
        # Test eukaryotyping with --weight-mode, --min-maf, and --no-psd
        sheet_path = os.path.join(self.test_out, "sheet_test.txt")
        sys.argv = [
            "pyeuk", "generate-sheet",
            "-s", "example_data/specimens",
            "-o", sheet_path
        ]
        main()

        out_matrix_raw = os.path.join(self.test_out, "dist_raw.csv")
        sys.argv = [
            "pyeuk", "eukaryotyping",
            "-i", sheet_path,
            "-o", out_matrix_raw,
            "--metric", "wibs",
            "--weight-mode", "heterozygosity",
            "--min-maf", "0.05",
            "--no-psd",
            "--ploidy", "2"
        ]
        main()
        self.assertTrue(os.path.exists(out_matrix_raw))

        out_matrix_king = os.path.join(self.test_out, "dist_king.csv")
        sys.argv = [
            "pyeuk", "eukaryotyping",
            "-i", sheet_path,
            "-o", out_matrix_king,
            "--metric", "wibs",
            "--weight-mode", "king",
            "--ploidy", "2"
        ]
        main()
        self.assertTrue(os.path.exists(out_matrix_king))

    def test_cli_clustering_k_min_and_relative_gap_floor(self):
        sheet_path = os.path.join(self.test_out, "sheet_cluster_test.txt")
        sys.argv = [
            "pyeuk", "generate-sheet",
            "-s", "example_data/specimens",
            "-o", sheet_path
        ]
        main()

        out_matrix = os.path.join(self.test_out, "dist_matrix_k.csv")
        sys.argv = [
            "pyeuk", "eukaryotyping",
            "-i", sheet_path,
            "-o", out_matrix,
            "--metric", "wibs"
        ]
        main()

        # Run cluster with --k-min, --k-max, and --relative-gap-floor
        cluster_out = os.path.join(self.test_out, "clusters_k_test")
        sys.argv = [
            "pyeuk", "cluster",
            "-m", out_matrix,
            "--k-min", "2",
            "--k-max", "10",
            "--relative-gap-floor", "0.10",
            "-o", cluster_out
        ]
        main()
        self.assertTrue(os.path.exists(cluster_out))

    def test_cli_clustering_distance_cut_and_single_linkage(self):
        sheet_path = os.path.join(self.test_out, "sheet_dist_test.txt")
        sys.argv = [
            "pyeuk", "generate-sheet",
            "-s", "example_data/specimens",
            "-o", sheet_path
        ]
        main()

        out_matrix = os.path.join(self.test_out, "dist_matrix_d.csv")
        sys.argv = [
            "pyeuk", "eukaryotyping",
            "-i", sheet_path,
            "-o", out_matrix,
            "--metric", "wibs"
        ]
        main()

        cluster_out = os.path.join(self.test_out, "clusters_dist_test")
        sys.argv = [
            "pyeuk", "cluster",
            "-m", out_matrix,
            "--cut", "distance",
            "--linkage-method", "single",
            "--linkage-threshold", "0.15",
            "-o", cluster_out
        ]
        main()
        self.assertTrue(os.path.exists(cluster_out))

    def test_cli_run_all_distance_cut(self):
        out_dir = os.path.join(self.test_out, "run_all_dist")
        sys.argv = [
            "pyeuk", "run-all",
            "-a", "example_data/cryptosporidium/cohort_contigs.fasta",
            "--de-novo",
            "--cut", "distance",
            "--linkage-method", "ward",
            "--linkage-threshold", "0.10",
            "--ploidy", "1",
            "--min-completeness", "0.0",
            "-o", out_dir
        ]
        main()
        self.assertTrue(os.path.exists(os.path.join(out_dir, "haplotype_data_sheet.txt")))
        self.assertTrue(os.path.exists(os.path.join(out_dir, "ensemble_distance_matrix.csv")))

    def test_cli_amplicon_build_sheet_subcommand(self):
        calls_dir = os.path.join(self.test_out, "mock_calls")
        os.makedirs(calls_dir, exist_ok=True)
        call_tsv = os.path.join(calls_dir, "SAMPLE_A.tsv")
        with open(call_tsv, "w") as f:
            f.write("specimen\tlocus\twindow\tstart\tend\thaplotype\treads\tfreq\tspanning\n")
            f.write("SAMPLE_A\tL1\tL1_W0001\t1\t100\t=\t50\t1.0\t50\n")

        sheet_out = os.path.join(self.test_out, "amplicon_sheet_out")
        sys.argv = [
            "pyeuk", "build-sheet",
            calls_dir,
            sheet_out
        ]
        main()
        self.assertTrue(os.path.exists(os.path.join(sheet_out, "sheet.tsv")))
        self.assertTrue(os.path.exists(os.path.join(sheet_out, "haplotype_map.tsv")))
        self.assertTrue(os.path.exists(os.path.join(sheet_out, "calls_long.tsv")))


if __name__ == "__main__":
    unittest.main()
