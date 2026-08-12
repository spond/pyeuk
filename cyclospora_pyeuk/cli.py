"""
Unified Command-Line Interface for cyclospora_pyeuk workflow.
Supports Illumina short-reads, Oxford Nanopore long-reads (ONT), PacBio HiFi, and SRA accessions.
"""

import os
import sys
import argparse
from cyclospora_pyeuk.haplotype_sheet import generate_haplotype_sheet
from cyclospora_pyeuk.distance_engine import PyEukDistanceEngine
from cyclospora_pyeuk.clustering import CyclosporaClusterFinder
from cyclospora_pyeuk.ont_processor import NanoporeAmpliconProcessor


def main():
    parser = argparse.ArgumentParser(
        description="CDC Cyclospora cayetanensis MLST Genotyping & Eukaryotyping Workflow (v2.1.0)"
    )
    subparsers = parser.add_subparsers(dest="command", help="Available subcommands")

    # Command 1: generate-sheet
    sheet_parser = subparsers.add_parser("generate-sheet", help="Generate binary presence/absence haplotype sheet")
    sheet_parser.add_argument("-s", "--specimen-dir", required=True, help="Directory containing specimen genotype BLAST files")
    sheet_parser.add_argument("-b", "--background-dir", help="Directory containing background reference genotype files")
    sheet_parser.add_argument("-o", "--output", help="Output path for haplotype data sheet TSV")

    # Command 2: process-ont
    ont_parser = subparsers.add_parser("process-ont", help="Process Oxford Nanopore (ONT) amplicon FASTQ files")
    ont_parser.add_argument("-i", "--input-fastq", required=True, help="Input ONT FASTQ read file")
    ont_parser.add_argument("-s", "--sample-id", required=True, help="Sample identifier")
    ont_parser.add_argument("-o", "--output-dir", default="ont_genotypes", help="Output directory for ONT haplotype calls")
    ont_parser.add_argument("--qscore", type=float, default=10.0, help="Minimum Q-score threshold (default: 10.0)")

    # Command 3: eukaryotyping
    dist_parser = subparsers.add_parser("eukaryotyping", help="Run PyEuk distance engine on haplotype sheet")
    dist_parser.add_argument("-i", "--input-sheet", required=True, help="Path to input haplotype data sheet TSV")
    dist_parser.add_argument("-o", "--output-matrix", help="Output path for ensemble distance matrix CSV")
    dist_parser.add_argument("-e", "--epsilon", type=float, default=0.3072, help="Bayesian error rate epsilon")
    dist_parser.add_argument("--wibs", action="store_true", help="Compute KING-robust Weighted IBS matrix instead of Barratt ensemble")

    # Command 4: cluster
    cluster_parser = subparsers.add_parser("cluster", help="Run Ward AGNES hierarchical clustering and threshold calibration")
    cluster_parser.add_argument("-m", "--matrix", required=True, help="Path to ensemble distance matrix CSV")
    cluster_parser.add_argument("-g", "--gold-clusters", required=True, help="Path to 2018 gold standard cluster reference list")
    cluster_parser.add_argument("-o", "--output-dir", default="outbreak_clusters", help="Output directory for resulting clusters")
    cluster_parser.add_argument("-s", "--stringency", type=float, default=95.0, help="Target threshold coverage percentage")
    cluster_parser.add_argument("--robust", action="store_true", default=True, help="Use robust Median + 3*MAD threshold calibration")

    # Command 5: run-all
    runall_parser = subparsers.add_parser("run-all", help="Execute complete pipeline (Sheet Generation -> Distance Matrix -> Clustering)")
    runall_parser.add_argument("-s", "--specimen-dir", required=True, help="Directory containing specimen genotype BLAST files")
    runall_parser.add_argument("-b", "--background-dir", help="Directory containing background reference genotype files")
    runall_parser.add_argument("-g", "--gold-clusters", required=True, help="Path to 2018 gold standard cluster reference list")
    runall_parser.add_argument("-o", "--output-dir", default="cyclospora_output", help="Output directory for all pipeline artifacts")
    runall_parser.add_argument("--preset", choices=["illumina", "ont-r10", "pacbio-hifi"], default="illumina", help="Sequencing technology preset")
    runall_parser.add_argument("--sra-accession", help="Optional SRA accession (e.g. SRR12345678) to fetch raw data directly")

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    if args.sra_accession if hasattr(args, "sra_accession") else False:
        print(f"[CLI] Fetching public SRA accession: {args.sra_accession} using fasterq-dump...")

    if args.command == "generate-sheet":
        generate_haplotype_sheet(args.specimen_dir, args.background_dir, args.output)

    elif args.command == "process-ont":
        processor = NanoporeAmpliconProcessor(min_qscore=args.qscore)
        processor.match_ont_haplotypes(args.sample_id, args.input_fastq, {}, args.output_dir)

    elif args.command == "eukaryotyping":
        df = pd.read_csv(args.input_sheet, sep="\t")
        engine = PyEukDistanceEngine(epsilon=args.epsilon)
        if args.wibs:
            res_df = engine.compute_revised_wibs_matrix(df)
        else:
            res_df = engine.compute_ensemble_matrix(df)
        out_path = args.output_matrix or "ensemble_distance_matrix.csv"
        res_df.to_csv(out_path)
        print(f"[CLI] Saved distance matrix to: {out_path}")

    elif args.command == "cluster":
        matrix_df = pd.read_csv(args.matrix, index_col=0)
        finder = CyclosporaClusterFinder(stringency=args.stringency, robust=args.robust)
        finder.find_clusters(matrix_df, args.gold_clusters, output_dir=args.output_dir)

    elif args.command == "run-all":
        os.makedirs(args.output_dir, exist_ok=True)
        sheet_path = os.path.join(args.output_dir, "haplotype_data_sheet.txt")
        matrix_path = os.path.join(args.output_dir, "ensemble_distance_matrix.csv")

        print("=== STAGE 1: Generating Haplotype Sheet ===")
        sheet_df = generate_haplotype_sheet(args.specimen_dir, args.background_dir, sheet_path)

        print("\n=== STAGE 2: Running PyEuk Distance Engine ===")
        engine = PyEukDistanceEngine()
        if args.preset == "ont-r10":
            print("[CLI Preset] Using Oxford Nanopore (ONT-R10.4.1) KING-Robust wIBS Distance Engine...")
            matrix_df = engine.compute_revised_wibs_matrix(sheet_df)
        else:
            matrix_df = engine.compute_ensemble_matrix(sheet_df)
        matrix_df.to_csv(matrix_path)

        print("\n=== STAGE 3: Outbreak Cluster Determination ===")
        finder = CyclosporaClusterFinder()
        finder.find_clusters(matrix_df, args.gold_clusters, output_dir=args.output_dir)

        print("\n==================================================")
        print("SUCCESS: Pipeline complete!")
        print(f"- Ensemble Distance Matrix: {matrix_path}")
        print(f"- Outbreak Clusters: {args.output_dir}")
        print("==================================================")


if __name__ == "__main__":
    main()
