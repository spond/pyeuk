"""
Unified Command-Line Interface for cyclospora_pyeuk workflow.
Supports Illumina short-reads, Oxford Nanopore long-reads (ONT), PacBio HiFi, and SRA accessions.
"""

import os
import sys
import argparse
import subprocess
import zipfile
from .haplotype_sheet import generate_haplotype_sheet
from .distance_engine import PyEukDistanceEngine
from .clustering import CyclosporaClusterFinder
from .ont_processor import NanoporeAmpliconProcessor


def fetch_test_data(target_dir: str = "./cdc_reference_data"):
    """
    Clones official CDC reference repository and automatically extracts benchmark test data.
    """
    os.makedirs(target_dir, exist_ok=True)
    repo_url = "https://github.com/Joel-Barratt/Complete-Cyclospora-typing-workflow.git"
    clone_dest = os.path.join(target_dir, "cdc_repo")

    if not os.path.exists(clone_dest):
        print(f"[FetchTestData] Cloning official CDC reference repository from {repo_url}...")
        res = subprocess.run(["git", "clone", "--depth", "1", repo_url, clone_dest])
        if res.returncode != 0:
            print("[FetchTestData Error] Failed to clone CDC repository. Please check your internet connection.")
            return

    specimens_out = os.path.join(target_dir, "specimens")
    os.makedirs(specimens_out, exist_ok=True)

    zip_found = None
    gold_found = None

    for root, _, files in os.walk(clone_dest):
        for fname in files:
            if fname == "SPECIMEN_GENOTYPES.zip":
                zip_found = os.path.join(root, fname)
            elif fname == "2018_gold_clusters.txt":
                gold_found = os.path.join(root, fname)

    if zip_found:
        print(f"[FetchTestData] Unpacking benchmark genotypes from {zip_found}...")
        with zipfile.ZipFile(zip_found, 'r') as zf:
            zf.extractall(specimens_out)
        print(f"[FetchTestData] Unpacked {len(os.listdir(specimens_out))} specimen genotype files into: {specimens_out}")

    gold_out = os.path.join(target_dir, "2018_gold_clusters.txt")
    if gold_found and not os.path.exists(gold_out):
        import shutil
        shutil.copy(gold_found, gold_out)

    print("\n==========================================================================")
    print("SUCCESS: CDC Benchmark Test Data Ready!")
    print(f"  • Specimen Genotypes Directory : {specimens_out}")
    print(f"  • Gold Standard Clusters File  : {gold_out}")
    print("\nQuick Run Example (Label-Free Unsupervised):")
    print(f"  cyclospora-typing run-all -s {specimens_out} -o ./outbreak_output")
    print("==========================================================================\n")


def main():
    parser = argparse.ArgumentParser(
        prog="pyeuk",
        description="PyEuk: High-Performance Molecular Typing, Genetic Distance Estimation, and Outbreak Clustering for Eukaryotic and Microbial Pathogens (Cyclospora, Cryptosporidium, MLST, cgMLST)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", help="Available sub-commands")

    # Command 0: fetch-test-data
    fetch_parser = subparsers.add_parser("fetch-test-data", help="Fetch validated public surveillance benchmark dataset")
    fetch_parser.add_argument("-o", "--output-dir", default="test_data", help="Directory to save test dataset")

    # Command 1: generate-sheet
    sheet_parser = subparsers.add_parser("generate-sheet", help="Generate binary presence/absence haplotype sheet (supports .zip, BLAST calls, or assembled FASTA)")
    sheet_parser.add_argument("-s", "--specimen-dir", help="Directory or .zip file containing specimen genotype BLAST files or FASTA contigs")
    sheet_parser.add_argument("-a", "--assembled-fasta", help="Directory or FASTA file containing externally assembled haplotype contigs")
    sheet_parser.add_argument("-b", "--background-dir", help="Directory or .zip file containing background reference genotype files")
    sheet_parser.add_argument("-o", "--output", help="Output path for haplotype data sheet TSV")
    sheet_parser.add_argument("--de-novo", action="store_true", help="Discover loci and haplotypes de novo without any reference database")

    # Command 2: process-ont
    ont_parser = subparsers.add_parser("process-ont", help="Process Oxford Nanopore (ONT) amplicon FASTQ files")
    ont_parser.add_argument("-i", "--input-fastq", required=True, help="Input ONT FASTQ read file")
    ont_parser.add_argument("-s", "--sample-id", required=True, help="Sample identifier")
    ont_parser.add_argument("-r", "--reference-fasta", help="Path to reference MLST FASTA database")
    ont_parser.add_argument("--de-novo", action="store_true", help="Perform de novo consensus haplotype discovery from reads")
    ont_parser.add_argument("-o", "--output-dir", default="ont_genotypes", help="Output directory for ONT haplotype calls")
    ont_parser.add_argument("--qscore", type=float, default=10.0, help="Minimum Q-score threshold (default: 10.0)")

    # Command 3: eukaryotyping
    dist_parser = subparsers.add_parser("eukaryotyping", help="Run PyEuk distance engine on haplotype sheet")
    dist_parser.add_argument("-i", "--input-sheet", required=True, help="Path to input haplotype data sheet TSV")
    dist_parser.add_argument("-o", "--output-matrix", help="Output path for ensemble distance matrix CSV")
    dist_parser.add_argument("-e", "--epsilon", type=float, default=0.3072, help="Bayesian error rate epsilon")
    dist_parser.add_argument("--wibs", action="store_true", help="Compute Weighted IBS matrix instead of Barratt ensemble")
    dist_parser.add_argument("--metric", choices=["wibs", "ensemble", "snp-wibs"], default=None, help="Distance metric to compute")
    dist_parser.add_argument("--ploidy", type=int, default=None, help="Organism base ploidy (e.g. 1 for haploid Cryptosporidium/bacteria, 2 for diploid)")
    dist_parser.add_argument("-f", "--fasta", help="Optional path to reference or learned FASTA file for sequence-weighted SNP-wIBS distance")
    dist_parser.add_argument("--min-completeness", type=float, default=0.10, help="Minimum locus completeness fraction (default: 0.10)")
    dist_parser.add_argument("--weight-mode", choices=["heterozygosity", "inverted-king", "king", "uniform"], default="heterozygosity", help="Allele weighting scheme for presence/absence indicator columns (default: heterozygosity)")
    dist_parser.add_argument("--min-maf", type=float, default=0.0, help="Minimum minor allele frequency threshold to filter rare/private singletons (default: 0.0)")
    dist_parser.add_argument("--no-psd", dest="project_psd", action="store_false", default=True, help="Skip Gram matrix PSD projection and return raw pairwise distances")

    # Command 4: cluster
    cluster_parser = subparsers.add_parser("cluster", help="Run Ward AGNES hierarchical clustering (unsupervised or supervised)")
    cluster_parser.add_argument("-m", "--matrix", required=True, help="Path to ensemble distance matrix CSV")
    cluster_parser.add_argument("-g", "--gold-clusters", required=False, default=None, help="Optional path to gold standard cluster reference list (for supervised mode)")
    cluster_parser.add_argument("-o", "--output-dir", default="outbreak_clusters", help="Output directory for resulting clusters")
    cluster_parser.add_argument("-s", "--stringency", type=float, default=95.0, help="Target threshold coverage percentage")
    cluster_parser.add_argument("--robust", action="store_true", default=True, help="Use robust Median + 3*MAD threshold calibration")
    cluster_parser.add_argument("--k-min", type=int, default=2, help="Minimum number of clusters to search (default: 2)")
    cluster_parser.add_argument("--k-max", type=int, default=50, help="Maximum number of clusters to search (default: 50)")
    cluster_parser.add_argument("--relative-gap-floor", type=float, default=0.2200, help="Minimum relative merge-height gap fraction of tree height required for unsupervised knee selection (default: 0.2200)")
    cluster_parser.add_argument(
        "--cut", choices=["count", "distance"], default="count",
        help="How to cut the tree. 'count' (default) chooses a number of clusters from the "
             "largest merge-height gap, subject to a minimum gap and a minimum cluster size -- "
             "the right question for a closed outbreak investigation. 'distance' cuts at a "
             "fixed dissimilarity instead, with no cluster count and no gap requirement, and "
             "returns specimens with no near neighbour as singletons -- the right question for "
             "surveillance, where most cases are unrelated. Choosing a count cannot represent a "
             "mostly-singleton structure and collapses to k=1 on such data.")
    cluster_parser.add_argument(
        "--linkage-method", choices=["ward", "single", "average", "complete"], default="ward",
        help="Linkage for the tree (default: ward). 'single' follows chains rather than "
             "compact blobs, which is what a transmission cluster is, and recovers the exact "
             "published cluster count on the P. vivax AmpliSeq cohort. It chains through noise "
             "at a loose threshold, so it is not the default.")
    cluster_parser.add_argument(
        "--linkage-threshold", type=float, default=None,
        help="Dissimilarity at which to cut in --cut distance mode. Omit to calibrate it from "
             "labelled pairs when --gold-clusters is given, or from the distance distribution "
             "otherwise; the provenance of the value is printed either way.")
    cluster_parser.add_argument(
        "--single-k", action="store_true", default=False,
        help="Report one cluster count instead of the default sweep diagnostic. The sweep "
             "reports the range of counts the data supports, how confident that is, and a "
             "per-branch confidence tree; a single k is only trustworthy when the count "
             "selectors agree, which the sweep tells you. Use this only when a downstream step "
             "strictly needs the legacy single-partition behaviour.")
    cluster_parser.add_argument(
        "--n-boot", type=int, default=200,
        help="Bootstrap resamples for branch support and stability in the sweep (default: 200).")

    # Commands 6-8: amplicon front end (BAMs -> haplotype sheet)
    # Each forwards straight to the module's own parser rather than restating its flags here.
    # Restating them is how the two copies drift, and the flags in question are the ones that
    # decide the result: window width, spanning floor, and the three gates.
    subparsers.add_parser(
        "define-windows", add_help=False,
        help="Choose analysis windows from the cohort's own reads (requires the amplicon extra)")
    subparsers.add_parser(
        "call-haplotypes", add_help=False,
        help="Read a haplotype off each single spanning read (requires the amplicon extra)")
    subparsers.add_parser(
        "build-sheet", add_help=False,
        help="Assemble the specimen x haplotype sheet from per-specimen calls")

    # Command 5: run-all
    runall_parser = subparsers.add_parser("run-all", help="Execute complete pipeline (Sheet Generation -> Distance Matrix -> Clustering)")
    runall_parser.add_argument("-s", "--specimen-dir", help="Directory or .zip file containing specimen genotype BLAST files or FASTA contigs")
    runall_parser.add_argument("-a", "--assembled-fasta", help="Directory or FASTA file containing externally assembled haplotype contigs")
    runall_parser.add_argument("-b", "--background-dir", help="Directory or .zip file containing background reference genotype files")
    runall_parser.add_argument("-g", "--gold-clusters", required=False, default=None, help="Optional path to gold standard cluster reference list (for supervised mode)")
    runall_parser.add_argument("-o", "--output-dir", default="outbreak_results", help="Output directory for all pipeline artifacts")
    runall_parser.add_argument("--preset", choices=["illumina", "ont-r10", "pacbio-hifi"], default="illumina", help="Sequencing technology preset")
    runall_parser.add_argument("--metric", choices=["wibs", "ensemble", "snp-wibs"], default=None, help="Distance metric to compute (defaults to wibs for de-novo/ont-r10, ensemble otherwise)")
    runall_parser.add_argument("-e", "--epsilon", type=float, default=0.3072, help="Bayesian error rate epsilon (default: 0.3072)")
    runall_parser.add_argument("--sra-accession", help="Optional SRA accession (e.g. SRR12345678) to fetch raw data directly")
    runall_parser.add_argument("--de-novo", action="store_true", help="Discover loci and haplotypes de novo without any reference database")
    runall_parser.add_argument("--ploidy", type=int, default=None, help="Organism base ploidy (e.g. 1 for haploid Cryptosporidium/bacteria, 2 for diploid)")
    runall_parser.add_argument("--min-completeness", type=float, default=0.10, help="Minimum locus completeness fraction (default: 0.10)")
    runall_parser.add_argument("--weight-mode", choices=["heterozygosity", "inverted-king", "king", "uniform"], default="heterozygosity", help="Allele weighting scheme for presence/absence indicator columns (default: heterozygosity)")
    runall_parser.add_argument("--min-maf", type=float, default=0.0, help="Minimum minor allele frequency threshold to filter rare/private singletons (default: 0.0)")
    runall_parser.add_argument("--no-psd", dest="project_psd", action="store_false", default=True, help="Skip Gram matrix PSD projection and return raw pairwise distances")
    runall_parser.add_argument("--k-min", type=int, default=2, help="Minimum number of clusters to search (default: 2)")
    runall_parser.add_argument("--k-max", type=int, default=50, help="Maximum number of clusters to search (default: 50)")
    runall_parser.add_argument("--relative-gap-floor", type=float, default=0.2200, help="Minimum relative merge-height gap fraction of tree height required for unsupervised knee selection (default: 0.2200)")
    runall_parser.add_argument(
        "--cut", choices=["count", "distance"], default="count",
        help="How to cut the tree. 'count' (default) chooses a number of clusters from the "
             "largest merge-height gap. 'distance' cuts at a fixed dissimilarity instead.")
    runall_parser.add_argument(
        "--linkage-method", choices=["ward", "single", "average", "complete"], default="ward",
        help="Linkage for the tree (default: ward).")
    runall_parser.add_argument(
        "--linkage-threshold", type=float, default=None,
        help="Dissimilarity at which to cut in --cut distance mode.")

    # The three amplicon subcommands own their own flags and are forwarded verbatim, so their
    # arguments must survive this parser rather than be rejected by it. Everything else is
    # parsed strictly, so a typo in an existing command still fails loudly.
    PASSTHROUGH = ("define-windows", "call-haplotypes", "build-sheet")
    if len(sys.argv) > 1 and sys.argv[1] in PASSTHROUGH:
        args, _ = parser.parse_known_args(sys.argv[1:2])
    else:
        args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    if args.command == "fetch-test-data":
        fetch_test_data(args.output_dir)

    elif args.command == "generate-sheet":
        input_target = args.assembled_fasta or args.specimen_dir
        if not input_target:
            print("[Error] Please specify either -s/--specimen-dir or -a/--assembled-fasta")
            sys.exit(1)
        generate_haplotype_sheet(
            specimen_dir=input_target,
            background_dir=args.background_dir,
            output_path=args.output,
            assembled_fasta=args.assembled_fasta,
            de_novo=args.de_novo
        )

    elif args.command == "process-ont":
        ref_db = {}
        if args.reference_fasta:
            if not os.path.exists(args.reference_fasta):
                print(f"[Error] Reference FASTA file not found: {args.reference_fasta}", file=sys.stderr)
                sys.exit(1)
            from .haplotype_sheet import parse_fasta
            ref_db = parse_fasta(args.reference_fasta)
        elif not args.de_novo:
            default_ref = os.path.join(os.path.dirname(__file__), "..", "references", "cyclopora_mlst_references.fasta")
            if os.path.exists(default_ref):
                from .haplotype_sheet import parse_fasta
                ref_db = parse_fasta(default_ref)
            else:
                print("[Error] process-ont requires a reference database. Please specify -r/--reference-fasta or --de-novo.", file=sys.stderr)
                sys.exit(1)

        processor = NanoporeAmpliconProcessor(min_qscore=args.qscore)
        processor.match_ont_haplotypes(args.sample_id, args.input_fastq, ref_db, args.output_dir)

    elif args.command == "eukaryotyping":
        import pandas as pd
        df = pd.read_csv(args.input_sheet, sep="\t")
        engine = PyEukDistanceEngine(
            epsilon=args.epsilon,
            min_completeness=args.min_completeness,
            ploidy=args.ploidy,
            weight_mode=args.weight_mode,
            min_maf=args.min_maf,
            project_psd=args.project_psd
        )
        if args.metric == "wibs" or args.wibs:
            res_df = engine.compute_revised_wibs_matrix(df)
        elif args.metric == "snp-wibs":
            res_df = engine.compute_snp_weighted_wibs_matrix(df, fasta_path=args.fasta)
        else:
            res_df = engine.compute_ensemble_matrix(df)
        out_path = args.output_matrix or "distance_matrix.csv"
        res_df.to_csv(out_path)
        print(f"[CLI] Saved distance matrix to: {out_path}")

    elif args.command == "cluster":
        import pandas as pd
        matrix_df = pd.read_csv(args.matrix, index_col=0)
        finder = CyclosporaClusterFinder(
            stringency=args.stringency,
            robust=args.robust,
            relative_gap_floor=args.relative_gap_floor
        )
        if args.single_k:
            # Legacy single-partition behaviour, opt-in.
            finder.find_clusters(
                matrix_df,
                args.gold_clusters,
                k_min=args.k_min,
                k_max=args.k_max,
                relative_gap_floor=args.relative_gap_floor,
                output_dir=args.output_dir,
                cut_mode=args.cut,
                linkage_threshold=args.linkage_threshold,
                linkage_method=args.linkage_method
            )
        else:
            # Default: the sweep diagnostic -- a count range, its confidence, and the
            # confidence tree. Also writes a representative partition for downstream tools.
            finder.cluster_sweep(
                matrix_df,
                k_min=args.k_min,
                k_max=args.k_max,
                n_boot=args.n_boot,
                linkage_method=args.linkage_method,
                output_dir=args.output_dir,
            )

    elif args.command in ("define-windows", "call-haplotypes", "build-sheet"):
        # argparse has already consumed the subcommand; hand the remainder to the module.
        from .amplicon import build_sheet, define_windows, window_haplotypes
        rest = sys.argv[2:]
        {"define-windows": define_windows.main,
         "call-haplotypes": window_haplotypes.main,
         "build-sheet": build_sheet.main}[args.command](rest)

    elif args.command == "run-all":
        os.makedirs(args.output_dir, exist_ok=True)
        sheet_path = os.path.join(args.output_dir, "haplotype_data_sheet.txt")
        matrix_path = os.path.join(args.output_dir, "ensemble_distance_matrix.csv")

        print("=== STAGE 1: Generating Haplotype Sheet ===")
        input_target = args.assembled_fasta or args.specimen_dir
        if not input_target:
            print("[Error] Please specify either -s/--specimen-dir or -a/--assembled-fasta for run-all")
            sys.exit(1)
        sheet_df = generate_haplotype_sheet(
            specimen_dir=input_target,
            background_dir=args.background_dir,
            output_path=sheet_path,
            assembled_fasta=args.assembled_fasta,
            de_novo=args.de_novo
        )
        all_specimens = sheet_df["Seq_ID"].tolist()

        print("\n=== STAGE 2: Running PyEuk Distance Engine ===")
        engine = PyEukDistanceEngine(
            epsilon=args.epsilon,
            min_completeness=args.min_completeness,
            ploidy=args.ploidy,
            weight_mode=args.weight_mode,
            min_maf=args.min_maf,
            project_psd=args.project_psd
        )
        if args.metric == "wibs" or (args.metric is None and (args.preset == "ont-r10" or args.de_novo)):
            print(f"[DistanceEngine] Using wIBS Distance Engine (weight_mode='{args.weight_mode}', min_maf={args.min_maf}, project_psd={args.project_psd})...")
            matrix_df = engine.compute_revised_wibs_matrix(sheet_df)
        elif args.metric == "snp-wibs":
            print(f"[DistanceEngine] Using SNP-Weighted wIBS Distance Engine (weight_mode='{args.weight_mode}', min_maf={args.min_maf}, project_psd={args.project_psd})...")
            matrix_df = engine.compute_snp_weighted_wibs_matrix(sheet_df, fasta_path=args.assembled_fasta)
        else:
            print(f"[DistanceEngine] Using Plucinski-Barratt Ensemble Distance Engine (project_psd={args.project_psd})...")
            matrix_df = engine.compute_ensemble_matrix(sheet_df)
        matrix_df.to_csv(matrix_path)

        print("\n=== STAGE 3: Outbreak Cluster Determination ===")
        finder = CyclosporaClusterFinder(
            relative_gap_floor=args.relative_gap_floor
        )
        finder.find_clusters(
            matrix_df,
            args.gold_clusters,
            k_min=args.k_min,
            k_max=args.k_max,
            relative_gap_floor=args.relative_gap_floor,
            output_dir=args.output_dir,
            all_input_ids=all_specimens,
            cut_mode=args.cut,
            linkage_threshold=args.linkage_threshold,
            linkage_method=args.linkage_method
        )

        print("\n==================================================")
        print("SUCCESS: Pipeline complete!")
        print(f"- Distance Matrix: {matrix_path}")
        print(f"- Outbreak Clusters: {args.output_dir}")
        print("==================================================")


if __name__ == "__main__":
    main()
