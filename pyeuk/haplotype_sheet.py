"""
Haplotype Sheet Generator module for CDC Cyclospora cayetanensis workflow.
Supports folders, raw text files, and compressed ZIP archives (.zip).
"""

import os
import glob
import zipfile
import io
import datetime
import pandas as pd
from typing import List, Optional, Dict, Tuple


def parse_fasta(fasta_path: str) -> Dict[str, str]:
    """
    Parses a FASTA file and returns a dictionary mapping sequence headers to uppercase nucleotide sequences.
    """
    records = {}
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        cur_header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur_header:
                    records[cur_header] = "".join(seq_lines).upper()
                cur_header = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if cur_header:
            records[cur_header] = "".join(seq_lines).upper()
    return records


def load_reference_haplotypes(ref_dirs: Optional[List[str]] = None) -> Dict[str, str]:
    """
    Loads reference haplotype sequences from reference directory or common paths.
    """
    import glob
    records = {}
    fpaths = []
    if ref_dirs:
        for rd in ref_dirs:
            if os.path.exists(rd):
                if os.path.isfile(rd):
                    fpaths.append(rd)
                else:
                    fpaths.extend(glob.glob(os.path.join(rd, "**/*.fasta"), recursive=True))
                    fpaths.extend(glob.glob(os.path.join(rd, "**/*.fa"), recursive=True))
    else:
        fpaths.extend(glob.glob("cdc_reference_data/**/*.fasta", recursive=True))
        fpaths.extend(glob.glob("cdc_reference_data/**/*.fa", recursive=True))

    for fp in fpaths:
        if "Illumina" in fp or "MAPPING" in fp:
            continue
        try:
            with open(fp, "r", encoding="utf-8", errors="ignore") as f:
                cur_header = None
                seq_lines = []
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        if cur_header:
                            records[cur_header] = "".join(seq_lines).upper()
                        cur_header = line[1:].strip().split()[0]
                        seq_lines = []
                    else:
                        seq_lines.append(line)
                if cur_header:
                    records[cur_header] = "".join(seq_lines).upper()
        except Exception:
            pass

    def norm_h(name):
        if "." in name:
            name = name.split(".")[-1]
        return name.replace("-", "_").strip()

    return {norm_h(k): seq for k, seq in records.items() if seq}


def match_assembled_contig(
    query_seq: str,
    ref_db: Dict[str, str],
    min_identity: float = 99.0,
    min_coverage: float = 85.0
) -> Optional[Tuple[str, float, float]]:
    """
    Matches an assembled sequence contig against reference MLST haplotype sequences.
    """
    if not query_seq or len(query_seq) < 20 or not ref_db:
        return None

    query = query_seq.upper().strip()
    qlen = len(query)

    best_match = None
    best_ident = 0.0
    best_cov = 0.0
    best_locus = None

    for hap_name, ref_seq in ref_db.items():
        rlen = len(ref_seq)
        if rlen == 0:
            continue

        # Exact substring or identity check
        if ref_seq in query or query in ref_seq:
            cov = min(qlen, rlen) / max(qlen, rlen) * 100.0
            if cov >= min_coverage:
                return (hap_name, 100.0, cov)

        # Semi-global alignment approximation
        min_len = min(qlen, rlen)
        max_len = max(qlen, rlen)
        cov = (min_len / max_len) * 100.0

        if cov < min_coverage:
            continue

        # Positional mismatch scan
        mismatches = sum(1 for a, b in zip(query[:min_len], ref_seq[:min_len]) if a != b)
        ident = ((min_len - mismatches) / min_len) * 100.0

        if ident > best_ident:
            best_ident = ident
            best_cov = cov
            best_match = hap_name
            # Extract locus prefix (e.g. Nu_378_PART_A)
            if "_Hap_" in hap_name:
                best_locus = hap_name.split("_Hap_")[0]
            elif "Junction" in hap_name:
                best_locus = "Mt_Cmt"
            else:
                best_locus = hap_name

    if best_match and best_ident >= min_identity and best_cov >= min_coverage:
        return (best_match, best_ident, best_cov)
    elif best_locus and best_cov >= min_coverage:
        # Novel haplotype allele for this locus window
        import hashlib
        short_hash = hashlib.md5(query.encode("utf-8")).hexdigest()[:6].upper()
        novel_id = f"{best_locus}_NOVEL_{short_hash}"
        return (novel_id, best_ident, best_cov)

    return None


def generate_haplotype_sheet_from_assemblies(
    assembled_input: str,
    reference_fasta_dir: Optional[str] = None,
    output_path: Optional[str] = None,
    min_identity: float = 99.0,
    min_coverage: float = 85.0
) -> pd.DataFrame:
    """
    Ingests externally assembled FASTA contigs (from directory, multi-FASTA, or ZIP),
    matches contigs against MLST reference haplotypes, and generates a binary presence/absence sheet.
    """
    ref_db = load_reference_haplotypes([reference_fasta_dir] if reference_fasta_dir else None)
    sample_contigs: Dict[str, List[str]] = {}

    def parse_fasta_content(content: str, default_sample: str):
        cur_sample = default_sample
        cur_seqs = []
        for line in content.splitlines():
            line = line.strip()
            if line.startswith(">"):
                if cur_seqs:
                    sample_contigs.setdefault(cur_sample, []).append("".join(cur_seqs))
                    cur_seqs = []
                header = line[1:].strip()
                if "|" in header:
                    cur_sample = header.split("|")[0].strip()
                elif "_" in header and not header.startswith("read"):
                    # e.g. SampleID_Contig1
                    parts = header.split("_")
                    if len(parts) > 1 and ("contig" in parts[-1].lower() or "hap" in parts[-1].lower()):
                        cur_sample = "_".join(parts[:-1])
                    else:
                        cur_sample = default_sample
                else:
                    cur_sample = default_sample
            else:
                cur_seqs.append(line)
        if cur_seqs:
            sample_contigs.setdefault(cur_sample, []).append("".join(cur_seqs))

    if os.path.isfile(assembled_input):
        if assembled_input.endswith(".zip"):
            with zipfile.ZipFile(assembled_input, "r") as zf:
                for member in zf.namelist():
                    basename = os.path.basename(member)
                    if not member.endswith("/") and basename and not basename.startswith(".") and not member.startswith("__MACOSX"):
                        sample_name = os.path.splitext(basename)[0]
                        with zf.open(member) as f:
                            text = f.read().decode("utf-8", errors="ignore")
                            parse_fasta_content(text, sample_name)
        else:
            sample_name = os.path.splitext(os.path.basename(assembled_input))[0]
            with open(assembled_input, "r", encoding="utf-8", errors="ignore") as f:
                parse_fasta_content(f.read(), sample_name)

    elif os.path.isdir(assembled_input):
        for root, _, files in os.walk(assembled_input):
            for fname in files:
                if fname.endswith((".fasta", ".fa", ".fna", ".fa.gz", ".fasta.gz")) and not fname.startswith("."):
                    sample_name = os.path.splitext(fname)[0]
                    fpath = os.path.join(root, fname)
                    import gzip
                    open_fn = gzip.open if fname.endswith(".gz") else open
                    try:
                        with open_fn(fpath, "rt", encoding="utf-8", errors="ignore") as f:
                            parse_fasta_content(f.read(), sample_name)
                    except Exception:
                        pass

    if not sample_contigs:
        raise ValueError(f"No valid assembled FASTA sequences found in: {assembled_input}")

    all_markers = set()
    sample_calls = {}

    for sample_id, contigs in sample_contigs.items():
        called = set()
        for cseq in contigs:
            match_res = match_assembled_contig(cseq, ref_db, min_identity=min_identity, min_coverage=min_coverage)
            if match_res:
                hap_id = match_res[0]
                called.add(hap_id)
                all_markers.add(hap_id)
            else:
                # If no reference match but valid sequence, use short hash for de novo locus
                import hashlib
                short_hash = hashlib.md5(cseq.upper().strip().encode("utf-8")).hexdigest()[:6].upper()
                novel_id = f"Novel_Marker_{short_hash}"
                called.add(novel_id)
                all_markers.add(novel_id)
        sample_calls[sample_id] = called

    sorted_markers = sorted(list(all_markers))
    sorted_samples = sorted(list(sample_calls.keys()))

    rows = []
    for seq_id in sorted_samples:
        present_markers = sample_calls[seq_id]
        row = {"Seq_ID": seq_id}
        for m in sorted_markers:
            row[m] = "X" if m in present_markers else ""
        rows.append(row)

    genotype_sheet = pd.DataFrame(rows, columns=["Seq_ID"] + sorted_markers)

    if output_path is None:
        today_str = datetime.date.today().strftime("%Y-%m-%d")
        os.makedirs("haplotype_sheets", exist_ok=True)
        output_path = os.path.join("haplotype_sheets", f"{today_str}_haplotype_data_sheet.txt")
    else:
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    genotype_sheet.to_csv(output_path, sep="\t", index=False)
    print(f"[HaplotypeSheet] Ingested external assemblies for {len(sorted_samples)} specimens ({len(sorted_markers)} markers called) -> {output_path}")
    return genotype_sheet


from .naming import (
    format_de_novo_haplotype_name,
    name_haplotype,
    parse_locus_name,
)


def learn_de_novo_haplotypes(
    assembled_input: str,
    output_path: Optional[str] = None,
    output_fasta: Optional[str] = None,
    min_identity: float = 1.0,
    locus_cluster_identity: float = 0.75,
    include_hash: bool = True
) -> Tuple[pd.DataFrame, Dict[str, str]]:
    """
    Reference-Free De Novo Haplotype Discovery & Typing (PyEuk De Novo Engine).
    Learns homologous locus windows and unique haplotypes directly from assembled FASTA contigs
    across specimens without requiring any pre-existing reference database.
    Applies the principled `<Locus>_L<Length>bp.H<Rank>_<Hash4>` naming scheme.

    Parameters
    ----------
    assembled_input : str
        Path to directory of FASTA files, multi-FASTA file, or .zip archive of assembled contigs.
    output_path : Optional[str]
        Path to save the generated binary presence/absence TSV sheet.
    output_fasta : Optional[str]
        Path to save the learned representative haplotype dictionary as a FASTA file.
    min_identity : float
        Identity threshold for grouping sequences into the same haplotype (default: 1.0 = exact 100%).
    locus_cluster_identity : float
        Sequence similarity threshold for grouping contigs into homologous locus windows (default: 0.75).
    include_hash : bool
        Whether to include 4-character MD5 content hash in haplotype identifiers.

    Returns
    -------
    Tuple[pd.DataFrame, Dict[str, str]]
        (Binary Haplotype Presence Matrix, Dictionary of {Haplotype_ID: Sequence})
    """
    raw_records: List[Tuple[str, Optional[str], str]] = []

    def parse_entry(content: str, default_sample: str):
        cur_sample = default_sample
        cur_locus = None
        cur_seqs = []
        for line in content.splitlines():
            line = line.strip()
            if line.startswith(">"):
                if cur_seqs:
                    s_clean = "".join(cur_seqs).upper()
                    if len(s_clean) >= 20:
                        raw_records.append((cur_sample, cur_locus, s_clean))
                    cur_seqs = []
                header = line[1:].strip()
                cur_locus = None
                if "|" in header:
                    parts = header.split("|")
                    cur_sample = parts[0].strip()
                    if len(parts) > 1:
                        cur_locus = parts[1].strip()
                elif "_" in header and not header.startswith("read"):
                    parts = header.split("_")
                    if len(parts) > 1 and ("contig" in parts[-1].lower() or "hap" in parts[-1].lower()):
                        cur_sample = "_".join(parts[:-1])
                    else:
                        cur_sample = default_sample
                else:
                    cur_sample = default_sample
            else:
                cur_seqs.append(line)
        if cur_seqs:
            s_clean = "".join(cur_seqs).upper()
            if len(s_clean) >= 20:
                raw_records.append((cur_sample, cur_locus, s_clean))

    if os.path.isfile(assembled_input):
        if assembled_input.endswith(".zip"):
            with zipfile.ZipFile(assembled_input, "r") as zf:
                for member in zf.namelist():
                    basename = os.path.basename(member)
                    if not member.endswith("/") and basename and not basename.startswith(".") and not member.startswith("__MACOSX"):
                        sample_name = os.path.splitext(basename)[0]
                        with zf.open(member) as f:
                            parse_entry(f.read().decode("utf-8", errors="ignore"), sample_name)
        else:
            sample_name = os.path.splitext(os.path.basename(assembled_input))[0]
            with open(assembled_input, "r", encoding="utf-8", errors="ignore") as f:
                parse_entry(f.read(), sample_name)

    elif os.path.isdir(assembled_input):
        for root, _, files in os.walk(assembled_input):
            for fname in files:
                if fname.endswith((".fasta", ".fa", ".fna", ".fa.gz", ".fasta.gz")) and not fname.startswith("."):
                    sample_name = os.path.splitext(fname)[0]
                    fpath = os.path.join(root, fname)
                    import gzip
                    open_fn = gzip.open if fname.endswith(".gz") else open
                    try:
                        with open_fn(fpath, "rt", encoding="utf-8", errors="ignore") as f:
                            parse_entry(f.read(), sample_name)
                    except Exception:
                        pass

    if not raw_records:
        raise ValueError(f"No valid assembled FASTA sequences found in: {assembled_input}")

    # Helper function for approximate alignment identity
    def approx_ident(s1: str, s2: str) -> float:
        min_l = min(len(s1), len(s2))
        max_l = max(len(s1), len(s2))
        if min_l / max_l < 0.60:
            return 0.0
        matches = sum(1 for a, b in zip(s1[:min_l], s2[:min_l]) if a == b)
        return float(matches) / float(max_l)

    # 1. Group sequences into homologous Loci
    locus_groups: Dict[str, List[Tuple[str, str]]] = {} # LocusName -> List of (sample_id, seq)
    locus_centroids: Dict[str, str] = {}
    locus_counter = 1

    for sample_id, locus_hint, seq in raw_records:
        assigned_locus = None
        if locus_hint:
            # Clean up locus hint
            l_name = locus_hint.split("_Hap_")[0].split(".")[0].strip()
            assigned_locus = l_name
            if assigned_locus not in locus_centroids:
                locus_centroids[assigned_locus] = seq
        else:
            # De novo clustering against existing locus centroids
            for loc_name, centroid_seq in locus_centroids.items():
                if approx_ident(seq, centroid_seq) >= locus_cluster_identity:
                    assigned_locus = loc_name
                    break

            if not assigned_locus:
                assigned_locus = f"Locus_{locus_counter:02d}"
                locus_counter += 1
                locus_centroids[assigned_locus] = seq

        locus_groups.setdefault(assigned_locus, []).append((sample_id, seq))

    # 2. Within each Locus, discover unique haplotypes and rank by cohort frequency
    learned_haplotype_dict: Dict[str, str] = {}
    sample_to_haplotypes: Dict[str, set] = {}
    all_all_samples = set(s for s, _, _ in raw_records)
    for s in all_all_samples:
        sample_to_haplotypes[s] = set()

    for loc_name, entries in sorted(locus_groups.items()):
        # Frequency of each unique sequence in this locus
        seq_counts: Dict[str, int] = {}
        for _, s in entries:
            seq_counts[s] = seq_counts.get(s, 0) + 1

        # Sort sequences by descending frequency, then length
        sorted_unique_seqs = sorted(seq_counts.keys(), key=lambda x: (-seq_counts[x], -len(x)))

        # Assign canonical, principled haplotype names (<Locus>_L<Len>bp.H<Rank>_<Hash4>)
        seq_to_hap_id = {}
        for h_idx, useq in enumerate(sorted_unique_seqs, start=1):
            hap_id = format_de_novo_haplotype_name(
                locus_name=loc_name,
                sequence=useq,
                rank=h_idx,
                include_length=True,
                include_hash=include_hash
            )
            seq_to_hap_id[useq] = hap_id
            learned_haplotype_dict[hap_id] = useq

        # Record calls for each specimen
        for sample_id, s in entries:
            hap_id = seq_to_hap_id[s]
            sample_to_haplotypes[sample_id].add(hap_id)

    sorted_markers = sorted(list(learned_haplotype_dict.keys()))
    sorted_samples = sorted(list(sample_to_haplotypes.keys()))

    rows = []
    for seq_id in sorted_samples:
        present_markers = sample_to_haplotypes[seq_id]
        row = {"Seq_ID": seq_id}
        for m in sorted_markers:
            row[m] = "X" if m in present_markers else ""
        rows.append(row)

    genotype_sheet = pd.DataFrame(rows, columns=["Seq_ID"] + sorted_markers)

    if output_path is None:
        today_str = datetime.date.today().strftime("%Y-%m-%d")
        os.makedirs("haplotype_sheets", exist_ok=True)
        output_path = os.path.join("haplotype_sheets", f"{today_str}_DeNovo_haplotype_data_sheet.txt")
    else:
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    genotype_sheet.to_csv(output_path, sep="\t", index=False)

    if output_fasta:
        os.makedirs(os.path.dirname(output_fasta) or ".", exist_ok=True)
        with open(output_fasta, "w") as out:
            for hid, hseq in sorted(learned_haplotype_dict.items()):
                out.write(f">{hid}\n{hseq}\n")
        print(f"[DeNovoLearner] Exported {len(learned_haplotype_dict)} learned haplotype reference sequences to: {output_fasta}")

    print(f"[DeNovoLearner] Reference-Free Discovery Complete: Learned {len(locus_groups)} loci and {len(sorted_markers)} unique haplotypes across {len(sorted_samples)} specimens -> {output_path}")
    return genotype_sheet, learned_haplotype_dict


def generate_haplotype_sheet(
    specimen_dir: str,
    background_dir: Optional[str] = None,
    output_path: Optional[str] = None,
    assembled_fasta: Optional[str] = None,
    de_novo: bool = False
) -> pd.DataFrame:
    """
    Parses specimen blast genotype files (or externally assembled FASTA contigs) and produces a binary haplotype matrix.
    If de_novo is True, discovers loci and haplotypes reference-free directly from sequence contigs.
    """
    if de_novo:
        target = assembled_fasta or specimen_dir
        output_fasta = None
        if output_path:
            out_dir = os.path.dirname(output_path) or "."
            output_fasta = os.path.join(out_dir, "learned_refs.fasta")
        sheet_df, _ = learn_de_novo_haplotypes(target, output_path=output_path, output_fasta=output_fasta)
        return sheet_df

    # If assembled_fasta is explicitly passed or specimen_dir is a FASTA file/folder of fastas
    if assembled_fasta:
        return generate_haplotype_sheet_from_assemblies(assembled_fasta, output_path=output_path)

    if specimen_dir.endswith((".fasta", ".fa", ".fna")):
        return generate_haplotype_sheet_from_assemblies(specimen_dir, output_path=output_path)

    specimen_map: Dict[str, str] = {}
    background_map: Dict[str, str] = {}
    is_fasta_dir = False

    def collect_source(path: str, target_dict: Dict[str, str], label: str):
        nonlocal is_fasta_dir
        if not path or not os.path.exists(path):
            return

        if os.path.isfile(path) and path.endswith(".zip"):
            try:
                with zipfile.ZipFile(path, 'r') as zf:
                    for member in zf.namelist():
                        basename = os.path.basename(member)
                        if not member.endswith("/") and basename and not basename.startswith(".") and not member.startswith("__MACOSX"):
                            if basename.endswith((".fasta", ".fa", ".fna")):
                                is_fasta_dir = True
                            with zf.open(member) as f:
                                content = f.read().decode('utf-8', errors='ignore')
                                sample_id = os.path.splitext(basename)[0]
                                if sample_id in target_dict:
                                    raise ValueError(f"Duplicate sample ID '{sample_id}' detected within {label} archive: {path}")
                                target_dict[sample_id] = content
            except Exception as e:
                if isinstance(e, ValueError):
                    raise
                print(f"[HaplotypeSheet Warning] Could not read zip archive {path}: {e}")
            return

        if os.path.isdir(path):
            for root, _, files in os.walk(path):
                for fname in files:
                    fpath = os.path.join(root, fname)
                    if fname.endswith((".fasta", ".fa", ".fna")):
                        is_fasta_dir = True
                    if fname.endswith(".zip"):
                        try:
                            with zipfile.ZipFile(fpath, 'r') as zf:
                                for member in zf.namelist():
                                    basename = os.path.basename(member)
                                    if not member.endswith("/") and basename and not basename.startswith(".") and not member.startswith("__MACOSX"):
                                        if basename.endswith((".fasta", ".fa", ".fna")):
                                            is_fasta_dir = True
                                        with zf.open(member) as f:
                                            content = f.read().decode('utf-8', errors='ignore')
                                            sample_id = os.path.splitext(basename)[0]
                                            if sample_id in target_dict:
                                                raise ValueError(f"Duplicate sample ID '{sample_id}' detected within {label} archive: {fpath}")
                                            target_dict[sample_id] = content
                        except Exception as e:
                            if isinstance(e, ValueError):
                                raise
                            print(f"[HaplotypeSheet Warning] Could not read zip {fpath}: {e}")
                    elif not fname.startswith(".") and not fname.endswith(".pdf") and not fname.endswith(".sh"):
                        try:
                            with open(fpath, 'r', encoding='utf-8', errors='ignore') as f:
                                content = f.read()
                                sample_id = os.path.splitext(fname)[0]
                                if sample_id in target_dict:
                                    raise ValueError(f"Duplicate sample ID '{sample_id}' detected within {label} directory: {fpath}")
                                target_dict[sample_id] = content
                        except Exception as e:
                            if isinstance(e, ValueError):
                                raise

    collect_source(specimen_dir, specimen_map, "specimen")
    if background_dir:
        collect_source(background_dir, background_map, "background")
        # Check collision between specimen and background
        collision = set(specimen_map.keys()) & set(background_map.keys())
        if collision:
            raise ValueError(f"Sample ID collision between specimen and background collections: {collision}")

    if not specimen_map:
        raise ValueError(f"No genotype files or valid .zip archives found in specimen directory: {specimen_dir}")

    # Primary case cohort is built from specimens (background controls are not emitted as patient case rows)
    file_map = specimen_map

    if is_fasta_dir:
        return generate_haplotype_sheet_from_assemblies(specimen_dir, output_path=output_path)

    all_markers = set()
    sample_data = {}

    for seq_id, text in file_map.items():
        text_clean = text.strip()
        if not text_clean:
            sample_data[seq_id] = set()
            continue

        try:
            df = pd.read_csv(io.StringIO(text_clean), sep="\t", header=None)
            if df.empty or 0 not in df.columns:
                sample_data[seq_id] = set()
                continue
            markers = set(df[0].dropna().astype(str).str.strip().tolist())
            sample_data[seq_id] = markers
            all_markers.update(markers)
        except Exception:
            sample_data[seq_id] = set()

    sorted_markers = sorted(list(all_markers))
    sorted_samples = sorted(list(sample_data.keys()))

    rows = []
    for seq_id in sorted_samples:
        present_markers = sample_data[seq_id]
        row = {"Seq_ID": seq_id}
        for m in sorted_markers:
            row[m] = "X" if m in present_markers else ""
        rows.append(row)

    genotype_sheet = pd.DataFrame(rows, columns=["Seq_ID"] + sorted_markers)

    if output_path is None:
        today_str = datetime.date.today().strftime("%Y-%m-%d")
        os.makedirs("haplotype_sheets", exist_ok=True)
        output_path = os.path.join("haplotype_sheets", f"{today_str}_haplotype_data_sheet.txt")

    genotype_sheet.to_csv(output_path, sep="\t", index=False)
    print(f"[HaplotypeSheet] Successfully generated haplotype sheet ({len(sorted_samples)} specimens, {len(sorted_markers)} markers) at: {output_path}")

    return genotype_sheet
