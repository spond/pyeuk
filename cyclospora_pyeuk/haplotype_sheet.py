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
from typing import List, Optional, Dict


def generate_haplotype_sheet(
    specimen_dir: str,
    background_dir: Optional[str] = None,
    output_path: Optional[str] = None
) -> pd.DataFrame:
    """
    Parses specimen blast genotype files (from directories or .zip archives) and produces a binary haplotype matrix.

    Parameters
    ----------
    specimen_dir : str
        Directory or .zip file containing tab-delimited genotype BLAST output files for newly processed specimens.
    background_dir : Optional[str]
        Directory or .zip file containing baseline/reference population genotype BLAST files.
    output_path : Optional[str]
        Path to save the generated TSV file. If None, saves to haplotype_sheets/YYYY-MM-DD_Cyclospora_haplotype_data_sheet.txt.

    Returns
    -------
    pd.DataFrame
        Genotype sheet matrix with Seq_ID as first column and presence marked by "X".
    """
    file_map: Dict[str, str] = {}

    def collect_source(path: str):
        if not path or not os.path.exists(path):
            return

        # If path is directly a .zip file
        if os.path.isfile(path) and path.endswith(".zip"):
            try:
                with zipfile.ZipFile(path, 'r') as zf:
                    for member in zf.namelist():
                        basename = os.path.basename(member)
                        if not member.endswith("/") and basename and not basename.startswith(".") and not member.startswith("__MACOSX"):
                            with zf.open(member) as f:
                                content = f.read().decode('utf-8', errors='ignore')
                                file_map[basename] = content
            except Exception as e:
                print(f"[HaplotypeSheet Warning] Could not read zip archive {path}: {e}")
            return

        # If path is a directory
        if os.path.isdir(path):
            for root, _, files in os.walk(path):
                for fname in files:
                    fpath = os.path.join(root, fname)
                    if fname.endswith(".zip"):
                        try:
                            with zipfile.ZipFile(fpath, 'r') as zf:
                                for member in zf.namelist():
                                    basename = os.path.basename(member)
                                    if not member.endswith("/") and basename and not basename.startswith(".") and not member.startswith("__MACOSX"):
                                        with zf.open(member) as f:
                                            content = f.read().decode('utf-8', errors='ignore')
                                            file_map[basename] = content
                        except Exception as e:
                            print(f"[HaplotypeSheet Warning] Could not read zip {fpath}: {e}")
                    elif not fname.startswith(".") and not fname.endswith(".pdf") and not fname.endswith(".sh"):
                        try:
                            with open(fpath, 'r', encoding='utf-8', errors='ignore') as f:
                                content = f.read()
                                file_map[fname] = content
                        except Exception as e:
                            pass

    collect_source(specimen_dir)
    if background_dir:
        collect_source(background_dir)

    if not file_map:
        raise ValueError(f"No genotype files or valid .zip archives found in: {specimen_dir}")

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
        output_path = os.path.join("haplotype_sheets", f"{today_str}_Cyclospora_haplotype_data_sheet.txt")

    genotype_sheet.to_csv(output_path, sep="\t", index=False)
    print(f"[HaplotypeSheet] Successfully generated haplotype sheet ({len(sorted_samples)} specimens, {len(sorted_markers)} markers) at: {output_path}")

    return genotype_sheet
