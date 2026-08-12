"""
Haplotype Sheet Generator module for CDC Cyclospora cayetanensis workflow.
Replaces legacy R script START_haplotype_sheet_generator.R.
"""

import os
import glob
import datetime
import pandas as pd
from typing import List, Optional


def generate_haplotype_sheet(
    specimen_dir: str,
    background_dir: Optional[str] = None,
    output_path: Optional[str] = None
) -> pd.DataFrame:
    """
    Parses specimen blast genotype files and produces a binary haplotype matrix dataframe.

    Parameters
    ----------
    specimen_dir : str
        Directory containing tab-delimited genotype BLAST output files for newly processed specimens.
    background_dir : Optional[str]
        Directory containing baseline/reference population genotype BLAST files (e.g. REFERENCE_POPULATION).
    output_path : Optional[str]
        Path to save the generated TSV file. If None, saves to haplotype_sheets/YYYY-MM-DD_Cyclospora_haplotype_data_sheet.txt.

    Returns
    -------
    pd.DataFrame
        Genotype sheet matrix with Seq_ID as first column and presence marked by "X".
    """
    file_map = {}

    def collect_files(directory: str):
        if not directory or not os.path.exists(directory):
            return
        for fpath in glob.glob(os.path.join(directory, "*")):
            if os.path.isfile(fpath) and not fpath.endswith(".zip") and not os.path.basename(fpath).startswith("."):
                seq_id = os.path.basename(fpath)
                file_map[seq_id] = fpath

    collect_files(specimen_dir)
    if background_dir:
        collect_files(background_dir)

    if not file_map:
        raise ValueError("No genotype files found in the specified directories.")

    all_markers = set()
    sample_data = {}

    for seq_id, fpath in file_map.items():
        if os.path.getsize(fpath) == 0:
            sample_data[seq_id] = set()
            continue

        try:
            df = pd.read_csv(fpath, sep="\t", header=None)
            if df.empty or 0 not in df.columns:
                sample_data[seq_id] = set()
                continue
            markers = set(df[0].dropna().astype(str).str.strip().tolist())
            sample_data[seq_id] = markers
            all_markers.update(markers)
        except Exception as e:
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
