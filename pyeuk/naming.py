"""
Standardized Haplotype and Locus Naming Contract for PyEuk.

Provides consistent, round-trip bidirectional naming and parsing of homologous
locus windows and haplotype identifiers across CDC legacy, de novo, and amplicon pipelines.
"""

import re
import hashlib
from typing import Optional, Union


def parse_locus_name(col: str) -> str:
    r"""
    Extracts the homologous locus window name from a marker or haplotype identifier.

    Handles:
    - CDC format: Nu_378_PART_A_Hap_4 -> Nu_378_PART_A
    - De novo format: gp60_L752bp.H01_9180 -> gp60
    - Content-hashed format: 18S_L830bp.H01_32A1 -> 18S
    - Sub-haplotype & window formats: Locus_01_L150bp.H02 -> Locus_01
    - Haplotype delimiters: _Hap_\d+, .H\d+, _H\d+, _NOVEL_\d+, .X_\d+, _X_\d+
    - Preserves locus names with underscores, hyphens, numbers (e.g. Cp_HSP70, beta-tubulin_ex2)
    - Strips amplicon length tags (_L\d+bp)

    Invariants:
    For any valid locus string L and parameters P:
        parse_locus_name(name_haplotype(L, **P)) == L
    """
    sub = str(col).strip()
    if not sub:
        return ""

    # Special CDC mitochondrial junction case normalization
    if "Junction" in sub or ("Mt_" in sub and "Cmt" in sub):
        # Check if it has Hap or H suffix
        sub_no_junc = re.sub(r"_Hap_\d+|\.H\d+.*|_H\d+.*", "", sub)
        if sub_no_junc == "Mt_Cmt" or "Junction" in sub:
            return "Mt_Cmt"

    patterns = [
        r"_Hap_\d+",
        r"\.(H|h)\d+(_[A-Fa-f0-9]+)?",
        r"_(H|h)\d+(_[A-Fa-f0-9]+)?",
        r"_(NOVEL|novel)_\d+",
        r"\.(X|x)_\d+",
        r"_(X|x)_\d+",
    ]
    for pat in patterns:
        sub = re.sub(pat, "", sub)

    # Strip length suffix like _L245bp
    sub = re.sub(r"_L\d+bp$", "", sub)
    sub = sub.rstrip("_")
    return sub


def name_haplotype(
    locus: str,
    hap_id: Union[int, str] = 1,
    sequence: Optional[str] = None,
    length_bp: Optional[int] = None,
    seq_hash: Optional[str] = None,
    hash_len: int = 4,
    style: str = "de_novo",
    include_length: bool = True
) -> str:
    """
    Generates a canonical, globally reproducible haplotype identifier.

    Parameters
    ----------
    locus : str
        Name of the locus or amplicon partition (e.g. 'Nu_378_PART_A', 'gp60', 'ITS-2').
    hap_id : Union[int, str]
        Haplotype integer rank/index (e.g. 1, 2) or identifier string ('01', 'H01').
    sequence : Optional[str]
        Nucleotide sequence of the haplotype. If provided, used to derive length and hash.
    length_bp : Optional[int]
        Sequence length in base pairs. Derived from `sequence` if not explicitly given.
    seq_hash : Optional[str]
        Sequence hash (e.g. 4-hex string). Derived from `sequence` MD5 if not explicitly given.
    hash_len : int
        Number of characters for MD5 sequence hash (default: 4, 0 to disable).
    style : str
        Naming style: 'de_novo' (default), 'cdc', 'novel', or 'compact'.
    include_length : bool
        Whether to include '_L<Length>bp' tag in de novo format (default: True).

    Returns
    -------
    str
        Standardized haplotype identifier string.
    """
    clean_locus = str(locus).strip().replace(" ", "_")

    if sequence:
        seq_clean = sequence.upper().strip()
        if length_bp is None:
            length_bp = len(seq_clean)
        if seq_hash is None and hash_len > 0:
            seq_hash = hashlib.md5(seq_clean.encode("utf-8")).hexdigest()[:hash_len].upper()

    if isinstance(hap_id, int):
        hap_str = f"{hap_id:02d}" if hap_id < 100 else str(hap_id)
    else:
        hap_str = str(hap_id).lstrip("Hh_")
        if hap_str.isdigit() and len(hap_str) == 1:
            hap_str = f"0{hap_str}"

    if style == "cdc":
        return f"{clean_locus}_Hap_{int(hap_str) if hap_str.isdigit() else hap_str}"
    elif style == "novel":
        return f"{clean_locus}_NOVEL_{int(hap_str) if hap_str.isdigit() else hap_str}"
    elif style == "compact":
        return f"{clean_locus}.H{hap_str}"

    # Default: de_novo style (<Locus>[_L<Len>bp].H<Rank>[_<Hash>])
    prefix = clean_locus
    if include_length and length_bp is not None:
        prefix = f"{clean_locus}_L{length_bp}bp"

    if seq_hash and hash_len > 0:
        return f"{prefix}.H{hap_str}_{seq_hash}"
    else:
        return f"{prefix}.H{hap_str}"


def format_de_novo_haplotype_name(
    locus_name: str,
    sequence: str,
    rank: int = 1,
    include_length: bool = True,
    include_hash: bool = True
) -> str:
    """
    Backward-compatibility wrapper for `name_haplotype`.
    """
    return name_haplotype(
        locus=locus_name,
        hap_id=rank,
        sequence=sequence,
        hash_len=4 if include_hash else 0,
        include_length=include_length,
        style="de_novo"
    )
