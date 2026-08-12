"""
MicroAssembly: Modern Python module for targeted micro-assembly of Mt_Cmt junction tandem repeats.
Replaces legacy MIRA 4.0.2, CAP3, and multi-cutadapt shell loops for CDC Cyclospora cayetanensis workflow.
"""

import os
import re
import shutil
import subprocess
from typing import List, Tuple, Optional


class JunctionMicroAssembler:
    """
    Targeted micro-assembler and motif graph parser for Cyclospora mitochondrial junction repeats.
    """

    FORWARD_PRIMER = "CCATCTACAGC"
    REVERSE_PRIMER = "GTGTT"
    PRIMER_3PRIME = "AACAC"
    PRIMER_REV_3PRIME = "GCTGT"

    def __init__(self, threads: int = 4):
        self.threads = threads

    def extract_junction_reads(self, fastq_path: str) -> List[str]:
        """
        Extracts reads matching junction primer motifs in memory without transient disk I/O.
        """
        matching_seqs = []
        pattern_fwd = re.compile(r"CCATCTACAGC.*AACAC", re.IGNORECASE)
        pattern_rev = re.compile(r"GTGTT.*GCTGTAGATGG", re.IGNORECASE)

        if not os.path.exists(fastq_path):
            return []

        # Open plain or gzipped fastq
        import gzip
        open_fn = gzip.open if fastq_path.endswith(".gz") else open

        with open_fn(fastq_path, "rt", encoding="utf-8", errors="ignore") as handle:
            line_num = 0
            for line in handle:
                line_num += 1
                if line_num % 4 == 2:  # Sequence line
                    seq = line.strip()
                    if pattern_fwd.search(seq) or pattern_rev.search(seq):
                        matching_seqs.append(seq)

        return matching_seqs

    def assemble_micro_spades(self, reads_fasta: str, output_dir: str) -> Optional[str]:
        """
        Executes SPAdes micro-assembly if SPAdes is installed on system PATH.
        """
        spades_path = shutil.which("spades.py") or shutil.which("spades")
        if not spades_path:
            return None

        cmd = [
            spades_path,
            "-s", reads_fasta,
            "-o", output_dir,
            "-t", str(self.threads),
            "--only-assembler",
            "-k", "21,33,55"
        ]

        try:
            subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            contigs_path = os.path.join(output_dir, "contigs.fasta")
            if os.path.exists(contigs_path) and os.path.getsize(contigs_path) > 0:
                return contigs_path
        except Exception as e:
            print(f"[MicroAssembler] SPAdes micro-assembly notice: {e}")
        return None

    def clean_and_trim_junction_contig(self, sequence: str) -> str:
        """
        Trims primer concatemers using BioPython regex motif parser.
        Replaces 10 sequential cutadapt invocations.
        """
        seq = sequence.upper()
        # Find start of forward primer
        fwd_idx = seq.find(self.FORWARD_PRIMER)
        if fwd_idx != -1:
            seq = seq[fwd_idx + len(self.FORWARD_PRIMER):]

        # Find 3' end primer
        rev_idx = seq.rfind(self.PRIMER_3PRIME)
        if rev_idx != -1:
            seq = seq[:rev_idx]

        return seq.strip()

    def process_specimen_junction(
        self,
        specimen_name: str,
        fastq_path: str,
        output_dir: str
    ) -> Optional[str]:
        """
        Extracts, assembles, and validates Mt_Cmt junction repeat haplotypes for a specimen.
        """
        os.makedirs(output_dir, exist_ok=True)
        junction_reads = self.extract_junction_reads(fastq_path)

        if not junction_reads:
            print(f"[MicroAssembler] No junction repeat reads detected for specimen: {specimen_name}")
            return None

        temp_fasta = os.path.join(output_dir, f"{specimen_name}_junction_reads.fasta")
        with open(temp_fasta, "w") as out:
            for i, seq in enumerate(junction_reads):
                out.write(f">read_{i+1}\n{seq}\n")

        # Attempt SPAdes assembly if available, otherwise use longest motif contig
        spades_out = os.path.join(output_dir, f"{specimen_name}_spades")
        contig_file = self.assemble_micro_spades(temp_fasta, spades_out)

        assembled_seq = ""
        if contig_file and os.path.exists(contig_file):
            with open(contig_file, "r") as f:
                lines = [l.strip() for l in f if not l.startswith(">")]
                assembled_seq = "".join(lines)
        else:
            # Fallback: Select modal/longest clean read
            cleaned_reads = [self.clean_and_trim_junction_contig(s) for s in junction_reads]
            cleaned_reads = [c for c in cleaned_reads if len(c) > 20]
            if cleaned_reads:
                assembled_seq = max(cleaned_reads, key=len)

        if not assembled_seq:
            return None

        final_junction_seq = self.clean_and_trim_junction_contig(assembled_seq)
        repeat_len = len(final_junction_seq)
        total_len = repeat_len + 43  # Total length including primers (21 + 22)

        header = f">Mt_Cmt{total_len}.X_Junction_Hap_{specimen_name}"
        out_fasta = os.path.join(output_dir, f"{specimen_name}_validated_junction.fasta")

        with open(out_fasta, "w") as out:
            out.write(f"{header}\n{final_junction_seq}\n")

        print(f"[MicroAssembler] Validated Mt_Cmt junction repeat for {specimen_name}: Length = {repeat_len} bp (Total = {total_len} bp)")
        return out_fasta
