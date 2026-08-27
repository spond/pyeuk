#!/usr/bin/env python3
"""Call haplotypes as read-level window strings, named by their own content.

Nothing here is species- or panel-specific: it takes a BAM and the reference the BAM was
aligned against, and works for any amplicon set. Windows are derived from the data rather
than from a curated BED.

Why windows at all, rather than whole amplicons: in a fragmented (Tn5/Nextera) library the
reads are pieces of the amplicon, not the amplicon. Linkage is directly observable only
across the span of one fragment. A window sized to the fragment distribution is therefore
the largest unit whose haplotype can be *observed* instead of inferred.

Why content-derived names: a haplotype is named by how its string differs from the
reference over that window (HGVS-like, 1-based within the window), so no curated
nomenclature file is needed and the name can be compared to another name. A read string
identical to the reference is named "=" -- an ordinary haplotype, which is what keeps
"amplified and reference-identical" distinguishable from "not amplified".

Why frequency is kept: it is the fraction of spanning reads carrying that string, i.e. a
property of one haplotype, not a sum over haplotypes. Thresholding per-site allele
frequency cannot separate a minor strain from a shared site; thresholding haplotype
frequency can.

Output is long-format, one row per called haplotype, so several haplotypes per window per
specimen (a mixture) are represented natively.

Usage:
  window_haplotypes.py --bam S.bam --ref panel.fa --specimen S [--bed w.bed]
                       [--window auto|INT] [--step INT] [--min-span 50]
                       [--min-reads 10] [--min-freq 0.05] [--out calls.tsv]
"""
import argparse
import os
import sys
from collections import Counter, defaultdict

# pysam is an OPTIONAL dependency, declared under the `amplicon` extra. Importing it at module
# scope would make `import pyeuk` fail for anyone who installed the core package to
# work on sheets from SeekDeep or DADA2 and has no interest in BAMs. Resolved on first use
# instead, with a message that says what to install.
_pysam = None


def _require_pysam():
    global _pysam
    if _pysam is None:
        try:
            import pysam as _p
        except ImportError as exc:
            raise SystemExit(
                "This step reads BAM files and needs pysam, which is an optional dependency.\n"
                "  pip:   pip install 'pyeuk[amplicon]'\n"
                "  conda: conda install -c bioconda pysam\n"
                "Conda has no concept of extras, so the bioconda `pyeuk` package carries pysam "
                "as an ordinary run dependency and this message should not appear there. If it "
                "does, the environment is not the one pyeuk was installed into."
            ) from exc
        _pysam = _p
    return _pysam


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def aligned_pairs(pos, cigar, seq):
    """Map reference position -> read base. Insertions are appended lowercase to the
    preceding reference position; deletions are '-'. Soft/hard clips consume no reference."""
    ref, read, aln = pos, 0, {}
    n = ""
    for ch in cigar:
        if ch.isdigit():
            n += ch
            continue
        c = int(n or 0)
        n = ""
        if ch in "M=X":
            for i in range(c):
                aln[ref + i] = seq[read + i]
            ref += c
            read += c
        elif ch == "I":
            if ref - 1 in aln:
                aln[ref - 1] += seq[read:read + c].lower()
            read += c
        elif ch == "D":
            for i in range(c):
                aln[ref + i] = "-"
            ref += c
        elif ch == "S":
            read += c
        elif ch == "N":
            ref += c
        # H consumes neither
    return aln


def denoise(obs, max_edits, ratio):
    """Fold rare strings into an abundant neighbour, UNOISE/DADA2-style.

    A window string is a whole molecule, so a single sequencing error creates a NEW
    haplotype rather than shifting one base of an existing one. At 250 bp most reads
    carry at least one error and a clean single-species control shatters: measured at
    62% of reads on its own exact string, against 98% at 90 bp. The residue is a cloud
    of low-abundance neighbours around each true haplotype.

    Two consequences, both bad. Quantitatively the true haplotype is under-counted, and
    unevenly so -- a haplotype containing a homopolymer indel loses far more of its
    reads than a clean one. Structurally every specimen acquires private haplotypes,
    which are rare by construction, and PyEuk weights a column by 1/sqrt(p(1-p)), which
    grows without bound as p approaches 0. Rare private strings therefore get the
    largest weights in the matrix.

    A string is folded into a more abundant one when it is within max_edits positions
    of it AND that neighbour is at least `ratio` times more abundant, which is the
    condition that distinguishes "error off a real haplotype" from "genuine minor
    variant". Folding is done most-abundant-first so chains collapse to their peak.
    """
    if max_edits <= 0 or len(obs) < 2:
        return obs
    order = [h for h, _ in obs.most_common()]
    keep, merged = [], {}
    for h in order:
        target = None
        for k in keep:
            if obs[k] < obs[h] * ratio:
                continue
            a, b = h.split("\t"), k.split("\t")
            if len(a) != len(b):
                continue
            d = 0
            for x, y in zip(a, b):
                if x != y:
                    d += 1
                    if d > max_edits:
                        break
            if d <= max_edits:
                target = k
                break
        if target is None:
            keep.append(h)
        else:
            merged[h] = target
    out = Counter()
    for h, c in obs.items():
        out[merged.get(h, h)] += c
    return out


def left_normalize(obs, ref):
    """Shift indels to their leftmost equivalent position, as variant callers do.

    `obs` is one entry per reference position in the window: entry[0] is the aligned base or
    '-' for a deletion, and any inserted bases follow in lowercase. The aligner is free to
    place an indel anywhere within a repeat, so the SAME molecule can be represented several
    ways and therefore acquire several haplotype names. Measured on Teladorsagia beta-tubulin,
    reads carrying the same call were aligned four different ways (300M, 172M8D128M,
    97M1D40M1D19M4I2M4D138M, 97M1D40M1D163M), and only 19.8% aligned without an indel at all.

    A deletion may move one position left when the deleted base equals the base before it;
    an insertion may move left when its last base equals the reference base it follows.
    Repeating until neither applies gives a canonical placement independent of the aligner.

    This does not merge genuinely different sequences -- it only removes representational
    freedom. Two reads with different bases stay different haplotypes."""
    if not obs:
        return obs
    o = list(obs)
    # deletions: '-' at i can swap with a non-deleted i-1 when ref[i] == ref[i-1]
    moved = True
    while moved:
        moved = False
        for i in range(1, len(o)):
            if o[i][0] == "-" and o[i - 1][0] != "-" and not o[i - 1][1:] \
                    and ref[i].upper() == ref[i - 1].upper():
                o[i - 1], o[i] = "-" + o[i - 1][1:], o[i - 1][0] + o[i][1:]
                moved = True
    # insertions: bases appended after i move to i-1 when the last inserted base == ref[i]
    moved = True
    while moved:
        moved = False
        for i in range(1, len(o)):
            ins = o[i][1:]
            if ins and o[i][0] != "-" and ins[-1].upper() == ref[i].upper() \
                    and not o[i - 1][1:] and o[i - 1][0] != "-":
                o[i] = o[i][0]
                o[i - 1] = o[i - 1][0] + (ref[i].lower() if False else ins[-1]) + ins[:-1]
                moved = True
    return o


def name_haplotype(obs, ref, first_pos=1):
    """HGVS-like description of obs relative to ref.

    Both strings are per-reference-position: obs[i] is what the read carried at the i-th
    reference base, '-' for a deletion, with any inserted bases appended lowercase.

    `first_pos` is the coordinate given to the window's first base:

      reference (default)  first_pos = the window's start on the contig, so a variant is
                           named by where it sits in the PANEL. The name is then a property
                           of the molecule, not of the tiling.
      window               first_pos = 1, naming positions within the window.

    Why reference is the default. Windows are derived from the reads, so they move between
    cohorts: across five disjoint 20-BAM subsets of one cohort, only 1 of 70 window names was
    shared by all five, because one base of difference in where coverage starts renames every
    window on that contig and shifts every coordinate inside it. Under window-relative naming
    the same substitution is therefore called 45T>A in one run and 46T>A in the next, and two
    runs' haplotype names cannot be compared at all. Under reference naming it is 350T>A in
    both, so names remain comparable even when the tiling does not match.

    This does NOT make two runs' haplotypes identical: a window covering 6-105 and one
    covering 5-104 span different sequence, so a variant at 105 appears in one and not the
    other. It makes the names DIFFER MEANINGFULLY (different content) instead of spuriously
    (same content, different coordinate origin)."""
    diffs = []
    for i, (o, r) in enumerate(zip(obs, ref), start=first_pos):
        base, ins = o[0], o[1:]
        if base == "-":
            diffs.append(f"{i}del{r}")
        elif base.upper() != r.upper():
            diffs.append(f"{i}{r}>{base.upper()}")
        if ins:
            diffs.append(f"{i}_{i+1}ins{ins.upper()}")
    return ",".join(diffs) if diffs else "="


def block_lengths(bam):
    with _require_pysam().AlignmentFile(bam, "rb") as fh:
        return sorted(r.reference_length or 0 for r in fh.fetch()
                      if not (r.is_unmapped or r.is_secondary or r.is_supplementary))


def covered_intervals(bam, min_depth):
    """Per contig, the first and last reference position covered at >= min_depth.

    Data-driven so no primer coordinates or panel BED are needed: whatever the reads
    reach is what gets tiled, for any panel in any species."""
    lo, hi = {}, {}
    with _require_pysam().AlignmentFile(bam, "rb") as fh:
        for contig in fh.references:
            cov = fh.count_coverage(contig, quality_threshold=0)
            depth = [a + c + g + t for a, c, g, t in zip(*cov)]
            hits = [i + 1 for i, d in enumerate(depth) if d >= min_depth]
            if hits:
                lo[contig], hi[contig] = hits[0], hits[-1]
    return {c: (lo[c], hi[c]) for c in lo}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam", required=True)
    p.add_argument("--ref", required=True, help="FASTA the BAM was aligned to")
    p.add_argument("--specimen", required=True)
    p.add_argument("--bed", help="explicit windows; otherwise the reference is tiled")
    p.add_argument("--window", default="auto",
                   help="window size, or 'auto' to size it from the fragment distribution")
    p.add_argument("--window-pct", type=float, default=25.0,
                   help="auto sizing: percentile of aligned block length to use (default 25, "
                        "so ~75%% of reads can span a window)")
    p.add_argument("--window-min", type=int, default=40)
    p.add_argument("--window-max", type=int, default=250)
    p.add_argument("--step", type=int, default=0, help="0 = non-overlapping tiling")
    p.add_argument("--min-span", type=int, default=50,
                   help="minimum spanning reads for the window to be called at all")
    p.add_argument("--min-reads", type=int, default=10, help="per haplotype")
    p.add_argument("--min-freq", type=float, default=0.05, help="per haplotype")
    p.add_argument("--normalize", choices=["left", "none"], default="left",
                   help="left-align indels before naming, so the aligner's freedom to place a "
                        "gap inside a repeat cannot split one molecule across several "
                        "haplotype names. 'none' reproduces the previous behaviour.")
    p.add_argument("--coord-system", choices=["reference", "window"], default="reference",
                   help="coordinate origin for haplotype names. 'reference' (default) names a "
                        "variant by its position on the panel contig, so the name survives a "
                        "change of windowing; 'window' numbers from 1 within each window, which "
                        "is the historic behaviour and is only comparable within one run.")
    p.add_argument("--denoise-edits", type=int, default=1,
                   help="fold a string into an abundant neighbour within this many differing "
                        "positions; 0 disables denoising")
    p.add_argument("--denoise-ratio", type=float, default=8.0,
                   help="the neighbour must be at least this many times more abundant")
    p.add_argument("--out", default="-")
    a = p.parse_args(argv)

    ref = read_fasta(a.ref)

    # Auto window sizing is only needed when no BED is supplied. Running it regardless
    # meant a specimen whose BAM held no aligned reads exited non-zero even though the
    # windows were already fully determined by the BED -- and in Galaxy one such job
    # pauses the downstream sheet builder and strands the whole cohort.
    if a.bed:
        w = a.window_min if a.window == "auto" else int(a.window)
    elif a.window == "auto":
        L = block_lengths(a.bam)
        if not L:
            sys.exit(f"{a.bam}: no aligned reads and no --bed, so windows cannot be defined")
        w = L[max(0, int(len(L) * a.window_pct / 100.0) - 1)]
        w = max(a.window_min, min(a.window_max, int(w)))
    else:
        w = int(a.window)
    step = a.step or w

    windows = []
    if a.bed:
        for line in open(a.bed):
            if not line.strip() or line.startswith(("#", "track")):
                continue
            f = line.split()
            windows.append((f[0], int(f[1]) + 1, int(f[2])))   # BED is 0-based half-open
    else:
        # Tile only the region reads actually reach. Amplicon reads start inside the
        # primer, so tiling from position 1 spends the first window on sequence no read
        # spans -- which silently drops whatever variation sits just past the primer.
        covered = covered_intervals(a.bam, a.min_span)
        for c, s in ref.items():
            lo, hi = covered.get(c, (None, None))
            if lo is None:
                continue
            for st in range(lo, hi + 1, step):
                en = min(st + w - 1, hi)
                if en - st + 1 >= a.window_min:
                    windows.append((c, st, en))

    by_contig = defaultdict(list)
    for c, s, e in windows:
        by_contig[c].append((s, e))

    # A specimen can legitimately have zero surviving reads -- a failed library, or
    # everything removed by the MAPQ/proper-pair filter. That is a NOT_CALLED result for
    # every window, not an error: the sheet builder needs a row saying "this specimen was
    # examined and nothing was called", and the cohort must not be held hostage to one
    # dud sample.
    rows = []
    for contig, wins in by_contig.items():
        if contig not in ref:
            continue
        # Stream per window rather than materialising the whole contig.
        #
        # Holding every read of a contig as a position->base dict costs
        # O(reads x read_length) Python objects. On a 20,000x amplicon that ran to
        # 815 MB for one specimen and got a Galaxy job killed outright, with exit code
        # None and empty stderr. Fetching per window touches only the reads overlapping
        # it, keeps one dict alive at a time, and is cheap because the fetch is indexed.
        with _require_pysam().AlignmentFile(a.bam, "rb") as fh:
            if contig not in fh.references:
                continue
            for (s, e) in wins:
                obs = Counter()
                refstr = ref[contig][s - 1:e]
                for r in fh.fetch(contig, s - 1, e):
                    # Secondary and supplementary alignments are the SAME molecule
                    # placed again. Counting them as separate spanning reads votes one
                    # molecule twice, inflates the depth a gate is judged against, and
                    # can promote a chimeric split alignment into its own haplotype.
                    if (r.is_unmapped or r.is_secondary or r.is_supplementary
                            or not r.cigarstring or not r.query_sequence):
                        continue
                    # reject non-spanning reads before building the dict at all
                    if r.reference_start + 1 > s or (r.reference_end or 0) < e:
                        continue
                    aln = aligned_pairs(r.reference_start + 1, r.cigarstring,
                                        r.query_sequence)
                    if all(pp in aln for pp in range(s, e + 1)):
                        parts = [aln[pp] for pp in range(s, e + 1)]
                        # normalise BEFORE counting: two reads whose indel the aligner placed
                        # differently are the same molecule and must land in the same bin.
                        # Normalising after counting only renames survivors and merges nothing.
                        if a.normalize == "left":
                            parts = left_normalize(parts, refstr)
                        obs["\t".join(parts)] += 1
                obs = denoise(obs, a.denoise_edits, a.denoise_ratio)
                n = sum(obs.values())
                wname = f"{contig}_W{s:04d}"
                if n < a.min_span:
                    rows.append([a.specimen, contig, wname, s, e, "NOT_CALLED", 0, 0.0, n])
                    continue
                emitted = 0
                for st, c in obs.most_common():
                    if c < a.min_reads or c / n < a.min_freq:
                        continue
                    first = s if a.coord_system == "reference" else 1
                    rows.append([a.specimen, contig, wname, s, e,
                                 name_haplotype(st.split("\t"), refstr, first), c,
                                 round(c / n, 5), n])
                    emitted += 1
                # A window can clear min_span and still have every haplotype fail
                # min_reads/min_freq -- an error-shattered window, a low-frequency mixture,
                # or simply aggressive gates. Emitting nothing there is not the same as
                # saying nothing happened: the sheet builder discovers specimens FROM these
                # rows, so a specimen whose every window is filtered this way disappears
                # from the cohort silently rather than being reported as excluded. Emit the
                # sentinel. The two NOT_CALLED cases stay distinguishable downstream because
                # the spanning column is always the true count: spanning < min_span means
                # the window was never called, spanning >= min_span means it was called and
                # nothing survived the gates.
                if emitted == 0:
                    rows.append([a.specimen, contig, wname, s, e, "NOT_CALLED", 0, 0.0, n])

    fh = sys.stdout if a.out == "-" else open(a.out, "w")
    print("specimen\tlocus\twindow\tstart\tend\thaplotype\treads\tfreq\tspanning", file=fh)
    for r in rows:
        print("\t".join(str(x) for x in r), file=fh)
    if fh is not sys.stdout:
        fh.close()
    called = len({r[2] for r in rows if r[5] != "NOT_CALLED"})
    print(f"[{a.specimen}] window={w} step={step}  windows={len(windows)}  "
          f"called={called}  rows={len(rows)}", file=sys.stderr)


if __name__ == "__main__":
    main()
