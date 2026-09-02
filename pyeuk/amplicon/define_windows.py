#!/usr/bin/env python3
"""Derive a window BED for a cohort, once, from the reads themselves.

Windows must be identical across every specimen or the resulting sheets have
non-corresponding columns, so window choice is a cohort-level decision, not a
per-specimen one. This scans a sample of BAMs and writes a BED that
window_haplotypes.py then applies unchanged to all of them.

Three things are measured, all panel- and species-agnostic:

  extent      the reference positions reads actually reach, at --min-depth. Amplicon reads
              start inside the primer, so tiling from position 1 wastes the first window on
              sequence nothing spans.

  placement   the SPANNABLE CORE inside that extent -- the widest interval that --min-spanning
              of its overlapping reads cross end to end. Measuring width alone is not enough:
              amplicon reads arrive as two strand populations whose starts differ, and tiling
              from the first covered base lets only the earlier one span, capping the reachable
              width before any width search runs. On Teladorsagia beta-tubulin (300 bp reads,
              amplicon 195-511, forward starts ~195 and reverse ~215) tiling from 195 gave a
              150 bp window at 32.9% spanning; the core 215-490 gives 276 bp at 91.3%, nearly
              doubling both width and depth, and 276 bp is the whole amplicon.

  size        the widest candidate width whose MEASURED spanning fraction inside that core
              clears --min-spanning. Measured, not inferred from a read-length percentile: a
              read spans a window only if it is long enough AND starts in the right place, so
              among overlapping reads the spanning share is roughly (L-W+1)/(L+W-1). The old
              --spanning-target 0.7 delivered 30-33% while claiming 70%. It survives only as
              an explicit opt-in for reproducing a pre-2026-08 run.

Secondary and supplementary alignments are excluded everywhere. They are the same molecule
placed again, and counting them votes one molecule more than once.

PERFORMANCE (issue #20). On a 66-library cohort this step was the single-threaded serial
bottleneck (~150 min). Three things made it slow, all fixed here without changing the chosen
windows (see the two correctness gates in the tests):

  1. spanning_profile() re-fetched every read from every BAM once per candidate width (~34
     widths). Both the placement search and the width search are now derived from ONE
     per-contig interval histogram (a Counter of (start,end)->count) built by read_intervals()
     -- a single pass over each contig's reads per BAM instead of ~34 re-fetches. See
     spanning_profile_hist(): summing histogram counts per window is provably equal to
     counting fetched reads per window, because the histogram keys ARE the (start,end) pairs
     fetch() would return and the same non-overlapping tiling is walked.
  2. block_lengths() re-read every read a second time only to take a per-contig read-length
     median; reference_length == end - (start-1), so that median is read straight off the same
     histogram and the extra full pass is gone (block_lengths() is kept for API compatibility).
  3. covered() ran count_coverage() over all 738 header contigs of every BAM even though only
     the read-bearing ones (23 here) can ever clear --min-depth. It now consults the BAM index
     (get_index_statistics) and scans only contigs that carry alignments -- a contig with zero
     mapped reads contributes zero coverage, so the covered bounds are byte-for-byte unchanged.

  Per-contig window selection is independent across reference contigs, so it is parallelised
  with a thread pool (--threads, honouring $GALAXY_SLOTS); pysam releases the GIL during BAM
  I/O. Results are keyed by contig and emitted in sorted order, so the BED does not depend on
  thread count or completion order.

  --max-reads-per-window caps how many reads are folded into a contig's histogram (spanning is
  a FRACTION, so a bounded sample estimates it), seeded by --seed for reproducibility. Default
  0 means no cap, i.e. behaviour identical to before this flag existed. Capping is NOT
  full-depth equivalent: it perturbs the set of read-start positions, so window PLACEMENT can
  shift (a ~1bp anchor move on deep contigs cascades a phase offset through the whole tiling)
  and borderline widths -- and therefore window counts -- can flip across --min-spanning.
  Seed-stability only sets in at a high cap (~100k reads on the reference cohort), and even
  there the placement anchor still differs from full depth. Treat capped runs as an estimate.

Usage:
  define_windows.py --ref panel.fa --bams a.bam b.bam ... [--sample 20]
                    [--spanning-target 0.7] [--min-depth 20] [--out windows.bed]
                    [--threads N] [--max-reads-per-window K --seed S]
"""
import argparse
import os
import random
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor

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


def block_lengths(bam):
    """Aligned block length per read, keyed by contig.

    Per contig, not pooled: amplicons differ in depth and in fragment length, so one
    cohort-wide window size is set by whichever amplicon is deepest and is then far too
    wide for the shallow ones, which lose their locus entirely. CDC's own PART windows
    differ in size per marker for the same reason.

    KEPT FOR API COMPATIBILITY. main() no longer calls this: reference_length equals
    end - (start - 1), so the per-contig read-length median is read straight off the
    interval histogram read_intervals() already builds (see _nth_length), sparing this
    second full pass over every read. The filter here is identical to read_intervals()'s,
    so the lengths recovered from the histogram are byte-for-byte the same multiset."""
    lens = defaultdict(list)
    with _require_pysam().AlignmentFile(bam, "rb") as fh:
        for r in fh.fetch():
            if (r.is_unmapped or r.is_secondary or r.is_supplementary
                    or not r.reference_length):
                continue
            lens[r.reference_name].append(r.reference_length)
    return lens


def read_intervals(bams, contig):
    """Counter over (start, end) reference intervals of primary alignments.

    A Counter rather than a list because amplicon reads take very few distinct intervals --
    they start at a primer -- so the placement search below runs over a few hundred keys
    instead of a few million reads.

    THIS IS THE SINGLE PASS. It is built once per contig and then feeds BOTH the placement
    search (spannable_core) AND the width search (spanning_profile_hist); neither re-reads the
    BAM. Keys are (reference_start + 1, reference_end), the exact 1-based-inclusive interval
    fetch() reports, so a window's overlapping/spanning read counts are recoverable from the
    histogram alone."""
    iv = Counter()
    for bam in bams:
        with _require_pysam().AlignmentFile(bam, "rb") as fh:
            if contig not in fh.references:
                continue
            for r in fh.fetch(contig):
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                if not r.reference_length:
                    continue
                iv[(r.reference_start + 1, r.reference_end)] += 1
    return iv


def subsample_counter(iv, cap, seed):
    """Deterministically cap the total reads in an interval histogram to `cap`.

    Spanning is a FRACTION, so a bounded, representative sample of a contig's reads estimates
    it without touching every read on a 3M-deep amplicon. `cap` is the ceiling on the total
    number of reads (summed counts) folded into the contig's histogram, not a per-key cap:
    reads sharing a (start,end) interval are interchangeable for both searches, so we sample
    `cap` reads uniformly at random from the whole contig and rebuild the histogram from the
    survivors.

    Reproducible given `seed`: the contig's reads are laid out in a fixed order (keys sorted),
    then `cap` global indices are drawn by a seeded RNG. Same seed -> same survivors -> same
    windows. Different seeds are what the subsampling-safety gate compares."""
    if cap <= 0:
        return iv
    total = sum(iv.values())
    if total <= cap:
        return iv
    items = sorted(iv.items())               # fixed layout so the draw is reproducible
    chosen = random.Random(seed).sample(range(total), cap)
    chosen.sort()
    out = Counter()
    ci = 0
    base = 0
    for key, c in items:
        end = base + c
        k = 0
        while ci < len(chosen) and chosen[ci] < end:
            k += 1
            ci += 1
        if k:
            out[key] = k
        base = end
    return out


def _nth_length(iv, k):
    """The k-th smallest aligned read length (0-based) implied by the interval histogram.

    reference_length == reference_end - reference_start == end - (start - 1) == end - start + 1,
    so the multiset of read lengths block_lengths() would collect for a contig is recoverable
    from its (start,end) histogram. Used for the informational per-contig median and for the
    deprecated --spanning-target percentile, both of which block_lengths() used to serve."""
    lengths = {}
    for (s, e), c in iv.items():
        L = e - s + 1
        lengths[L] = lengths.get(L, 0) + c
    if not lengths:
        return 0
    acc = 0
    for L in sorted(lengths):
        acc += lengths[L]
        if k < acc:
            return L
    return max(lengths)


def spannable_core(iv, lo, hi, floor, step):
    """Widest interval within [lo, hi] that `floor` of its overlapping reads span end to end.

    THE REASON THIS EXISTS. Choosing a width is not enough; where the tiling STARTS decides
    what widths are even reachable. Amplicon reads arrive as two strand populations whose
    starts differ, and tiling from the first covered base lets only one of them span.

    Measured on Teladorsagia beta-tubulin, 300 bp reads, amplicon 195-511, forward reads
    starting ~195 and reverse ~215:

        tiling from 195   best width 150   spanning 32.9%   25,423 spanning reads
        core 215-490      width      276   spanning 91.3%   49,656 spanning reads

    Moving the origin 20 bp nearly doubles both the width and the depth, and 276 bp is the
    whole amplicon -- the same unit the source publication analysed, with all three
    resistance codons in one haplotype string instead of split across two windows.

    SCORED BY EXPECTED SPANNED BASES, width x spanning fraction, not by width alone.
    Maximising width subject to a floor picks a barely-passing wide interval over a slightly
    narrower near-perfect one, and that is the wrong trade. Measured on the same amplicon:

        213-510   298 bp x 0.3995 spanning  =  119 expected spanned bases
        215-491   277 bp x 0.9712 spanning  =  269 expected spanned bases

    The second is 7% narrower and carries 2.4x the spanning reads, and it is also the exact
    277 bp frame the source publication analysed. Width alone chose the first.

    Searched over a coarse grid of starts and ends so the cost stays linear in the number of
    distinct read intervals. Consumes the interval histogram directly (no BAM I/O)."""
    pairs = list(iv.items())
    if not pairs:
        return lo, hi
    starts = sorted({s for (s, _), _ in pairs if lo <= s <= hi})
    ends = sorted({e for (_, e), _ in pairs if lo <= e <= hi})
    if not starts or not ends:
        return lo, hi
    cs = [x for i, x in enumerate(starts) if i % max(1, len(starts) // 60) == 0] or starts[:1]
    ce = [x for i, x in enumerate(ends) if i % max(1, len(ends) // 60) == 0] or ends[-1:]
    best = None
    for st in cs:
        for en in ce:
            w = en - st + 1
            if w < step:
                continue
            sp = ov = 0
            for (rs, re), c in pairs:
                if re < st or rs > en:
                    continue
                ov += c
                if rs <= st and re >= en:
                    sp += c
            if not ov:
                continue
            f = sp / ov
            if f < floor:
                continue
            score = w * f
            if best is None or score > best[0]:
                best = (score, st, en, w, f)
    if not best:
        return lo, hi
    return best[1], best[2]


def spanning_profile_hist(iv, lo, hi, widths):
    """For each candidate width, the fraction of reads OVERLAPPING a window that cross it end
    to end -- computed from the interval histogram instead of re-fetching the BAMs.

    EQUIVALENCE TO THE OLD FETCH LOOP. The previous implementation, for each width, tiled the
    core into non-overlapping abutting windows `range(lo, hi - w + 2, w)` and, per window,
    called fetch(contig, st-1, en) and counted every returned primary read into `ovl`, plus
    those with start<=st and end>=en into `span`. pysam's fetch returns exactly the reads whose
    1-based interval [start, end] overlaps [st, en] (start <= en and end >= st) -- and those
    intervals, with multiplicity, ARE the keys of this histogram. So summing histogram counts
    over the same windows, with the same overlap and spanning tests, reproduces `ovl` and
    `span` exactly, including the fact that a read straddling a window boundary is counted once
    in each window it overlaps (the tiling is walked identically). No BAM is opened here."""
    pairs = list(iv.items())
    out = {}
    for w in widths:
        if hi - lo + 1 < w:
            continue
        ovl = span = 0
        for st in range(lo, hi - w + 2, w):
            en = st + w - 1
            for (rs, re), c in pairs:
                if re < st or rs > en:          # not overlapping this window
                    continue
                ovl += c
                if rs <= st and re >= en:       # spans the window end to end
                    span += c
        if ovl:
            out[w] = span / ovl
    return out


def spanning_profile(bams, contig, lo, hi, widths):
    """For each candidate width, the fraction of reads OVERLAPPING a window that actually
    cross it end to end.

    This replaces inferring spanning from a read-length percentile, which does not work: a
    read spans a window only if it is long enough AND starts in the right place, so among
    overlapping reads the spanning share is roughly (L-W+1)/(L+W-1). Asking for
    --spanning-target 0.7 at W=100 would need ~570 bp reads; measured, the old rule
    delivered 30-33% while claiming 70%.

    Measuring it directly also removes the need for a --window-max cap. On Cyclospora
    spanning fell below any sensible floor at 120 bp, before a 100 bp cap could bind. On
    Teladorsagia beta-tubulin the reads are 300 bp and 70% span the whole amplicon, so the
    cap was destroying real linkage across resistance codons 167/198/200.

    Kept as the public entry point. It now builds the contig's interval histogram once and
    defers to spanning_profile_hist(), so it no longer re-fetches every read per width; the
    returned fractions are identical (see spanning_profile_hist for the equivalence)."""
    iv = read_intervals(bams, contig)
    return spanning_profile_hist(iv, lo, hi, widths)


def covered(bam, min_depth, contigs=None):
    """Per-contig [lo, hi] reference span reaching >= min_depth, keyed by contig.

    `contigs` restricts the scan; when None it is read from the BAM index and is exactly the
    set of contigs that carry alignments. count_coverage() over a contig with no mapped reads
    returns all zeros and yields no positions at any depth, so omitting those contigs leaves
    the covered bounds byte-for-byte unchanged while skipping the hundreds of empty header
    contigs the reference carries (738 here, only 23 with reads)."""
    lo, hi = {}, {}
    with _require_pysam().AlignmentFile(bam, "rb") as fh:
        if contigs is None:
            contigs = [s.contig for s in fh.get_index_statistics() if s.mapped > 0]
        for contig in contigs:
            cov = fh.count_coverage(contig, quality_threshold=0)
            depth = [a + c + g + t for a, c, g, t in zip(*cov)]
            hits = [i + 1 for i, d in enumerate(depth) if d >= min_depth]
            if hits:
                lo[contig], hi[contig] = hits[0], hits[-1]
    return lo, hi


# NOTES -- optimization #4 (Rust kernel), a documented follow-up, NOT implemented here.
# No Rust toolchain is present in this environment (`which cargo`/`which rustc` are empty), so
# a native kernel would not build. Were one added later, the profitable surface is small and
# self-contained: the doubly-nested count in spannable_core() and spanning_profile_hist(), both
# of which reduce the same per-contig (start,end)->count histogram to per-window overlap/span
# tallies. A PyO3 extension exporting two functions --
#     core_scan(pairs, cs, ce, floor, step) -> (best_start, best_end)
#     width_profile(pairs, lo, hi, widths)  -> {width: fraction}
# taking the histogram as a pair of parallel i64 arrays (starts, ends) plus a counts array,
# would keep the Python I/O layer (pysam) untouched and move only the arithmetic. It must
# reproduce the tie-break here exactly (strict `>` on score, grid walked in sorted order) or
# the chosen window can shift. Given the histogram is already tiny (a few hundred keys per
# contig) after single-pass + subsampling, the Python cost is no longer the bottleneck, so
# this is low priority next to parallelising across libraries at the workflow level.


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ref", required=True)
    p.add_argument("--bams", nargs="+", required=True)
    p.add_argument("--sample", type=int, default=20, help="how many BAMs to scan")
    p.add_argument("--min-spanning", type=float, default=0.30,
                   help="minimum MEASURED fraction of overlapping reads that must span a "
                        "window. The widest candidate width meeting this is chosen. Replaces "
                        "--spanning-target, which inferred spanning from a read-length "
                        "percentile and delivered 30-33%% while claiming 70%%.")
    p.add_argument("--spanning-target", type=float, default=None,
                   help="deprecated: the old read-length percentile rule. Set it to restore "
                        "the previous behaviour exactly; otherwise --min-spanning is used.")
    p.add_argument("--min-depth", type=int, default=20,
                   help="depth required for a position to count as covered")
    p.add_argument("--window-min", type=int, default=40)
    p.add_argument("--window-max", type=int, default=0,
                   help="hard ceiling on window width. 0 = none, which is the default now that "
                        "width is chosen from measured spanning: the ceiling was redundant on "
                        "cohorts where spanning bound first, and harmful where reads span a "
                        "whole amplicon.")
    p.add_argument("--width-step", type=int, default=10,
                   help="granularity of the candidate widths searched")
    p.add_argument("--step", type=int, default=0, help="0 = non-overlapping")
    p.add_argument("--max-reads-per-window", type=int, default=0,
                   help="cap on the reads folded into each contig's interval histogram "
                        "(spanning is a fraction, so a bounded sample estimates it). 0 = no "
                        "cap, the default, which reproduces the un-subsampled windows exactly. "
                        "Subsampling is deterministic given --seed but is NOT full-depth "
                        "equivalent: it perturbs the set of read-start positions, so chosen "
                        "windows can shift -- the placement anchor moves ~1bp on deep contigs "
                        "(cascading a phase offset through the tiling) and borderline amplicon "
                        "widths, and hence window counts, can flip across the --min-spanning "
                        "floor. Seed-stability only sets in at a high cap (~100k on the "
                        "reference cohort). Treat capped output as an estimate, not a "
                        "replacement for the default full-depth run.")
    p.add_argument("--seed", type=int, default=0,
                   help="seed for --max-reads-per-window subsampling; fixes which reads are "
                        "kept so a capped run is reproducible.")
    p.add_argument("--threads", type=int, default=None,
                   help="worker threads for the per-contig window search (contigs are "
                        "independent). Defaults to $GALAXY_SLOTS, else 1. The chosen windows do "
                        "not depend on thread count.")
    p.add_argument("--out", default="-")
    a = p.parse_args(argv)

    # honour $GALAXY_SLOTS when --threads is not given, else fall back to single-threaded
    threads = a.threads or int(os.environ.get("GALAXY_SLOTS", "0") or 0) or 1

    # sample evenly across the cohort rather than taking the first N, so the estimate is
    # not dominated by whatever happens to sort first
    bams = a.bams
    if len(bams) > a.sample:
        stride = len(bams) / a.sample
        bams = [bams[int(i * stride)] for i in range(a.sample)]

    # Covered bounds per BAM, restricted to contigs that actually carry reads (see covered()).
    # Independent per BAM, so scanned in parallel; count_coverage releases the GIL.
    lo_all, hi_all = defaultdict(list), defaultdict(list)

    def scan(b):
        return covered(b, a.min_depth)

    if threads > 1 and len(bams) > 1:
        with ThreadPoolExecutor(max_workers=threads) as ex:
            scanned = list(ex.map(scan, bams))
    else:
        scanned = [scan(b) for b in bams]
    for lo, hi in scanned:
        for c in lo:
            lo_all[c].append(lo[c])
            hi_all[c].append(hi[c])

    if not lo_all:
        # No position in any sampled BAM reached --min-depth. The old code exited here only
        # when there were no aligned reads at all; in practice a cohort with reads but none
        # deep enough produced an empty BED. Keep the informative message.
        sys.exit("no aligned reads in the sampled BAMs")

    print(f"[define_windows] scanned {len(bams)} BAMs", file=sys.stderr)

    # covered bounds first -- the width search needs the interval it will tile
    bounds = {}
    for c in sorted(lo_all):
        L = sorted(lo_all[c]); H = sorted(hi_all[c])
        bounds[c] = (L[len(L) // 2], H[len(H) // 2])

    def process(c):
        """All window arithmetic for one contig, from a single per-contig histogram.

        Returns (contig, width, anchor, log). Pure w.r.t. shared state -- it reads only `bams`,
        `bounds[c]` and this contig's reads, and returns everything the parent assembles, so it
        is safe to run concurrently for different contigs."""
        lo_c, hi_c = bounds[c]
        iv = read_intervals(bams, c)                      # THE single pass over this contig
        iv = subsample_counter(iv, a.max_reads_per_window, a.seed)
        n = sum(iv.values())
        med = _nth_length(iv, n // 2)

        if a.spanning_target is not None:
            # deprecated path, kept so an old run can be reproduced exactly
            idx = max(0, int(n * (1.0 - a.spanning_target)) - 1)
            cap = a.window_max or 10 ** 9
            pv = _nth_length(iv, idx)
            w = max(a.window_min, min(cap, pv))
            pct = int((1 - a.spanning_target) * 100)
            log = (f"[define_windows]   {c}: n={n} median={med} p{pct}={pv} "
                   f"-> window {w} (legacy percentile rule)")
            # anchor defaults to the covered lo, which equals lo_c, so emission is unchanged
            return c, w, lo_c, log

        # PLACEMENT FIRST, then width. Tiling from the first covered base lets only the
        # earlier-starting strand span, which caps the reachable width before any width
        # search runs -- see spannable_core().
        lines = []
        s_core, e_core = spannable_core(iv, lo_c, hi_c, a.min_spanning, a.window_min)
        if (s_core, e_core) != (lo_c, hi_c):
            lines.append(f"[define_windows]   {c}: covered {lo_c}-{hi_c} -> spannable core "
                         f"{s_core}-{e_core} ({e_core - s_core + 1} bp)")

        top = min(hi_c - lo_c + 1, a.window_max or (hi_c - lo_c + 1))
        widths = list(range(a.window_min, top + 1, a.width_step))
        if top >= a.window_min and top not in widths:
            widths.append(top)                            # the core itself is always a candidate
        prof = spanning_profile_hist(iv, s_core, e_core, widths)   # measured in-core
        ok = [w for w, f in prof.items() if f >= a.min_spanning]
        w = max(ok) if ok else a.window_min
        shown = ", ".join(f"{ww}:{100*prof[ww]:.0f}%" for ww in sorted(prof)[:9])
        lines.append(
            f"[define_windows]   {c}: n={n} median_read={med} "
            f"spanning[{shown}] -> window {w}"
            f"{'' if ok else ' (NOTHING met --min-spanning; fell back to --window-min)'}")
        # The core sets the tiling PHASE, not the extent. Restricting emission to the core
        # throws away everything outside it: on Cyclospora that collapsed 28 tiled windows to
        # 8 and cost most of each amplicon. Keep the full covered interval (bounds[c]) and
        # simply start the tiling where reads begin to span.
        return c, w, s_core, "\n".join(lines)

    # Per-contig window selection is independent across contigs; run it in parallel. ex.map
    # preserves input order, so both the assembled dicts and the stderr log are in sorted
    # contig order regardless of how many threads ran -- the BED cannot depend on thread count.
    keys = sorted(bounds)
    if threads > 1 and len(keys) > 1:
        with ThreadPoolExecutor(max_workers=threads) as ex:
            results = list(ex.map(process, keys))
    else:
        results = [process(c) for c in keys]

    win = {}
    anchor = {}
    for c, w, anch, log in results:
        win[c] = w
        anchor[c] = anch
        print(log, file=sys.stderr)

    # Emit over the SPANNABLE CORE, not the covered interval. bounds[c] holds the core
    # once the placement search has run; it falls back to the median of the per-specimen
    # covered bounds, which is robust to one deep or one shallow BAM.
    fh = sys.stdout if a.out == "-" else open(a.out, "w")
    n = 0
    for c in sorted(lo_all):
        L = sorted(lo_all[c]); H = sorted(hi_all[c])
        lo, hi = bounds.get(c, (L[len(L) // 2], H[len(H) // 2]))
        w = win.get(c, a.window_min)
        step = a.step or w
        # Phase the tiling so a window boundary lands on the first position reads actually
        # span from. Walk back from there to `lo` in whole steps, then tile forward to `hi`.
        first = anchor.get(c, lo)
        while first - step >= lo:
            first -= step
        # Whole steps rarely land exactly on `lo`, and whatever is left over in front of the
        # first tile is simply dropped. On Cyclospora that cost the 5' end of four of seven
        # amplicons -- Nu_360i2 began at 97 instead of 7, losing 90 bp of real sequence and
        # three windows across the panel. Emit the remainder as a leading window whenever it
        # is at least --window-min wide.
        if first - lo >= a.window_min:
            print(f"{c}\t{lo-1}\t{first-1}\t{c}_W{lo:04d}", file=fh)
            n += 1
        for st in range(first, hi + 1, step):
            en = min(st + w - 1, hi)
            if en - st + 1 >= a.window_min:
                print(f"{c}\t{st-1}\t{en}\t{c}_W{st:04d}", file=fh)   # BED is 0-based
                n += 1
    if fh is not sys.stdout:
        fh.close()
    print(f"[define_windows] wrote {n} windows over {len(lo_all)} contigs", file=sys.stderr)


if __name__ == "__main__":
    main()
