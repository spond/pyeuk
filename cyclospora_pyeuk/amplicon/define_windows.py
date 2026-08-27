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

Usage:
  define_windows.py --ref panel.fa --bams a.bam b.bam ... [--sample 20]
                    [--spanning-target 0.7] [--min-depth 20] [--out windows.bed]
"""
import argparse
import sys
from collections import defaultdict

# pysam is an OPTIONAL dependency, declared under the `amplicon` extra. Importing it at module
# scope would make `import cyclospora_pyeuk` fail for anyone who installed the core package to
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
                "  pip install 'cyclospora_pyeuk[amplicon]'"
            ) from exc
        _pysam = _p
    return _pysam


def block_lengths(bam):
    """Aligned block length per read, keyed by contig.

    Per contig, not pooled: amplicons differ in depth and in fragment length, so one
    cohort-wide window size is set by whichever amplicon is deepest and is then far too
    wide for the shallow ones, which lose their locus entirely. CDC's own PART windows
    differ in size per marker for the same reason."""
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
    instead of a few million reads."""
    from collections import Counter
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
    distinct read intervals."""
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
    cap was destroying real linkage across resistance codons 167/198/200."""
    out = {}
    for w in widths:
        if hi - lo + 1 < w:
            continue
        ovl = span = 0
        for bam in bams:
            with _require_pysam().AlignmentFile(bam, "rb") as fh:
                if contig not in fh.references:
                    continue
                for st in range(lo, hi - w + 2, w):
                    en = st + w - 1
                    for r in fh.fetch(contig, st - 1, en):
                        if (r.is_unmapped or r.is_secondary or r.is_supplementary
                                or not r.reference_length):
                            continue
                        ovl += 1
                        if r.reference_start + 1 <= st and (r.reference_end or 0) >= en:
                            span += 1
        if ovl:
            out[w] = span / ovl
    return out


def covered(bam, min_depth):
    lo, hi = {}, {}
    with _require_pysam().AlignmentFile(bam, "rb") as fh:
        for contig in fh.references:
            cov = fh.count_coverage(contig, quality_threshold=0)
            depth = [a + c + g + t for a, c, g, t in zip(*cov)]
            hits = [i + 1 for i, d in enumerate(depth) if d >= min_depth]
            if hits:
                lo[contig], hi[contig] = hits[0], hits[-1]
    return lo, hi


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
    p.add_argument("--out", default="-")
    a = p.parse_args(argv)

    # sample evenly across the cohort rather than taking the first N, so the estimate is
    # not dominated by whatever happens to sort first
    bams = a.bams
    if len(bams) > a.sample:
        stride = len(bams) / a.sample
        bams = [bams[int(i * stride)] for i in range(a.sample)]

    lens = defaultdict(list)
    lo_all, hi_all = defaultdict(list), defaultdict(list)
    for b in bams:
        for c, v in block_lengths(b).items():
            lens[c] += v
        lo, hi = covered(b, a.min_depth)
        for c in lo:
            lo_all[c].append(lo[c])
            hi_all[c].append(hi[c])
    if not lens:
        sys.exit("no aligned reads in the sampled BAMs")

    print(f"[define_windows] scanned {len(bams)} BAMs", file=sys.stderr)

    # covered bounds first -- the width search needs the interval it will tile
    bounds = {}
    for c in sorted(lo_all):
        L = sorted(lo_all[c]); H = sorted(hi_all[c])
        bounds[c] = (L[len(L) // 2], H[len(H) // 2])

    win = {}
    anchor = {}
    for c, v in sorted(lens.items()):
        v.sort()
        med = v[len(v) // 2]
        if c not in bounds:
            win[c] = a.window_min
            continue
        lo_c, hi_c = bounds[c]

        if a.spanning_target is not None:
            # deprecated path, kept so an old run can be reproduced exactly
            idx = max(0, int(len(v) * (1.0 - a.spanning_target)) - 1)
            cap = a.window_max or 10 ** 9
            win[c] = max(a.window_min, min(cap, v[idx]))
            pct = int((1 - a.spanning_target) * 100)
            print(f"[define_windows]   {c}: n={len(v)} median={med} p{pct}={v[idx]} "
                  f"-> window {win[c]} (legacy percentile rule)", file=sys.stderr)
            continue

        # PLACEMENT FIRST, then width. Tiling from the first covered base lets only the
        # earlier-starting strand span, which caps the reachable width before any width
        # search runs -- see spannable_core().
        iv = read_intervals(bams, c)
        s_core, e_core = spannable_core(iv, lo_c, hi_c, a.min_spanning, a.window_min)
        if (s_core, e_core) != (lo_c, hi_c):
            print(f"[define_windows]   {c}: covered {lo_c}-{hi_c} -> spannable core "
                  f"{s_core}-{e_core} ({e_core - s_core + 1} bp)", file=sys.stderr)
        # The core sets the tiling PHASE, not the extent. Restricting emission to the core
        # throws away everything outside it: on Cyclospora that collapsed 28 tiled windows to
        # 8 and cost most of each amplicon. Keep the full covered interval and simply start
        # the tiling where reads begin to span.
        bounds[c] = (lo_c, hi_c)
        anchor[c] = s_core

        top = min(hi_c - lo_c + 1, a.window_max or (hi_c - lo_c + 1))
        widths = list(range(a.window_min, top + 1, a.width_step))
        if top >= a.window_min and top not in widths:
            widths.append(top)          # the core itself is always a candidate
        prof = spanning_profile(bams, c, s_core, e_core, widths)  # measured in-core
        ok = [w for w, f in prof.items() if f >= a.min_spanning]
        win[c] = max(ok) if ok else a.window_min
        shown = ", ".join(f"{w}:{100*prof[w]:.0f}%" for w in sorted(prof)[:9])
        print(f"[define_windows]   {c}: n={len(v)} median_read={med} "
              f"spanning[{shown}] -> window {win[c]}"
              f"{'' if ok else ' (NOTHING met --min-spanning; fell back to --window-min)'}",
              file=sys.stderr)

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
