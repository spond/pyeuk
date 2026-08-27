#!/usr/bin/env python3
"""Turn per-specimen window-haplotype calls into a PyEuk sheet.

Species- and panel-agnostic: it only reads the long-format TSVs that
window_haplotypes.py writes.

A column is one observed haplotype string in one window. Its name is
<window>_Hap_<i>, because PyEuk's parse_locus_name() derives locus windows by
stripping exactly that suffix, so each window becomes its own locus window --
matching how CDC's PART_A / PART_B behave. The content-derived name lives in
haplotype_map.tsv.

Two properties fall out of naming by content rather than by catalogue:

  * "amplified and reference-identical" is the ordinary haplotype "=", so it is a
    positive observation. No special CALLED column is needed, and a specimen matching
    the reference is not confusable with a dropout, which has no row at all.
  * a mixture is several haplotypes in one window, which the long format and the
    presence/absence sheet both represent natively.

Usage: build_window_sheet.py <calls_dir> <outdir> [--min-freq F] [--min-reads N]
"""
import argparse
import csv
import glob
import os
from collections import defaultdict

def _parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("calls_dir")
    p.add_argument("outdir")
    p.add_argument("--min-freq", type=float, default=0.0,
                   help="extra haplotype-frequency filter on top of whatever the caller used")
    p.add_argument("--min-reads", type=int, default=0)
    p.add_argument("--min-span", type=int, default=0,
                   help="minimum total fully-spanning reads for a window to count as CALLED for a "
                        "specimen. 0 disables (the caller's own --min-span already applied). Set >0 "
                        "only when the caller was run permissively, so the gate can be swept here "
                        "without re-reading a BAM -- spanning is a per-window total the caller "
                        "already emits, so every gate is a pure post-filter on its output.")
    p.add_argument("--min-maf", type=float, default=0.0,
                   help="drop a haplotype column unless it is carried by at least this fraction "
                        "of the specimens called at its window (and at most 1-min_maf). 0 disables.")
    return p


def main(argv=None):
    a = _parser().parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)

    # specimen -> window -> {haplotype: (reads, freq)}
    data = defaultdict(lambda: defaultdict(dict))
    spanning = defaultdict(dict)
    specimens, windows = [], {}

    for path in sorted(glob.glob(f"{a.calls_dir}/*.tsv")):
        for r in csv.DictReader(open(path), delimiter="\t"):
            s, win, hap = r["specimen"], r["window"], r["haplotype"]
            if s not in specimens:
                specimens.append(s)
            windows[win] = (r["locus"], int(r["start"]), int(r["end"]))
            spanning[s][win] = int(r["spanning"])
            if hap == "NOT_CALLED":
                continue
            # Window-level gate first: below min_span the window is NOT CALLED for this specimen,
            # so it must not contribute a haplotype AND must not count toward n_called below.
            # Applying it here rather than in the caller is what makes a gate sweep cheap -- one
            # permissive caller pass over the BAMs serves the entire grid.
            if a.min_span and int(r["spanning"]) < a.min_span:
                continue
            reads, freq = int(r["reads"]), float(r["freq"])
            if reads < a.min_reads or freq < a.min_freq:
                continue
            data[s][win][hap] = (reads, freq)

    # stable column order: window by locus/start, haplotype by cohort frequency then name
    order = sorted(windows, key=lambda w: (windows[w][0], windows[w][1]))
    seen = defaultdict(lambda: defaultdict(int))
    for s in specimens:
        for win, haps in data[s].items():
            for h in haps:
                seen[win][h] += 1

    # Minor-allele filtering, before any distance is computed.
    #
    # PyEuk weights a column w = 1/sqrt(p(1-p)) where p is the fraction of specimens CALLED at
    # that window which carry the haplotype. That function is minimised at p=0.5 and diverges as
    # p approaches 0 or 1, so it weights RARE alleles most heavily. On a 153-specimen Cyclospora
    # sheet, 17 columns carried by <=2 specimens held 79.6% of all weight, while the ten balanced
    # columns carrying the entire outbreak signal held 8.7%. A private singleton outweighed a
    # perfectly discriminating allele 6:1.
    #
    # The damage is not merely dilution. When a private singleton is the ONLY weight-bearing
    # column two specimens share, the numerator equals the denominator and their distance
    # saturates at exactly 1.0 -- observed between two specimens of the SAME outbreak. Capping the
    # weight cannot fix that, because the weight cancels; the column has to go.
    #
    # 0.05 is the conventional minor-allele cutoff in population genetics, chosen a priori rather
    # than fitted. Note p is a fraction of SPECIMENS, not of reads, so this does not touch
    # within-specimen minor haplotypes -- a 25% minor component in a mixture is untouched.
    n_called = {}
    for win in order:
        n_called[win] = sum(1 for s in specimens if data[s].get(win))

    dropped = []
    if a.min_maf > 0:
        for win in order:
            n = n_called[win]
            if not n:
                continue
            for h in list(seen[win]):
                frac = seen[win][h] / n
                if frac < a.min_maf or frac > 1 - a.min_maf:
                    dropped.append((win, h, seen[win][h], round(frac, 4)))
                    del seen[win][h]

    cols, mapping = [], []
    idx = {}
    for win in order:
        haps = sorted(seen[win], key=lambda h: (-seen[win][h], h))
        idx[win] = {h: i + 1 for i, h in enumerate(haps)}
        for h in haps:
            c = f"{win}_Hap_{idx[win][h]}"
            cols.append((win, h, c))
            mapping.append([c, win, windows[win][0], windows[win][1], windows[win][2],
                            h, seen[win][h]])

    with open(f"{a.outdir}/sheet.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["Seq_ID"] + [c for _, _, c in cols])
        for s in specimens:
            w.writerow([s] + ["X" if h in data[s].get(win, {}) else "" for win, h, _ in cols])

    with open(f"{a.outdir}/haplotype_map.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["column", "window", "locus", "start", "end", "haplotype", "n_specimens"])
        w.writerows(mapping)

    # frequency is the reason for doing any of this, so it is kept, not discarded
    with open(f"{a.outdir}/calls_long.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["Seq_ID", "window", "column", "haplotype", "reads", "freq", "spanning"])
        for s in specimens:
            for win in order:
                for h, (reads, freq) in sorted(data[s].get(win, {}).items(),
                                               key=lambda kv: -kv[1][1]):
                    # the long format is the record of what was OBSERVED, so MAF-filtered
                    # haplotypes stay here and are simply marked as having no sheet column
                    col = f"{win}_Hap_{idx[win][h]}" if h in idx.get(win, {}) else "FILTERED"
                    w.writerow([s, win, col, h, reads, freq, spanning[s].get(win, 0)])

    with open(f"{a.outdir}/dropped_columns.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["window", "haplotype", "n_specimens", "fraction_of_called"])
        w.writerows(dropped)

    mixed = sum(1 for s in specimens for win in data[s] if len(data[s][win]) > 1)
    if a.min_maf > 0:
        print(f"minor-allele filter at {a.min_maf}: dropped {len(dropped)} columns "
              f"(see dropped_columns.tsv)")
    print(f"specimens {len(specimens)}   windows {len(order)}   columns {len(cols)}")
    print(f"windows with >1 haplotype in a specimen (mixtures): {mixed}")
    print(f"{'window':<18}{'haps':>6}{'called':>8}{'ref-identical':>15}")
    for win in order:
        ncall = sum(1 for s in specimens if data[s].get(win))
        nref = sum(1 for s in specimens if "=" in data[s].get(win, {}))
        print(f"{win:<18}{len(seen[win]):>6}{ncall:>8}{nref:>15}")


if __name__ == "__main__":
    main()
