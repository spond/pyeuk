"""
Amplicon front end: BAMs in, haplotype sheet out.

This is the half of the chain that runs before the distance engine. It turns aligned reads
into the specimen x haplotype sheet the rest of PyEuk consumes, without a curated haplotype
catalogue.

    define_windows      choose the analysis windows from the cohort's own reads
    window_haplotypes   read a haplotype off each single spanning read
    build_sheet         assemble the presence/absence sheet and the long-format calls

The unit is a WINDOW HAPLOTYPE: the string a single read carries across an interval, with the
read required to cover every base of that interval or be discarded. Linkage between positions
is therefore observed on one molecule rather than inferred across molecules. That distinction
is the point of the approach. Per-site variant calling records each mutation independently of
the molecule carrying it, and so cannot separate one strain carrying N mutations from a mixture
in which a second strain contributes some of them.

Why these live beside the distance engine rather than in a separate project: the haplotype
names produced here are parsed by pyeuk.naming on the way in, and that contract has
already broken once when the two halves were versioned separately.

pysam is required here and only here. It is declared as the `amplicon` extra and imported on
first use, so the core package installs and imports without it:

    pip install 'pyeuk[amplicon]'
"""

from . import build_sheet, define_windows, window_haplotypes  # noqa: F401

__all__ = ["define_windows", "window_haplotypes", "build_sheet"]
