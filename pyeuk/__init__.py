"""
pyeuk: Modern Python package for eukaryotic amplicon MLST genotyping,
Eukaryotyping ensemble distance engine, and outbreak cluster finder.
Supports Illumina short-reads, Oxford Nanopore long-reads (ONT), and PacBio HiFi.
"""

from .haplotype_sheet import (
    generate_haplotype_sheet,
    generate_haplotype_sheet_from_assemblies,
    learn_de_novo_haplotypes,
    format_de_novo_haplotype_name,
)
from .naming import (
    name_haplotype,
    parse_locus_name,
)
from .distance_engine import PyEukDistanceEngine
from .clustering import CyclosporaClusterFinder
from .micro_assembly import JunctionMicroAssembler
from .ont_processor import NanoporeAmpliconProcessor

__version__ = "0.6.0"
__all__ = [
    "generate_haplotype_sheet",
    "generate_haplotype_sheet_from_assemblies",
    "learn_de_novo_haplotypes",
    "format_de_novo_haplotype_name",
    "name_haplotype",
    "parse_locus_name",
    "PyEukDistanceEngine",
    "CyclosporaClusterFinder",
    "JunctionMicroAssembler",
    "NanoporeAmpliconProcessor",
]
