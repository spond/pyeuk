"""
cyclospora_pyeuk: Modern Python package for CDC Cyclospora cayetanensis MLST Genotyping Workflow.
Supports Illumina short-reads, Oxford Nanopore long-reads (ONT), and PacBio HiFi.
"""

from .haplotype_sheet import (
    generate_haplotype_sheet,
    generate_haplotype_sheet_from_assemblies,
    learn_de_novo_haplotypes,
    format_de_novo_haplotype_name,
)
from .distance_engine import PyEukDistanceEngine
from .clustering import CyclosporaClusterFinder
from .micro_assembly import JunctionMicroAssembler
from .ont_processor import NanoporeAmpliconProcessor

__version__ = "0.4.0"
__all__ = [
    "generate_haplotype_sheet",
    "generate_haplotype_sheet_from_assemblies",
    "learn_de_novo_haplotypes",
    "format_de_novo_haplotype_name",
    "PyEukDistanceEngine",
    "CyclosporaClusterFinder",
    "JunctionMicroAssembler",
    "NanoporeAmpliconProcessor",
]
