"""
Deprecated backward-compatibility shim for cyclospora_pyeuk.
Use `pyeuk` instead.
"""

import warnings
from pyeuk import (
    __version__,
    __all__,
    generate_haplotype_sheet,
    generate_haplotype_sheet_from_assemblies,
    learn_de_novo_haplotypes,
    format_de_novo_haplotype_name,
    name_haplotype,
    parse_locus_name,
    PyEukDistanceEngine,
    CyclosporaClusterFinder,
    JunctionMicroAssembler,
    NanoporeAmpliconProcessor,
)
from pyeuk import amplicon

warnings.warn(
    "Importing from 'cyclospora_pyeuk' is deprecated and will be removed in a future release. "
    "Please import from 'pyeuk' instead.",
    DeprecationWarning,
    stacklevel=2,
)
