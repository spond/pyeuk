# syntax=docker/dockerfile:1
# The syntax directive is required: the verification step below uses a heredoc, which is a
# BuildKit feature. The CI workflow uses buildx, so BuildKit is active there.
# Container image for PyEuk.
#
# Exists so a Galaxy tool, a Nextflow process or a Snakemake rule can declare a portable
# container instead of a local path. Until the Bioconda package is published, this is the only
# way to run PyEuk on a machine that did not build it. Once Bioconda publishes, the BioContainers
# image built from that recipe is the better choice for tool wrappers, and this image remains
# useful for anyone pinning a specific commit rather than a release.
#
# Dependency versions are pinned to the same values used in the Apptainer image that produced
# the published validation results, so a result reproduced from this container is comparable to
# those numbers rather than merely similar.

FROM python:3.11-slim-bookworm

LABEL org.opencontainers.image.title="PyEuk"
LABEL org.opencontainers.image.description="MLST/cgMLST typing, dropout-robust genetic distance estimation, and outbreak clustering for eukaryotic and microbial pathogens"
LABEL org.opencontainers.image.source="https://github.com/spond/pyeuk"
LABEL org.opencontainers.image.licenses="Apache-2.0"

RUN set -eu; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
        ca-certificates zlib1g libbz2-1.0 liblzma5 libcurl4 libssl3; \
    rm -rf /var/lib/apt/lists/*

RUN python -m pip install --no-cache-dir --upgrade \
        "pip==25.3" "setuptools==80.9.0" "wheel==0.45.1" && \
    python -m pip install --no-cache-dir \
        "numpy==2.1.3" "scipy==1.14.1" "pandas==2.2.3" \
        "scikit-learn==1.5.2" "numba==0.61.0"

WORKDIR /src
COPY . /src
RUN python -m pip install --no-cache-dir ".[amplicon]" && rm -rf /src

# Fail the build rather than ship an image that cannot do what it claims.
RUN python - <<'PYEOF'
import inspect
import cyclospora_pyeuk
from cyclospora_pyeuk.clustering import CyclosporaClusterFinder
from cyclospora_pyeuk.distance_engine import PyEukDistanceEngine
import cyclospora_pyeuk.amplicon as amplicon
import pysam

fp = inspect.signature(CyclosporaClusterFinder.find_clusters).parameters
ep = inspect.signature(PyEukDistanceEngine.__init__).parameters
dp = inspect.signature(PyEukDistanceEngine.compute_revised_wibs_matrix).parameters

assert "cut_mode" in fp, "distance-cut mode missing"
assert "linkage_method" in fp, "linkage_method missing"
assert hasattr(CyclosporaClusterFinder, "suggest_linkage_threshold")
assert "weight_mode" in ep or "weight_mode" in dp, "weight_mode missing"
assert "project_psd" in ep or "project_psd" in dp, "project_psd missing"
assert set(amplicon.__all__) == {"define_windows", "window_haplotypes", "build_sheet"}
print("pyeuk", cyclospora_pyeuk.__version__, "| pysam", pysam.__version__, "| all checks passed")
PYEOF

WORKDIR /data
ENTRYPOINT ["pyeuk"]
CMD ["--help"]
