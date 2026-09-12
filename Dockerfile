# Single image carrying every tool mpCGC calls.
#
#   docker build -t mpcgc:1.0.0 .
#   nextflow run . -profile docker --container mpcgc:1.0.0 ...
#
# The docker and singularity profiles set process.container to --container, so
# one image serves every process. Per-process biocontainer images are still
# declared in modules/local/*.nf for anyone who prefers granular images.
#
# The reference databases are NOT baked in: they are ~8 GB and are versioned
# independently of the code. Build them once with --step download_db and mount
# the directory.
#
# R is deliberately absent. The only R step, GGTREE_FIGURE, is optional and runs
# from the stock bioconductor-ggtree biocontainer (--container_ggtree), which
# keeps an R and Bioconductor stack out of this image.

FROM condaforge/miniforge3:24.9.2-0

LABEL org.opencontainers.image.title="mpCGC" \
      org.opencontainers.image.description="marine polysaccharide CAZyme gene clusterer" \
      org.opencontainers.image.source="https://github.com/AaronAOliver/mpCGC" \
      org.opencontainers.image.licenses="MIT"

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update \
    && apt-get install -y --no-install-recommends procps ca-certificates curl \
    && rm -rf /var/lib/apt/lists/*

COPY env/ /tmp/env/

# One environment holding all of it, installed into the base env so nothing
# needs activating inside a Nextflow task.
RUN mamba install -y -n base \
        -c conda-forge -c bioconda \
        python=3.12 \
        dbcan=5.2.9 \
        pyrodigal=3.7.1 \
        pyhmmer=0.12.0 \
        diamond=2.1.9 \
        muscle=5.1 \
        fasttree=2.1.11 \
        numpy'>=1.26' \
        scipy'>=1.11' \
        pandas'>=2.1' \
        matplotlib-base'>=3.8' \
        seaborn'>=0.13' \
        python-igraph'>=0.11' \
        leidenalg'>=0.10' \
    && mamba clean -afy \
    && rm -rf /tmp/env

# FastTreeMP is the OpenMP build the pipeline asks for by default; some conda
# packages ship only FastTree, so provide the name either way.
RUN if [ ! -x /opt/conda/bin/FastTreeMP ] && [ -x /opt/conda/bin/FastTree ]; then \
        ln -s /opt/conda/bin/FastTree /opt/conda/bin/FastTreeMP; \
    fi

ENV PATH=/opt/conda/bin:$PATH \
    PYTHONDONTWRITEBYTECODE=1 \
    MPLBACKEND=Agg

RUN run_dbcan --help > /dev/null \
    && diamond --version \
    && muscle -version \
    && python -c "import numpy, scipy, pandas, igraph, leidenalg, matplotlib"

CMD ["/bin/bash"]
