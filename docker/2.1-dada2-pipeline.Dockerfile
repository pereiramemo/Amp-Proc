FROM mambaorg/micromamba:1.5.10-noble

COPY --chown=$MAMBA_USER:$MAMBA_USER docker/resources/1.2-quality-check.requirements.yml /tmp/conda.yml
RUN micromamba install -y -n base -f /tmp/conda.yml \
    && micromamba install -y -n base conda-forge::procps-ng \
    && micromamba run -n base Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")' \
    && micromamba run -n base Rscript -e 'BiocManager::install("GenomeInfoDbData", update=FALSE, ask=FALSE)' \
    && micromamba run -n base Rscript -e 'library(dada2); library(tidyverse)' \
    && micromamba env export --name base --explicit > environment.lock \
    && echo ">> CONDA_LOCK_START" \
    && cat environment.lock \
    && echo "<< CONDA_LOCK_END" \
    && micromamba clean -a -y

USER root
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"