FROM rocker/r-ver:4.5.2

RUN apt-get update && apt-get install -y \
    git \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libglpk-dev \
    poppler-utils \
    && rm -rf /var/lib/apt/lists/*

ENV RENV_PATHS_LIBRARY=/renv/library
ENV RENV_CONFIG_AUTOLOADER_ENABLED=FALSE
ENV R_LIBS_SITE=/renv/library

WORKDIR /workflow
COPY renv.lock renv.lock
# renv::restore() resolves a GitHub-sourced package's install-order dependencies by
# fetching its DESCRIPTION live from GitHub, ignoring whatever this lockfile records
# for it -- so WaveICA (v1)'s own (upstream, not ours) DESCRIPTION under-declaring
# pROC/plsdepot/etc. as Imports (its NAMESPACE imports them regardless) means restore()
# can schedule WaveICA's install before theirs are present. Installing them explicitly
# first guarantees they're already on the library path by the time WaveICA builds,
# regardless of restore()'s internal ordering for that one broken graph edge.
RUN Rscript -e "install.packages('renv'); options(BiocManager.version = '3.22'); renv::install(c('pROC','plsdepot','fdrtool','scatterplot3d','ggfortify','RColorBrewer','ggplot2','gridExtra'), prompt = FALSE); renv::restore(prompt = FALSE)"

COPY . .
ENTRYPOINT ["Rscript", "notame-workflow.r"]