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
    && rm -rf /var/lib/apt/lists/*

ENV RENV_PATHS_LIBRARY=/renv/library

WORKDIR /workflow
COPY renv.lock renv.lock
RUN Rscript -e "install.packages('renv'); options(BiocManager.version = '3.22'); renv::restore(prompt = FALSE)"

COPY . .
ENTRYPOINT ["Rscript", "notame-workflow.r"]