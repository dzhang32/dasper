FROM bioconductor/bioconductor_docker:RELEASE_3_13

LABEL authors="dyzhang32@gmail.com" \
      maintainer="dyzhang32@gmail.com" \
      description="Docker image based on the bioconductor/bioconductor_docker:RELEASE_3_13 with dasper installed"

WORKDIR /home/rstudio/dasper

COPY --chown=rstudio:rstudio . /home/rstudio/dasper

RUN Rscript -e "devtools::install('.', dependencies = TRUE, repos = BiocManager::repositories(), build_vignettes = TRUE)"
