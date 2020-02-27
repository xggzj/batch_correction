FROM r-base:latest

RUN apt-get -y install libcurl4-gnutls-dev libxml2-dev libssl-dev

RUN mkdir -p /code
RUN R -e "install.packages(c('BiocManager', 'devtools', 'dplyr', 'ggplot2', 'umap'))"
RUN R -e "BiocManager::install('flowCore')"
RUN R -e "devtools::install_github('immunogenomics/harmony')"

COPY . /code
WORKDIR /code

CMD ["Rscript", "R/pipeline.R"]






