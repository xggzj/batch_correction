FROM rocker/tidyverse:3.6.1

RUN mkdir -p /code
RUN R -e "devtools::install_github('immunogenomics/harmony')"
RUN R -e "install.packages(c('BiocManager', 'dplyr', 'ggplot2', 'umap'))"
RUN R -e "BiocManager::install('flowCore')"

COPY . /code
WORKDIR /code

#CMD ["Rscript", "R/pipeline.R"]






