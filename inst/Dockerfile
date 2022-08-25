FROM bioconductor/bioconductor_docker:latest
ADD . ~/ADTnorm/
WORKDIR ~/
RUN R -q -e 'install.packages("devtools")'
RUN R -q -e 'devtools::install("./ADTnorm", dependencies = TRUE)'
WORKDIR /home/rstudio
