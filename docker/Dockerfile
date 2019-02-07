FROM rocker/r-base

# Metadata
LABEL description="Docker image containing all dependencies for PheWAS pipeline"

# Maintainer
LABEL Phil Palmer <phil@lifebit.ai>

# Install plink
RUN wget http://s3.amazonaws.com/plink1-assets/dev/plink_linux_x86_64.zip && \
    unzip *.zip && \
    mv plink /usr/local/bin/

# Install R dependencies
RUN apt-get update
RUN apt-get install -t unstable -y libcurl4-openssl-dev
RUN apt-get install -t unstable -y libssl-dev
RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "library(devtools)"
RUN Rscript -e "devtools::install_github('PheWAS/PheWAS')"
RUN apt-get install -y procps

RUN mkdir /data && \
    mkdir /data/scripts/

COPY ./R/ /data/scripts/

RUN chmod -R u+x /data/scripts

ENV PATH=$PATH:/data/scripts/

WORKDIR /data/scripts/

