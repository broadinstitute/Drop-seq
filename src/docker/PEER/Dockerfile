# NOTE: Apple Silicon Docker runs may show small numeric differences relative to native linux x64 systems.
# Factors are highly correlated (most > 0.995), in tests the least correlated factor was 0.985.

# Assume that the docker file is built from the root of the public repository:
# docker build -f src/docker/PEER/Dockerfile -t PEER .

FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    software-properties-common \
    wget \
    curl \
    cmake \
    r-base-core \
    build-essential \
    gcc-5 \
    g++-5 \
    gfortran-5 \
    && rm -rf /var/lib/apt/lists/*

RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 50 \
    --slave /usr/bin/g++ g++ /usr/bin/g++-5

RUN curl -L -o /tmp/R_peer_source_1.3.tgz \
    https://raw.githubusercontent.com/broadinstitute/gtex-pipeline/49e544617123a9e81170f25b278c85aa122db890/qtl/R_peer_source_1.3.tgz

RUN ls -lh /tmp/R_peer_source_1.3.tgz

RUN R CMD INSTALL /tmp/R_peer_source_1.3.tgz

RUN rm /tmp/R_peer_source_1.3.tgz

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

RUN R -e "install.packages(c('argparser','qvalue'), repos='http://cran.us.r-project.org')"

COPY src/R/PEER/run_PEER.R /scripts/run_PEER.R

ENTRYPOINT ["Rscript", "/scripts/run_PEER.R"]