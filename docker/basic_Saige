FROM ubuntu:16.04

LABEL maintainer="wzhou@broadinstitute.org"

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ENV R_VERSION 3.5.1

# Install dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    curl \
    libboost-all-dev \
    libbz2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    libopenblas-dev \
    libpcre3-dev \
    libpng-dev \
    libreadline-dev \
    libssl-dev \
    libz-dev \
    python3-pip \
    python3-setuptools \
    unzip \
&& pip3 install cget \
&& rm -rf /var/lib/apt/lists/*

WORKDIR /tmp

# Install R from source
RUN curl -O https://cloud.r-project.org/src/base/R-3/R-${R_VERSION}.tar.gz \
    && tar xvzf R-${R_VERSION}.tar.gz

WORKDIR /tmp/R-${R_VERSION}

RUN ./configure --with-x=no --with-blas="-lopenblas" \
    && make \
    && mkdir -p /usr/local/lib/R/lib \
    && make install

WORKDIR /tmp

# Cleanup R source code
RUN rm -rf /tmp/R-${R_VERSION}*

# Download latest commit of SAIGE from GitHub
RUN curl -OL https://github.com/weizhouUMICH/SAIGE/archive/master.zip \
    && unzip master.zip \
    && mv SAIGE* SAIGE

# Install R packages and SAIGE
RUN Rscript SAIGE/extdata/install_packages.R \
    && R CMD INSTALL SAIGE

# Make R scripts to /usr/local/bin to make executable and cleanup SAIGE source code
RUN mv SAIGE/extdata/step1_fitNULLGLMM.R SAIGE/extdata/step2_SPAtests.R SAIGE/extdata/createSparseGRM.R /usr/local/bin/ \
    && rm -rf SAIGE \
    && rm master.zip

WORKDIR /
