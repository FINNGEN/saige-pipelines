FROM eu.gcr.io/finngen-refinery-dev/bioinformatics:0.7

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8


RUN apt-get update &&  apt-get install wget make curl gcc g++ cmake gfortran libreadline6-dev libboost-all-dev libopenblas-dev \
default-jre lsb-core libcairo2-dev libpcre2-dev libpcre3-dev libssl-dev libcurl4-openssl-dev libz-dev libbz2-dev liblzma-dev zlib1g-dev --yes && apt-get install python3-pip -y

RUN apt-get install moreutils  zip unzip -y

RUN wget https://github.com/rgcgithub/regenie/releases/download/v3.0/regenie_v3.0.gz_x86_64_Linux_mkl.zip && unzip  regenie_v3.0.gz_x86_64_Linux_mkl.zip && mv regenie_v3.0.gz_x86_64_Linux_mkl /usr/bin/regenie

RUN pip3 install --no-cache numpy scipy pandas cget

RUN pip3 install csvkit
ADD scripts ../scripts