FROM ubuntu:18.04

MAINTAINER Veronika Schusterbauer v.schusterbauer@gmx.at

WORKDIR /home

RUN apt-get update && apt install -y wget build-essential cmake make libz-dev libbz2-dev

RUN apt-get install -y liblzma-dev libcurl4-openssl-dev

RUN wget https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2 && \
  tar -vxjf bcftools-1.11.tar.bz2 && \
  rm bcftools-1.11.tar.bz2 && \
  cd bcftools-1.11 && make

RUN wget https://github.com/iqbal-lab/cortex/archive/refs/tags/v1.0.5.21.tar.gz && \
    tar xzvf v1.0.5.21.tar.gz && rm v1.0.5.21.tar.gz && \
    cd cortex-1.0.5.21 && bash install.sh

RUN  cd cortex-1.0.5.21 && make MAXK=31 cortex_var && make MAXK=63 cortex_var

