# This dockerfile uses the Ubuntu image
# For install breakdancer automatically
# Author: Veronika Schusterbauer

FROM ubuntu:14.04

MAINTAINER Veronika Schusterbauer v.schusterbauer@gmx.at

EXPOSE 22

RUN apt-get update \
    && apt-get -y install vim python perl git libncurses5-dev libncursesw5-dev  gcc g++ make cmake zlib1g-dev openssh-server wget

WORKDIR /root
RUN wget https://github.com/genome/breakdancer/archive/refs/tags/v1.4.5.tar.gz && \
     tar xzvf v1.4.5.tar.gz && rm v1.4.5.tar.gz && \
     mkdir breakdancer-1.4.5/build

RUN cd /root/breakdancer-1.4.5/build && cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/usr/local && make && make install

RUN apt-get install -y autoconf libncurses5-dev libbz2-dev liblzma-dev

RUN cd /root && \
  wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && \
  tar vxjf samtools-1.11.tar.bz2 && \
  rm samtools-1.11.tar.bz2 && \
  cd samtools-1.11 && autoheader && autoconf -Wno-syntax && \
  ./configure && make

ENV PATH=/root/breakdancer-1.4.5/build/bin/:/root/samtools-1.11/:$PATH

 ## Add for /root/breakdancer/perl/bam2cfg.pl ##

RUN apt-get -y install libstatistics-descriptive-perl libgd3 libgd-perl cpanminus && cpanm GD::Graph GD::Graph::histogram
WORKDIR /root
