# use the ubuntu base image
FROM ubuntu:22.04

MAINTAINER Veronika Schusterbauer veronika.schusterbauer@bisy.at

# install required packages using apt
RUN apt-get update && apt-get install -y \
#    autoconf \
#    build-essential \
#    cmake \
#    g++ \
#    gfortran \
#    git \
    libcurl4-gnutls-dev \
    hdf5-tools \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libbz2-dev \
    libdeflate-dev \
    libhdf5-dev \
    libncurses-dev \
    liblzma-dev \
    pkg-config \
    zlib1g-dev \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# set environment
ENV BOOST_ROOT /usr

# install delly_1.1.6
RUN cd /bin && \
    wget https://github.com/dellytools/delly/releases/download/v1.1.6/delly_v1.1.6_linux_x86_64bit && \
    chmod 777 delly_v1.1.6_linux_x86_64bit \

# Workdir
WORKDIR /root/

# Add Delly to PATH
ENV PATH="/opt/delly/bin:${PATH}"

# by default /bin/sh is executed
CMD ["/bin/sh"]
