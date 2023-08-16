# This dockerfile uses the Ubuntu image
# For install breakdancer automatically
# Author: Veronika Schusterbauer

FROM ubuntu:18.04

MAINTAINER Veronika Schusterbauer v.schusterbauer@gmx.at

RUN apt-get update \
    && apt-get -y install vim python perl git libncurses5-dev libncursesw5-dev  gcc g++ make cmake zlib1g-dev openssh-server wget

RUN apt-get -y install default-jre

WORKDIR /bin
#RUN wget https://github.com/PapenfussLab/clove/releases/download/v0.17/clove-0.17-jar-with-dependencies.jar
RUN wget https://github.com/PapenfussLab/clove/releases/download/v0.16/clove-0.16-jar-with-dependencies.jar
