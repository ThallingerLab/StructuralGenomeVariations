FROM ubuntu:16.04

MAINTAINER Veronika Schusterbauer v.schusterbauer@gmx.at

RUN apt-get update && \
    apt-get install -y wget

WORKDIR /root

RUN wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz && \
     tar xzvf artbinmountrainier2016.06.05linux64.tgz && rm artbinmountrainier2016.06.05linux64.tgz
