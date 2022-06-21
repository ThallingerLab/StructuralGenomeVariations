FROM alpine:3.16.0

MAINTAINER Veronika Schusterbauer v.schusterbauer@gmx.at

RUN mkdir -p /opt/delly/bin
WORKDIR /opt/delly/bin

RUN cd /opt/delly/bin&& \
  wget https://github.com/dellytools/delly/releases/download/v1.0.3/delly_v1.0.3_linux_x86_64bit && \
  chmod +rwx delly_v1.0.3_linux_x86_64bit
