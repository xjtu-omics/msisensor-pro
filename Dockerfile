FROM ubuntu:18.04

RUN apt-get update && apt-get install -y git make build-essential liblzma-dev libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev

RUN cd /tmp \
  && git clone https://github.com/xjtu-omics/msisensor-pro.git \
  && cd msisensor-pro \
  && ./INSTALL \
  && cp -r /tmp/msisensor-pro/binary/msisensor-pro /usr/bin/