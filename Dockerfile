FROM ubuntu:24.04 AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    autoconf \
    automake \
    make \
    g++ \
    liblzma-dev \
    libbz2-dev \
    zlib1g-dev \
    libcurl4-gnutls-dev \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY Makefile VERSION /app/
COPY vendor /app/vendor/
COPY cpp /app/cpp/
RUN make

FROM ubuntu:24.04 AS final

RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 \
    libcurl4-gnutls-dev \
    libbz2-1.0 \
    liblzma5 \
    zlib1g \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /app/msisensor-pro /usr/local/bin/msisensor-pro

CMD ["msisensor-pro"]
