FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y \
    wget \
    bzip2 \
    gcc \
    make \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    python3 \
    python3-pip && \
    rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2 && \
    tar -xvjf htslib-1.12.tar.bz2 && \
    cd htslib-1.12 && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -rf htslib-1.12 htslib-1.12.tar.bz2

RUN pip3 install pysam
CMD ["bash"]
