FROM ubuntu:20.04

ARG DEBIAN_FRONTEND="noninteractive"
ENV LANG en_US.utf8
ENV MAMBA_EXE=/opt/micromamba/bin/micromamba
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV PATH /opt/micromamba/bin:$PATH
ENV TZ=America/LosAngeles

RUN apt-get -y update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        default-jdk \
        default-jre \
        gcc g++ \
        gfortran \
        libarchive13 \
        libbz2-dev \
        libcurl4-openssl-dev \
        libgfortran4 \
        liblzma-dev \
        libpcre2-dev \
        libpng-dev \
        libreadline-dev \
        libssl-dev \
        libxml2-dev \
        libxt-dev libcairo2-dev \
        locales \
        tzdata \
        unzip \
        wget \
        zlib1g-dev && \
    locale-gen en_US.UTF-8 && \
    apt-get clean

# Install non-OS dependencies as their own layer
COPY conda.yml /tmp
RUN cd /tmp && \
    wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba && \
    mkdir -p /opt/micromamba/bin && \
    mv bin/micromamba /opt/micromamba/bin/ && \
    echo "y" | micromamba shell init -s bash -p /opt/micromamba && \
    . /opt/micromamba/etc/profile.d/micromamba.sh && \
    micromamba activate && \
    micromamba install -y python=3.9 mamba -c defaults -c conda-forge -c bioconda && \
    mamba env update --name base --file conda.yml --prune && \
    mamba clean -afy && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /opt/askcell
COPY ./ ./

RUN pip install -U pip && \
    pip install poetry && \
    poetry config virtualenvs.create false && \
    poetry install && \
    pip install setuptools
# poetry install removes setuptools, although it is set as dependency.
# To fix this, doing pip install setuptools after poetry install. Will have to fix this!

CMD ["pytest"]
