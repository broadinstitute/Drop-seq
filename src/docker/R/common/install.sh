#!/usr/bin/env bash
# MIT License
#
# Copyright 2025 Broad Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Installs the Drop-seq R common dependencies shared during build and runtime

# Based on / hat tip:
#   - https://cloud.r-project.org/bin/linux/ubuntu/
#   - https://github.com/rocker-org/rocker-versioned2/blob/77e4b1c/scripts/install_R_source.sh

set -euo pipefail

export DEBIAN_FRONTEND=noninteractive

# Install R public key.
apt-get -qq update
apt-get -qq install wget

source /etc/lsb-release
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc >> /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
echo "deb https://cloud.r-project.org/bin/linux/ubuntu ${DISTRIB_CODENAME}-cran40/" >> /etc/apt/sources.list

# Install R along with tools required to install the Drop-seq R_LIBS and their dependencies.
r_version=4.4.3
r_base_version=${r_version}-1.${DISTRIB_RELEASE/./}.0

apt-get -qq update
apt-get -qq install --no-install-recommends \
    automake \
    build-essential \
    cmake \
    file \
    gfortran \
    libblas-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
    libgeos-dev \
    libglpk40 \
    libhdf5-dev \
    libjpeg-dev \
    liblapack-dev \
    liblzma-dev \
    libopenblas-dev \
    libpng-dev \
    libpq-dev \
    libssl-dev \
    libtool \
    libxml2-dev \
    libz-dev \
    pkg-config \
    tk \
    zlib1g-dev \
    r-base-core="${r_base_version}"
rm -rf /var/lib/apt/lists/*
