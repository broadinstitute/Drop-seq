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

# Builds and installs the Drop-seq Java programs

# NOTE: As of May 2025, this crashes when building this image on Apple Silicon with the following error:
#   [Too many errors, abort]
#   Aborted (core dumped)
# Tried adding -XX:UseSVE but that doesn't seem to be supported in the ubuntu 24.04 openjdk-21.
# Likely related to the following issues:
# - https://bugs.openjdk.org/browse/JDK-8345296
# - https://forums.docker.com/t/image-builds-fail-on-new-macbook-despite-working-fine-on-prior-apple-silicon/145772/6

set -euo pipefail

export DEBIAN_FRONTEND=noninteractive

apt-get -qq update
apt-get -qq install unzip

"$SRCDIR"/gradlew -Pmanifest_version="$SRCVERSION" clean distZip
unzip "$SRCDIR"/dropseq/build/distributions/dropseq*.zip -d /tmp/unzip/
mv /tmp/unzip/dropseq* "$BASEDIR"
