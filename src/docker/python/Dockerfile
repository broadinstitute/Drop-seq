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

# After podman catches up on ARG scoping, remove the duplicated variables and
# only use the first introduction!
# https://github.com/containers/buildah/issues/5762

FROM ubuntu:24.04 AS common
ARG BASEDIR=/usr/local/dropseq

FROM common AS build
ARG BASEDIR=/usr/local/dropseq
ARG SRCVERSION=unknown
ARG SRCDIR=/src/Drop-seq
COPY . $SRCDIR/
WORKDIR $SRCDIR/
RUN $SRCDIR/src/docker/python/build/install.sh

FROM common AS runtime
ARG BASEDIR=/usr/local/dropseq
ENV PATH=$BASEDIR:$PATH
COPY --from=build $BASEDIR $BASEDIR
