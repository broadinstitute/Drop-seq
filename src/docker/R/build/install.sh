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

# Builds and installs the Drop-seq R packages

# Downloads binary R packages by using the Posit Package Manager and specifying the user agent
# - https://docs.posit.co/rspm/admin/serving-binaries.html#using-linux-binary-packages
# - https://docs.posit.co/rspm/admin/serving-binaries.html#binary-user-agents

# Installs into an empty R_LIBS directory, and explicitly sets upgrade=FALSE
# to avoid upgrading the DropSeq packages from GitHub.

set -euo pipefail

pkg_dir="$SRCDIR"/src/R/packages
source /etc/lsb-release
ppm_repo=https://packagemanager.posit.co/cran/__linux__/$DISTRIB_CODENAME/latest

mkdir -p "$BASERLIBS"

Rscript -e '
options(
  HTTPUserAgent = sprintf(
    "R/%s R (%s)",
    getRversion(),
    paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])
  )
)

install.packages(
  pkgs = c("withr", "remotes"),
  repos = "'"$ppm_repo"'"
)

withr::with_libpaths(
  new = "'"$BASERLIBS"'",
  code = remotes::install_local(
    path = c(
      "'"$pkg_dir"'/DropSeq.utilities",
      "'"$pkg_dir"'/DropSeq.dropulation"
    ),
    upgrade = FALSE,
    repos = "'"$ppm_repo"'"
  )
)
'
