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

# Builds DropSeq.utilities first, then all other packages. If the dependency tree gets more involved we may need to
# use an external file to specify the order of installation, or use something like remotes::local_package_deps()
# https://remotes.r-lib.org/reference/package_deps.html

# Overrides the BioConductor version as a workaround for https://github.com/r-lib/remotes/issues/798

set -euo pipefail

pkg_dir="$SRCDIR"/src/R/packages
source /etc/lsb-release
ppm_repo=https://packagemanager.posit.co/cran/__linux__/$DISTRIB_CODENAME/latest

mkdir -p "$BASERLIBS"

pkg_names="'DropSeq.utilities'"
for pkg_description in "$pkg_dir"/*/DESCRIPTION; do
    pkg_name=$(basename "$(dirname "$pkg_description")")
    if [[ "$pkg_name" == "DropSeq.utilities" ]]; then continue; fi
    pkg_names+=", '$pkg_name'"
done

Rscript -e '
options(
  HTTPUserAgent = sprintf(
    "R/%s R (%s)",
    getRversion(),
    paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])
  )
)

base_rlibs <- "'"$BASERLIBS"'"
ppm_repo <- "'"$ppm_repo"'"
pkg_dir <- "'"$pkg_dir"'"
pkg_names <- c('"$pkg_names"')

init_withr <- function() {
  message("Initializing withr")
  install.packages(pkgs = "withr", lib = base_rlibs, repos = ppm_repo)
  library(withr, lib.loc = base_rlibs)
}

init_remotes <- function() {
  message("Initializing remotes")
  install.packages(pkgs = "remotes", lib = base_rlibs, repos = ppm_repo)
  library(remotes, lib.loc = base_rlibs)
}

init_bioconductor <- function() {
  message("Initializing BioConductor")
  install.packages("BiocManager", lib = base_rlibs, repos = ppm_repo)
  library(BiocManager, lib.loc = base_rlibs)
  BiocManager::install(lib = base_rlibs)

  if (getRversion() >= "4.4.0") {
    bioc_manager_version <- BiocManager::version()
    remotes_package_version <- packageVersion("remotes", lib.loc = base_rlibs)
    message(sprintf(
      "NOTE: %s, remotes version %s, forcing R_BIOC_VERSION to %s to work around https://github.com/r-lib/remotes/issues/798",
      R.version.string, remotes_package_version, bioc_manager_version
    ))

    Sys.setenv(R_BIOC_VERSION = as.character(bioc_manager_version))
  }
}

install_dropseq_package <- function(pkg_name) {
  message(sprintf("Installing Drop-seq package: %s", pkg_name))
  withr::local_libpaths(base_rlibs)
  remotes::install_local(
    path = file.path(pkg_dir, pkg_name),
    upgrade = FALSE,
    repos = ppm_repo,
    lib = base_rlibs
  )
  message(sprintf("Verifying installation of Drop-seq package: %s", pkg_name))
  library(pkg_name, character.only = TRUE, lib.loc = base_rlibs)
}

init_withr()
init_remotes()
init_bioconductor()
lapply(pkg_names, install_dropseq_package)
'
