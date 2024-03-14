
# DropSeq.dropulation

<!-- badges: start -->
<!-- badges: end -->

Dropulation and Dropulation census functions.

## Installation

The preferred installation method is to go to the latest release in 
[Releases](https://github.com/broadinstitute/Drop-seq/releases) and follow the
 instructions there.

It is not recommended to install the latest version of DropSeq.dropulation from 
[GitHub](https://github.com/), because it is not necessarily stable, but if you choose to do so:

``` r
# install.packages("devtools")
devtools::install_github("broadinstitute/Drop-seq", subdir="src/R/packages/DropSeq.dropulation")
```

Note that dependent package DropSeq.utilities is Suggests: rather than Imports: in order
to make it less likely that a user will install a possibly-unstable version of that package
rather than a released version.