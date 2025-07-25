% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DonorAssignmentStandardAnalysis.R
\name{getSingletDoubletDF}
\alias{getSingletDoubletDF}
\title{Evaluate Donor Assignment}
\usage{
getSingletDoubletDF(
  likelihoodSummaryFile,
  doubletLikelihoodFile,
  expectedSamplesFile = NULL,
  doubletPvalue = 0.9,
  bestPairPvalue = 0.9,
  ignoreDiffuseDoublets = FALSE,
  fdrThreshold = 0.05,
  anonymizeDonors = FALSE
)
}
\arguments{
\item{likelihoodSummaryFile}{The dropulation single donor likelihood file.}

\item{doubletLikelihoodFile}{The dropulation doublet likelihood file.}

\item{expectedSamplesFile}{A file containing the expected list of donors in the experiment. The file has no header,
and a single donor ID on each line.}

\item{doubletPvalue}{The p-value threshold for calling a cell a doublet}

\item{bestPairPvalue}{The p-value threshold for calling a cell a diffuse doublet}

\item{ignoreDiffuseDoublets}{If TRUE, diffuse doublets are labeled}

\item{fdrThreshold}{The false discovery rate threshold for labeling singlets}

\item{anonymizeDonors}{If set to true, donor IDs are changed to DONOR_1, DONOR_2, etc.  Purely for sharing 
visualizations of not-yet-public data.}
}
\value{
A dataframe containing the cell barcodes, donor assignments, and classification of each cell barcode
}
\description{
Given the donor assignment and doublet results, construct a data frame of cell barcodes 
and their donor assignments, along with classification of each cell barcode
as a singlet, doublet, or diffuse doublet.
}
