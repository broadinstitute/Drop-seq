% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geneLeveleQTLPlots.R
\name{plotGeneQTLs}
\alias{plotGeneQTLs}
\title{Plot many eQTL expression vs genotype plots on one page.}
\usage{
plotGeneQTLs(
  snps_location_file_name,
  expression_file_name,
  SNP_file_name,
  eQTLPermutationResultFile,
  outPDF = NULL,
  outIndexSNPs = NULL,
  qValueThreshold = 0.05,
  geneList = NULL
)
}
\arguments{
\item{snps_location_file_name}{The location of each SNP.  Same number of lines as the SNP_file_name, encodes the snp name, chromosome, position.}

\item{expression_file_name}{A matrix of expression data, 1 row per gene, 1 column per sample.}

\item{SNP_file_name}{A matrix of genotypes, 1 row per SNP, 1 column per sample.  Encoded as 0,1,2 copies of the alternate allele.}

\item{eQTLPermutationResultFile}{The permuted eQTL result file}

\item{outPDF}{The file location to output the PDF of SNP/Expression plots.}

\item{outIndexSNPs}{The file location of the gene/SNPs with a q-value <= qValueThreshold}

\item{qValueThreshold}{Genes with a score <= this threshold will be emitted in the outIndexSNPs and outPDF.}

\item{geneList}{Override the qValueThreshold and plot these genes instead.}
}
\description{
Plot many eQTL expression vs genotype plots on one page.
}
