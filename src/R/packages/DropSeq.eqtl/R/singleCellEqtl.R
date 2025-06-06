
###############################################################################
# Model eQTL results at the single cell level.
# Run eQTL with and without interaction terms to one or more latent factors.
###############################################################################

#https://github.com/dynverse/anndata

#######################
# SeuratDisk install on OSX:
# if (!requireNamespace("remotes", quietly = TRUE)) {
# install.packages("remotes")
# }
# remotes::install_github("mojaveazure/seurat-disk")
########################

# library (Seurat)
# library (SeuratDisk)
# library (logger)
# library (lme4)
# library (ggpubr)
# library (ggplot2)
# library(cowplot)
# library (RColorBrewer)
# library (viridis)
# library (qvalue)
# library (qqman)
# library (muHVT)
# library (fitdistrplus)

#source("/Users/nemesh/dropseqrna/transcriptome/R/packages/DropSeq.eqtl/R/Matrix_eQTL_CommonFunctions.R")

# 20 genes that have high LRT pvalues for factor 2 from the factor 2 only test
#eQTLGeneList20=c("LINC01515","HSP90AA1","DPP10","TFRC","DPP6","NR4A3","IDI1","PRKG1","LINC01094","CACNB2","UAP1","EMP1","LRMDA","FOS","DGKB","HPSE2","PRRX1","FTH1","FSTL5","TPD52L1")

###############################
#OLD INPUTS
###############################
# donorCovariatesFile="/downloads/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.int.covariates_peer.txt"
# permutedPseudobulkEqtlResultsFile="/downloads/single_cell_eQTL/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.cis_qtl.txt.gz"

###############################
# REPLACEMENT INPUTS (Emi's PEER factors and eQTL discovery data)
###############################
# donorCovariatesFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10/maf_0.05_cisDist_1000kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.covariates_peer.txt"
# permutedPseudobulkEqtlResultsFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10/maf_0.05_cisDist_1000kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.eQTL_eigenMT_results.txt"


#################################
# INITIAL DISCOVERY PARAMETERS
# This is the "old" method before data compression was extracted into a separate process
# This runs discovery on all of the individual cells, and is quite slow.
#################################

# snpFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10/maf_0.05_cisDist_1000kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.genotype_matrix.txt.gz"
# h5SeuratFile="/broad/mccarroll/dropulation/analysis/latent_factor/BA46_n180_Astrocytes/raw_data/BA46_n180_astrocyte.h5seurat"
# cellMetaDataFile="/broad/mccarroll/dropulation/analysis/latent_factor/metadata/BA46/BA46.n191/BA46_cell_level_metadata.txt"
# donorCovariatesFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10/maf_0.05_cisDist_1000kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.covariates_peer.txt"
# permutedPseudobulkEqtlResultsFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10/maf_0.05_cisDist_1000kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.eQTL_eigenMT_results.txt"
# usageFile="/broad/mccarroll/dropulation/analysis/latent_factor/BA46_n180_Astrocytes/BA46_n180_Astrocytes/BA46_n180_Astrocytes.usages.k_11.dt_0_05.consensus.txt"
# usePeerFactors=F;
# factorName="FACTOR_1";
# #eQTLGeneList=NULL; chromosome="chr1";qvalueThreshold=0.05; nPerm=0;
# eQTLGeneList=c("RERE", "CUL3", "CERS5", "BAALC-AS2", "NPIPB12", "NMRAL1", "CLP1", "FAM120AOS"); chromosome=NULL;qvalueThreshold=NULL; nPerm=0;
# outDetailsFile="/downloads/single_cell_eQTL/8genes/BA46_Astrocyte.single_cell_eQTLs.baseline.factor1.details.txt"
# outSummaryFile="/downloads/single_cell_eQTL/8genes/BA46_Astrocyte.single_cell_eQTLs.baseline.factor1.summary.txt"

# singleCellEqtl(snpFile, h5SeuratFile,cellMetaDataFile, donorCovariatesFile, permutedPseudobulkEqtlResultsFile, usageFile, qvalueThreshold, usePeerFactors, eQTLGeneList, chromosome, factorName, nPerm, outDetailsFile, outSummaryFile)

###################################################################
# INITIAL DISCOVERY PARAMETERS - multiple factors simultaneously
# This is the "old" method before data compression was extracted into a separate process
# This runs discovery on all of the individual cells, and is quite slow.
####################################################################

# snpFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10/maf_0.05_cisDist_1000kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.genotype_matrix.txt.gz"
# h5SeuratFile="/broad/mccarroll/dropulation/analysis/latent_factor/BA46_n180_Astrocytes/raw_data/BA46_n180_astrocyte.h5seurat"
# cellMetaDataFile="/broad/mccarroll/dropulation/analysis/latent_factor/metadata/BA46/BA46.n191/BA46_cell_level_metadata.txt"
# donorCovariatesFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10/maf_0.05_cisDist_1000kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.covariates_peer.txt"
# permutedPseudobulkEqtlResultsFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10/maf_0.05_cisDist_1000kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.eQTL_eigenMT_results.txt"
# usageFile="/broad/mccarroll/dropulation/analysis/latent_factor/BA46_n180_Astrocytes/BA46_n180_Astrocytes/BA46_n180_Astrocytes.usages.k_11.dt_0_05.consensus.txt"
# usePeerFactors=F;
# factorName=NULL;
# # eQTLGeneList=NULL; chromosome="chr1";qvalueThreshold=0.05; nPerm=0;
# eQTLGeneList=eQTLGeneList20; chromosome=NULL;qvalueThreshold=NULL; nPerm=0;
# #eQTLGeneList=c("RERE", "CUL3", "CERS5", "BAALC-AS2", "NPIPB12", "NMRAL1", "CLP1", "FAM120AOS"); chromosome=NULL;qvalueThreshold=NULL; nPerm=0;
# outDetailsFile="/downloads/single_cell_eQTL/20genes/BA46_Astrocyte.single_cell_eQTLs.baseline_all_factors.details.txt"
# outSummaryFile="/downloads/single_cell_eQTL/20genes/BA46_Astrocyte.single_cell_eQTLs.baseline_all_factors.summary.txt"

# debug (parseCnmfFactor); debug (prepareDataForSingleCellEqtl)
# singleCellEqtl(snpFile, h5SeuratFile,cellMetaDataFile, donorCovariatesFile, permutedPseudobulkEqtlResultsFile, usageFile, qvalueThreshold, usePeerFactors, eQTLGeneList, chromosome, factorName, nPerm, outDetailsFile, outSummaryFile)

##################################
# DATA COMPRESSION PARAMETERS
##################################

# FACTOR 2 only

# snpFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10/maf_0.05_cisDist_1000kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.genotype_matrix.txt.gz"
# h5SeuratFile="/broad/mccarroll/dropulation/analysis/latent_factor/BA46_n180_Astrocytes/raw_data/BA46_n180_astrocyte.h5seurat"
# cellMetaDataFile="/broad/mccarroll/dropulation/analysis/latent_factor/metadata/BA46/BA46.n191/BA46_cell_level_metadata.txt"
# donorCovariatesFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10/maf_0.05_cisDist_1000kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.covariates_peer.txt"
# usageFile="/broad/mccarroll/dropulation/analysis/latent_factor/BA46_n180_Astrocytes/BA46_n180_Astrocytes/BA46_n180_Astrocytes.usages.k_11.dt_0_05.consensus.txt"
# usePeerFactors=F; useAllLatentFactors=F; factorName="FACTOR_2"
# numClustersPerDonor=100
# permutedPseudobulkEqtlResultsFile=NULL;
# qvalueThreshold=NULL; eQTLGeneList=NULL; chromosome=NULL;
# outCellFeaturesFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k100.factor2.cell_features.txt"
# outExpressionFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k100.factor2.h5seurat"
# outUsageFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k100.factor2.usages.txt"
# reduceDataPDF="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k100.factor2.compression_report.pdf"

#All latent factors

# snpFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10/maf_0.05_cisDist_1000kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.genotype_matrix.txt.gz"
# h5SeuratFile="/broad/mccarroll/dropulation/analysis/latent_factor/BA46_n180_Astrocytes/raw_data/BA46_n180_astrocyte.h5seurat"
# cellMetaDataFile="/broad/mccarroll/dropulation/analysis/latent_factor/metadata/BA46/BA46.n191/BA46_cell_level_metadata.txt"
# donorCovariatesFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10/maf_0.05_cisDist_1000kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.covariates_peer.txt"
# usageFile="/broad/mccarroll/dropulation/analysis/latent_factor/BA46_n180_Astrocytes/BA46_n180_Astrocytes/BA46_n180_Astrocytes.usages.k_11.dt_0_05.consensus.txt"
# usePeerFactors=F; useAllLatentFactors=T; factorName="FACTOR_2"
# numClustersPerDonor=50
# permutedPseudobulkEqtlResultsFile=NULL;
# qvalueThreshold=NULL; eQTLGeneList=NULL; chromosome=NULL; quantizationErrorThreshold=NULL
# outCellFeaturesFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k50.all_factor.cell_features.txt"
# outExpressionFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k50.all_factor.h5seurat"
# outUsageFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k50.all_factor.usages.txt"
# reduceDataPDF="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k50.all_factor.compression_report.pdf"

# reduceDataCmdLine(snpFile, h5SeuratFile, cellMetaDataFile, donorCovariatesFile, usageFile, factorName, useAllLatentFactors=useAllLatentFactors, usePeerFactors=usePeerFactors, numClustersPerDonor=numClustersPerDonor, permutedPseudobulkEqtlResultsFile=NULL, qvalueThreshold=NULL, eQTLGeneList=NULL, chromosome=NULL, reduceDataPDF, outCellFeaturesFile, outExpressionFile, outUsageFile)

#######################################
# GUIDED DISCOVERY WITH COMPRESSED DATA
# The compressed data is substituted in for the original
########################################

# Latent factor 2 data compression
# h5SeuratFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k100.factor2.h5seurat"
# usageFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k100.factor2.usages.txt"
# cellMetaDataFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k100.factor2.cell_features.txt"
# snpFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10/maf_0.05_cisDist_1000kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.genotype_matrix.txt.gz"
# donorCovariatesFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10/maf_0.05_cisDist_1000kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.covariates_peer.txt"
# permutedPseudobulkEqtlResultsFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10/maf_0.05_cisDist_1000kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.eQTL_eigenMT_results.txt"
# usePeerFactors=F
# factorName="FACTOR_2"; nPerm=0;
# eQTLGeneList=NULL; chromosome="chr1";qvalueThreshold=0.05; quantizationErrorThreshold=NULL
# eQTLGeneList=c("RERE", "CUL3", "CERS5", "BAALC-AS2", "NPIPB12", "NMRAL1", "CLP1", "FAM120AOS"); chromosome=NULL;qvalueThreshold=NULL; quantizationErrorThreshold=NULL
# outDetailsFile=NULL
# outSummaryFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k100.factor2.summary.txt"

# all latent factors data compression
# h5SeuratFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k50.all_factor.h5seurat"
# usageFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k50.all_factor.usages.txt"
# cellMetaDataFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k50.all_factor.cell_features.txt"
# snpFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10/maf_0.05_cisDist_1000kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.genotype_matrix.txt.gz"
# donorCovariatesFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10/maf_0.05_cisDist_1000kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.covariates_peer.txt"
# permutedPseudobulkEqtlResultsFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10/maf_0.05_cisDist_1000kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.eQTL_eigenMT_results.txt"
# usePeerFactors=F; nPerm=0;
# factorName=NULL;
# eQTLGeneList=NULL; chromosome="chr1";qvalueThreshold=0.05; quantizationErrorThreshold=NULL; nPerm=0;
# outDetailsFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k20.all_factor.details.txt"
# outSummaryFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k20.all_factor.summary.txt"


#all latent factors selelect QTL list with permutation
# h5SeuratFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k50.all_factor.h5seurat"
# usageFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k50.all_factor.usages.txt"
# cellMetaDataFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k50.all_factor.cell_features.txt"
# snpFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10/maf_0.05_cisDist_1000kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.genotype_matrix.txt.gz"
# donorCovariatesFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10/maf_0.05_cisDist_1000kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.covariates_peer.txt"
# permutedPseudobulkEqtlResultsFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10/maf_0.05_cisDist_1000kb/BA46.n191.raw.dge.astrocyte.All.noOutliers.sct.n.weightedMean.NUM_UMIS.peer10.maf_0.05_cisDist_1000kb.eQTL_eigenMT_results.txt"
# usePeerFactors=F; factorName=NULL;
# chromosome=NULL ;qvalueThreshold=NULL; quantizationErrorThreshold=NULL; nPerm=20;
# outDetailsFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k20.all_factor.details.txt"
# outSummaryFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k20.all_factor.summary.txt"
# eQTLGeneList=c("RERE", "CUL3", "CERS5", "BAALC-AS2", "NPIPB12", "NMRAL1", "CLP1", "FAM120AOS"); chromosome=NULL;qvalueThreshold=NULL

#singleCellEqtl(snpFile, h5SeuratFile,cellMetaDataFile, donorCovariatesFile, permutedPseudobulkEqtlResultsFile, usageFile, qvalueThreshold, usePeerFactors, eQTLGeneList, chromosome, factorName, useAllLatentFactors, quantizationErrorThreshold, nPerm=nPerm, outDetailsFile, outSummaryFile)



#####################################
# Guided discovery plotting
#####################################

# uMAPCoordinatesFile="/broad/mccarroll/dropulation/analysis/latent_factor/metadata/BA46/BA46.n191/uMAP/BA46.astrocyte.sct.scaled.centered.bbknn.donors.uMap_coordinates.txt"

# outDetailsFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.details.txt"
# outSummaryFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.summary.txt"
# outPermutedDetailsFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.permuted.details.txt"
# outPermutedSummaryFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.permuted.summary.txt"
# interactionEqtlSummaryFile="/broad/mccarroll/nemesh/single_cell_eQTL/BA46.astrocyte.noOutliers.single_cell_eQTL_interaction.summary.txt"
# outInteractionEqtlPDF="/broad/mccarroll/nemesh/single_cell_eQTL/BA46.astrocyte.noOutliers.single_cell_eQTL_interaction.report2.pdf"

# outDetailsFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k10.details.txt"
# outSummaryFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k10.summary.txt"
# reduceDataPDF="/broad/mccarroll/nemesh/single_cell_eQTL/BA46.astrocyte.noOutliers.single_cell_eQTL.k10.compression_report.pdf"


#################################
# TODO:
# 1) decouple data preparation/collapse from analysis [done]
# 2) Genome wide testing
# 3) Permutation of index SNP - maybe 1000 permutations, fit beta distribution, estimate pvalue?
#################################

#density plot of expression by factor, one data series per allele.
#maybe regress out the other data features, so we're plotting the residuals of expression per cell vs genotype?

#' Discover eQTL genotype interactions with latent factors
#'
#' For a given set of eQTLs, fit those SNP/gene pairs using a poisson regression with random effects for the donor
#' and batch.  The set of SNP/gene pairs tested is defined by the permutedPseudobulkEqtlResultsFile, the qvalue threshold to filter
#' that file, the chromosome to optionally further filter that gene list to a single chromosome, and optionally an eQTLGeneList
#' that contains an explicit set of genes to test.
#'
#' If usageFile is null, run the standard eQTL analysis without the latent factor interaction.
#' If the usage file is provided, runs the genotype:latent factor interaction analysis
#' @param snpFile A SNP genotype matrix as generated by the eQTL pipeline
#' @param h5SeuratFile A file containing single cell expression data in h5 Seruat format
#' @param cellMetaDataFile A file containing a matrix of cell level metadata
#' @param donorCovariatesFile A file containing donor known covariates
#' @param permutedPseudobulkEqtlResultsFile The result file from the eQTL pipeline containg the permuted eQTL results generated by a pseudobulk analysis.
#' @param qvalueThreshold Filter snp/gene interactions in the permutedPseudobulkEqtlResultsFile
#' @param usageFile A file containing a matrix of latent factor scores
#' @param usePeerFactors If true, include the donor covariates from that match "PEER_"
#' @param eQTLGeneList An explicit list of genes to further filter the set of SNP/Gene interactions being tested
#' @param chromosome A single chromosome to filter the SNP/Gene interactions on
#' @param factorName Selects the factor(s) from the usageFile to analyze.  If null, regress all latent factors
#' @param quantizationErrorThreshold If the cell level metadata contains the column QUANT_ERROR, the data has been previously compressed, and will be filtered to cells with values less than this score.
#' @param nPerm The number of permutations to run to generate a permuted pvalue.  This has two functions:
#' 1) It generates an permuted pvalue where the (number of permuted results > empiric result)+1 / total results +1 = permuted pvalue.
#' 2) The distributon of permuted pvalues is used to fit a beta distribution to compute a final permuted pvalue, in the
#' same manner as fastQTL.
#' @param outDetailsFile For each SNP/Gene pair, ouputs all the coefficents in the fitted main model,
#' as well as the pvalue from the comparison of the fitted to the null model. One coefficient per line.
#' @param outSummaryFile For each SNP/Gene pair, outputs the genotype, factor, and genotype:factor interaction
#' as well as the pvalue from the comparison of the fitted to the null model.  One result per line.
#' @export
#TODO: rename to directedSingleCellEqtlDiscovery
singleCellEqtl<-function (snpFile, h5SeuratFile,
                    cellMetaDataFile, donorCovariatesFile, permutedPseudobulkEqtlResultsFile,
                    usageFile=NULL, qvalueThreshold=0.05, usePeerFactors=T, eQTLGeneList=NULL, chromosome=NULL,
                    factorName=c("FACTOR_2"), nPerm=0, outDetailsFile=NULL, outSummaryFile=NULL) {

    d=prepareDataForSingleCellEqtl(snpFile, h5SeuratFile,
                                   cellMetaDataFile, donorCovariatesFile, permutedPseudobulkEqtlResultsFile,
                                   usageFile, qvalueThreshold=qvalueThreshold, eQTLGeneList=eQTLGeneList, chromosome=chromosome,
                                   factorList=factorName, quantizationErrorThreshold=NULL)

    metricsDF=d$metricsDF
    bulkResults=d$bulkResults
    genotypeDF=d$genotypeDF
    expData=d$expData
    usageDF=d$usage

    log_info("Number of cells in analysis [", dim (expData)[2], "]")
    if (is.null(usageFile)) {
        result=runManyGenotypeRegressions (metricsDF, bulkResults, genotypeDF, expData, usePeerFactors)
    } else {
        system.time(result<-runManyGenotypeSingleFactorRegressions (metricsDF, bulkResults, genotypeDF, expData, usePeerFactors, usageDF, factorName, nPerm))
    }

    summarizedResult=buildSummaryReportSingleFactor(result)
    if (!is.null(outSummaryFile))
        write.table(summarizedResult, outSummaryFile, row.names=F, col.names = T, quote=F, sep="\t")


    if (!is.null(outDetailsFile))
        write.table(result, outDetailsFile, row.names=F, col.names = T, quote=F, sep="\t")

}

#' Get a list of cell level features to use in data reduction and regression analysis
#'
#' @param fixedFeatures The known cell level metadata
#' @param factorName A latent factor name to use in downstream analysis (optional)
#' @param usePeerFactors Should PEER factors be included in analysis?
#' @param metricsDF The cell level metadata that contains fixed features and PEER factors
#' @param useAllLatentFactors If true and usageDF is not null, extract all the latent factor names to include in analysis. (optional)
#' @param usageDF A matrix of latent factor by cell scores (optional)
#'
#' @return A list of features from the metricsDF and possibly usageDF that will be used by downstream analysis
getFeatureList<-function (fixedFeatures=c("SCALED_LOG_UMIS", "pct_mt"), factorName=NULL, usePeerFactors=FALSE, metricsDF, useAllLatentFactors=F, usageDF=NULL) {

    #start with the given factor name.
    factorNames=factorName

    #if useAllLatentFactors is true and usageDF exists, add in all all latent factors.
    if (useAllLatentFactors && !is.null(usageDF))
        factorNames=sort(unique (c(factorName, colnames(usageDF))))

    featureList=c(fixedFeatures, factorNames)
    if (usePeerFactors) {
        peerFactors=colnames(metricsDF)[grep ("PEER", colnames (metricsDF))]
        if (length(peerFactors)>0) featureList=c(featureList, peerFactors)
    }
    return (featureList)
}

#' Permutes data to calculate an FDR of the non-permuted results
#'
#' If no usageFile is provided, test the effect of genotype on expression.  The permutation
#' test permutes the donor genotype, so all cells that have a particular donor genotype will
#' retain the same genotype, but that genotype may be different.
#'
#' If a usage file is provided, test the effect of the genotype:factor interaction.  The
#' permutation test permutes the factor score for each cell.
#' @param outPermutedDetailsFile For each SNP/Gene pair, ouputs all the coefficents in the fitted main model,
#' as well as the pvalue from the comparison of the fitted to the null model. One coefficient per line.
#' @param outPermutedSummaryFile For each SNP/Gene pair, outputs the genotype, factor, and genotype:factor interaction
#' as well as the pvalue from the comparison of the fitted to the null model.  One result per line.
#' @export
#' @inheritParams singleCellEqtl
singleCellEqtlPermuted<-function (snpFile, h5SeuratFile,
                          cellMetaDataFile, donorCovariatesFile, permutedPseudobulkEqtlResultsFile,
                          usageFile=NULL, qvalueThreshold=0.05, usePeerFactors=T, eQTLGeneList, chromosome=NULL,
                          factorName=c("FACTOR_2"), numClustersPerDonor=NULL,
                          outPermutedDetailsFile=NULL, outPermutedSummaryFile=NULL) {

    d=prepareDataForSingleCellEqtl(snpFile, h5SeuratFile,
                                   cellMetaDataFile, donorCovariatesFile, permutedPseudobulkEqtlResultsFile,
                                   usageFile, qvalueThreshold=qvalueThreshold, eQTLGeneList=eQTLGeneList, chromosome=chromosome,
                                   factorList=factorName)

    metricsDF=d$metricsDF
    bulkResults=d$bulkResults
    genotypeDF=d$genotypeDF
    expData=d$expData
    usageDF=d$usage

    if (!is.null(numClustersPerDonor)) {
        featureList=getFeatureList (fixedFeatures=c("SCALED_LOG_UMIS", "pct_mt"), factorName, usePeerFactors, metricsDF)
        z=reduceData(numClustersPerDonor=numClustersPerDonor, featureName=factorName, featureList, metricsDF, expData, usageDF, reduceDataPDF)
        metricsDF=z$metricsDF
        expData=z$expData
    }

    if (is.null(usageFile)) {
        #permute the donor genotypes - here we permute the column labels to swap genotypes.
        colnames(genotypeDF)[2:dim(genotypeDF)[2]]=sample (colnames (genotypeDF)[-1])
        result=runManyGenotypeRegressions (metricsDF, bulkResults, genotypeDF, expData, usePeerFactors)
    } else {
        #Permute the cell state score - here we permute the cell barcode labels.
        #It's assumed that the usageDF has the same cell barcode label order as the metricsDF,
        #so permute the cell barcode names and then resort.
        rownames(usageDF)=sample (rownames(usageDF))
        usageDF=usageDF[match(metricsDF$CELL_BARCODE_FINAL, rownames(usageDF)),]
        result=runManyGenotypeSingleFactorRegressions (metricsDF, bulkResults, genotypeDF, expData, usePeerFactors, usageDF, factorName)
    }

    summarizedResult=buildSummaryReportSingleFactor(factorName, result)

    if (!is.null(outPermutedDetailsFile))
        write.table(result, outPermutedDetailsFile, row.names=F, col.names = T, quote=F, sep="\t")

    if (!is.null(outPermutedSummaryFile))
        write.table(summarizedResult, outPermutedSummaryFile, row.names=F, col.names = T, quote=F, sep="\t")

}

baselineStdAnalysis<-function (baselineEqtlSummaryFile, outBaselineEqtlPDF=NULL, snpFile, h5SeuratFile,
                               cellMetaDataFile, donorCovariatesFile, permutedPseudobulkEqtlResultsFile,
                               usageFile, usePeerFactors=T, eQTLGeneList, factorName=c("FACTOR_2")) {

    d=prepareDataForSingleCellEqtl(snpFile, h5SeuratFile,
                                   cellMetaDataFile, donorCovariatesFile, permutedPseudobulkEqtlResultsFile,
                                   usageFile, qvalueThreshold=1, eQTLGeneList=NULL, chromosome=NULL, factorList=factorName)

    metricsDF=d$metricsDF
    bulkResults=d$bulkResults
    genotypeDF=d$genotypeDF
    expData=d$expData
    usageDF=d$usage

    r=read.table(baselineEqtlSummaryFile, header=T, stringsAsFactors = F, sep="\t")
    r=r[match(bulkResults$gene, r$gene),]
    if (!all (r$gene==bulkResults$gene))
        stop ("genes don't match in the bulk and single cell results")

    df=data.frame(gene=bulkResults$gene, snp=bulkResults$SNP, bulkPval=bulkResults$`p-value`, bulkQval=bulkResults$qvalue,
                  singleCellPval=r$lrt_pval, bulkEffectSize=bulkResults$beta, singleCellEffectSize=r$genotype_estimate,
                  bulkZscore=bulkResults$beta/bulkResults$beta_se, singleCellZscore=r$genotype_estimate/r$genotype_std_error)

    if (!is.null(outBaselineEqtlPDF)) pdf(outBaselineEqtlPDF)


    plotBaselineEqtlSignTest(df, 1, "All Genes \n ", ylim=c(-2,2))
    plotBaselineEqtlSignTest(df, 0.05, "Genes significant in bulk [qvalue<=0.05] \n ")

    plotBaseLineZscore (df, qvalueThreshold=1, titleStrPrefix="All Genes")
    plotBaseLineZscore (df, qvalueThreshold=0.05, titleStrPrefix="Genes significant in bulk [qvalue<=0.05] \n ")

    plotBaselineEqtlPvalues(df, 1, "All genes (qvalue <=1)\n")
    plotBaselineEqtlPvalues(df, 0.05, "Significant genes (qvalue <0.05)\n")
    plotBaselineEqtlPvalues(df[df$bulkPval>1e-15,], 0.05, "Less significant genes (qvalue <0.05, pvalue > 1e-15)\n")
    #plotBaselineEqtlPvalues(df[df$bulkPval>1e-10,], 0.05, "Less significant genes (qvalue <0.05, pvalue > 1e-10)\n")

    if (!is.null(outBaselineEqtlPDF)) dev.off()

}

#' Generate eQTL interaction report
#'
#' Generate a qq plot of eGene : factor interactions, followed by a plot of each eGene eQTL effect.
#'
#' @param interactionEqtlSummaryFile The result of running singleCellEqtl
#' @param outInteractionEqtlPDF The location of the output PDF.
#' @param compressedImages Use some semi-fancy ggplot stuff to compress each image as a PNG then write it to the PDF.
#' Causes no end of misery on UNIX, but works on OSX?
#' @inheritParams singleCellEqtl
#' @export
interactionStdAnalysis<-function (interactionEqtlSummaryFile, outInteractionEqtlPDF,snpFile, h5SeuratFile,
                                  cellMetaDataFile, donorCovariatesFile, permutedPseudobulkEqtlResultsFile,
                                  usageFile, uMAPCoordinatesFile, qvalueThreshold=0.05, factorName=c("FACTOR_2"),
                                  compressedImages=FALSE) {

    d=prepareDataForSingleCellEqtl(snpFile, h5SeuratFile,
                                   cellMetaDataFile, donorCovariatesFile, permutedPseudobulkEqtlResultsFile,
                                   usageFile, qvalueThreshold=qvalueThreshold, factorList=factorName)

    metricsDF=d$metricsDF
    bulkResults=d$bulkResults
    genotypeDF=d$genotypeDF
    expData=d$expData
    usageDF=d$usage

    summarizedResult=read.table(interactionEqtlSummaryFile, header=T, stringsAsFactors = F, sep="\t")
    bulkResults=bulkResults[match(summarizedResult$gene, bulkResults$gene),]
    if (!all (summarizedResult$gene==bulkResults$gene))
        stop ("genes don't match in the bulk and single cell results")

    #correct for the number of genes tested.
    summarizedResult$qvalue=qvalue(summarizedResult$lrt_pval)$qvalues

    #load up the uMAP coordinates
    uMapCoordinates=parseUmapCoordinates(uMAPCoordinatesFile, cellBarcodeList=metricsDF$CELL_BARCODE_FINAL)

    #results to plot
    summarizedResultFiltered=summarizedResult[summarizedResult$qvalue<qvalueThreshold,]
    #sort the results by qvalue.
    summarizedResultFiltered=summarizedResultFiltered[order(summarizedResultFiltered$qvalue),]

    if (!is.null(outInteractionEqtlPDF)) pdf(outInteractionEqtlPDF)
    plotEmpiricPvalueDistribution(summarizedResult)

    for (i in 1:dim(summarizedResultFiltered)[1]) {
        r=summarizedResultFiltered[i,]
        log_info("Plottting interaction eGene [", r$gene,"] [", i, "/", dim(summarizedResultFiltered)[1], "]")
        plotEgeneSummary (r$gene, r$snp, factorName, expData, metricsDF, usageDF, genotypeDF,
                          method="avgLog2PlusOne", uMapCoordinates, summarizedResult, compressedImages)
    }

    if (!is.null(outInteractionEqtlPDF)) dev.off()


}

# Construct the model forumulas used in regressions
# If the factorList is NULL, test the effect of the genotype
# If the factorList is not null, test the interaction of the genotype:factor for each factor.
getForumlas<-function (df, usePeerFactors=T, factorList=NULL, useRandomFactors=T) {
    #Construct a model string that is the same for both tests - this is the fixed and random effects
    #that are in both models.
    jointModelStr=NULL

    if (useRandomFactors) {
        jointModelStr="(1 | DONOR)"
        #add the prefix random effect
        if ("PREFIX" %in% colnames(df) && useRandomFactors)
            jointModelStr=paste(jointModelStr, "(1 | PREFIX)", sep=" + ")
    }

    #add fixed effects
    if (length(jointModelStr)==0) sep="" else sep=" + "
    jointModelStr= paste(jointModelStr, "SCALED_LOG_UMIS + pct_mt", sep=sep)

    if (usePeerFactors) {
        peerFactors=colnames(df)[grep ("PEER", colnames (df))]
        log_info("Using [", length(peerFactors), "] PEER factors")
        jointModelStr=jointModelStr(fullModelStr, paste(peerFactors, collapse=" + "), sep=" + ")
    }

    #if there are no factors, test the genotype.
    if (is.null(factorList)) {
        #the null model is the full model, but without genotype.
        nullModelStr=paste("expression ~ ", jointModelStr, sep=" ")
        nullModelFormula=as.formula(nullModelStr)
        #full model includes the genotype
        fullModelStr=paste("expression ~ genotype", jointModelStr,  sep=" + ")
        fullModelFormula=as.formula(fullModelStr)
    } else {

        jointModelStr=paste("expression ~ genotype", jointModelStr , sep=" + ")
        #add factors
        jointModelStr=paste(jointModelStr, paste0(factorList, collapse = " + "), sep=" + ")
        #the null model is the genotype and the factor, but no interactions
        nullModelStr=jointModelStr
        nullModelFormula=as.formula(nullModelStr)

        #the full model includes the interactions
        #add in cNMF latent factor(s) being tested by the model
        latentFactorStr=paste0("genotype:", factorList, collapse = " + ")

        fullModelStr=paste(jointModelStr, latentFactorStr, sep=" + ")
        fullModelFormula=as.formula(fullModelStr)
    }

    r=list(fullModelFormula=fullModelFormula, nullModelFormula=nullModelFormula)
    return (r)

}

#run single cell eQTL regressions.
#This looks only at the effect of genotype on the counts, after controlling for other factors.
#the independent variable tested is genotype.
runManyGenotypeRegressions<-function (metricsDF, bulkResults, genotypeDF, expData, usePeerFactors=T) {

    #set up a dataframe with slots for all the analysis.
    #this may be faster to run many times if we use data table to change column values per snp/gene tested.
    df=data.table::copy(metricsDF)
    df[,genotype:=NA]
    df[,expression:=NA]

    idxToGenotype=match(df$DONOR, colnames(genotypeDF))

    #set up the formulas
    formulas=getForumlas (df, usePeerFactors, factorList=NULL)
    nullModelFormula=formulas$nullModelFormula
    fullModelFormula=formulas$fullModelFormula

    runOne<-function (index, bulkResults) {
        gene=bulkResults[index,]$gene
        snp=bulkResults[index,]$SNP
        log_info("Genotype only poisson regression gene [", gene, "] SNP [", snp, "]")
        runOneRegression (gene, snp, df, expData, genotypeDF, idxToGenotype, fullModelFormula, nullModelFormula)
    }

    log_info("Null Model: ", deparse1(nullModelFormula))
    log_info("Full Model: ", deparse1(fullModelFormula))

    r=do.call(rbind, lapply(1:dim(bulkResults)[1], runOne, bulkResults))
    r=cbind(regression_type="genotype", r)
    return (r)
}

testFixedModel<-function (metricsDF, bulkResults, genotypeDF, expData, usePeerFactors=T, usageDF, factorList=NULL) {
    factorList=colnames(usageDF)

    df=data.table::copy(metricsDF)
    idxToGenotype=match(df$DONOR, colnames(genotypeDF))

    gene=bulkResults[1,]$gene
    snp=bulkResults[1,]$SNP

    df[,genotype:=as.numeric (genotypeDF[genotypeDF$id==snp,idxToGenotype, with=F])]
    df[,expression:=as.numeric (expData[gene,])]

    #build model matrix and check for colinearity?
    fixedFormula=getForumlas (df, usePeerFactors, factorList, FALSE)
    formula=fixedFormula$nullModelFormula

    mm=model.matrix(formula, df)
    r=caret::findLinearCombos(mm)

    validateModel(formula, df)
    #should weights sum to the total number of observations to keep likelihoods similar to unweighted model?

    #df$weight=(df$NUM_CELLS/sum(df$NUM_CELLS))*dim(df)[1]
    #df$weight=(df$NUM_CELLS/sum(df$NUM_CELLS))
    #browser()



}

runManyGenotypeSingleFactorRegressions<-function (metricsDF, bulkResults, genotypeDF, expData, usePeerFactors=T, usageDF, factorList=NULL, nPerm=0) {

    #if no factor is selected use all factors.
    if (is.null(factorList)) {
        #factorList=colnames(usageDF)[-dim(usageDF)[2]]
        factorList=colnames(usageDF)
    }

    #set up a dataframe with slots for all the analysis.
    #this may be faster to run many times if we use data table to change column values per snp/gene tested.
    df=data.table::copy(metricsDF)
    df[,genotype:=NA]
    df[,expression:=NA]

    idxToGenotype=match(df$DONOR, colnames(genotypeDF))

    #set up the formulas
    formulas=getForumlas (df, usePeerFactors, factorList)
    nullModelFormula=formulas$nullModelFormula
    fullModelFormula=formulas$fullModelFormula

    fixedFormula=getForumlas (df, usePeerFactors, factorList, FALSE)

    runOne<-function (index, bulkResults) {
        gene=bulkResults[index,]$gene
        snp=bulkResults[index,]$SNP
        log_info("Genotype factor interaction poisson regression gene [", gene, "] SNP [", snp, "]")
        runOneRegression (gene, snp, df, expData, genotypeDF, idxToGenotype, fullModelFormula, nullModelFormula, factorList, nPerm)
    }

    log_info("Null Model: ", deparse1(nullModelFormula))
    log_info("Full Model: ", deparse1(fullModelFormula))

    r=do.call(rbind, lapply(1:dim(bulkResults)[1], runOne, bulkResults))
    r=cbind(regression_type="genotype_factor_interaction", r)
    return (r)

}

#snp="", genotypeDF
#gene="CERS5";snp="chr12:50071572:A:G";
#gene="RERE";snp="chr1:8388317:C:CTT";
#gene="RERE";snp="chr1:8358584:C:A"
runOneRegression<-function (gene, snp, df, expData, genotypeDF, idxToGenotype, fullModelFormula, nullModelFormula, factorList=NULL, nPerm=0) {

    df[,genotype:=as.numeric (genotypeDF[genotypeDF$id==snp,idxToGenotype, with=F])]
    df[,expression:=as.numeric (expData[gene,])]

    #TODO: REMOVE THIS!
    # if (TRUE) {
    #     outFile=paste("/downloads/single_cell_eQTL/20genes/", gene, "_data_", paste(factorList, collapse=":"), ".txt", sep="")
    #     write.table(df, outFile, row.names=F, col.names = T, quote=F, sep="\t")
    # }

    #build model matrix and check for colinearity?
    #model.matrix(forumla, data)

    #should weights sum to the total number of observations to keep likelihoods similar to unweighted model?

    #df$weight=(df$NUM_CELLS/sum(df$NUM_CELLS))*dim(df)[1]
    #df$weight=(df$NUM_CELLS/sum(df$NUM_CELLS))
    #browser()

    full_model <- lme4::glmer(formula=fullModelFormula, data=df, control = lme4::glmerControl(optimizer = "nloptwrap"), family = "poisson", nAGQ=0)
    null_model <- lme4::glmer(formula=nullModelFormula, data=df, control = lme4::glmerControl(optimizer = "nloptwrap"), family = "poisson", nAGQ=0)

    #experimental weighted models.
    #full_model <- lme4::glmer(formula=fullModelFormula, data=df, control = lme4::glmerControl(optimizer = "nloptwrap"), family = "poisson", nAGQ=0, weights=weight)
    #null_model <- lme4::glmer(formula=nullModelFormula, data=df, control = lme4::glmerControl(optimizer = "nloptwrap"), family = "poisson", nAGQ=0, weights=weight)

    model_lrt <- anova(null_model, full_model)

    out <- summary (full_model)$coefficients
    colnames(out) <- c("Estimate","StdError","zvalue","pval")

    out=cbind (Feature=rownames(out), data.frame(out, stringsAsFactors = F))
    rownames (out)=NULL

    #bind the full model as a feature.
    full_model_result=data.frame(Feature="FULL_MODEL_LRT", Estimate=NA, StdError=NA, zvalue=NA, pval=model_lrt$`Pr(>Chisq)`[2])
    out=rbind(out, full_model_result)

    #out=cbind(gene=gene, snp=snp, lrt_pval=model_lrt$`Pr(>Chisq)`[2], out)
    out=cbind(gene=gene, snp=snp,out)
    #permute if requested.
    if (nPerm>0) {
        out=permuteFactorScoresOneRegression(null_model, fullModelFormula, df, out, factorList, nPerm)
    }

    return (out)
}

#given the empiric full model results, the data and the model forumula, permute the factor scores and retest the model.
#Important: the relationship of scores within a cell must be maintained, so you want to shuffle all factor values for a single cell to a different random cell
#Rerun the full model many times to generate a distribution of p-values for each coefficient of interest.
#out is the results from the empiric permutation.
permuteFactorScoresOneRegression<-function (null_model, fullModelFormula, df, out, factorList, nPerm=0) {

    #to get more precise pvalues, need to test with the tail false when Z > 0 and tail TRUE when Z<0.
    # getPFromZ<-function (z) {
    #     2*pnorm(q=z, lower.tail=z<0)
    # }

    onePermutationPvalue<-function () {
        data.table::set(df, j = factorList, value = df[sample(1:dim(df)[1]),..factorList])
        permFull <- lme4::glmer(formula=fullModelFormula, data=df, control = lme4::glmerControl(optimizer = "nloptwrap"), family = "poisson", nAGQ=0)
        permNull <- lme4::glmer(formula=nullModelFormula, data=df, control = lme4::glmerControl(optimizer = "nloptwrap"), family = "poisson", nAGQ=0)
        model_lrt <- anova(permNull, permFull)
        full_model_p=model_lrt$`Pr(>Chisq)`[2]
        r=summary (p)$coefficients
        #r[,4]=sapply(r[,3], getPFromZ)
        c(r[,4], FULL_MODEL_LRT=full_model_p)
    }

    log_info(paste("Running [", nPerm, "] permutations"))

    r=replicate(nPerm, onePermutationPvalue(), simplify = FALSE)
    r=do.call(rbind, r)

    getOne<-function (index, out, r) {
        log_info("Index [", index, "]")
        empiricP=out[index,]$pval
        permutedPDist=r[,index]
        betaApproximationPvalue(empiricP, permutedPDist)
    }

    out$pval_perm=NA
    out$pval_beta=NA

    #which values should be calculated?
    idxFeatures=c(grep(":", out$Feature),which(out$Feature=="FULL_MODEL_LRT"))
    z=lapply(idxFeatures, getOne, out, r)
    z=do.call(rbind, z)
    out[idxFeatures,c("pval_perm", "pval_beta")]=z
    return (out)
}

#factorResult=read.table("/broad/mccarroll/nemesh/single_cell_eQTL/BA46.astrocyte.noOutliers.single_cell_eQTL.details.txt", header=T, stringsAsFactors=F, sep="\t")
#factorResult=read.table("/broad/mccarroll/nemesh/single_cell_eQTL/BA46.astrocyte.noOutliers.single_cell_eQTL_interaction.details.txt", header=T, stringsAsFactors=F, sep="\t")
buildSummaryReportSingleFactor<-function (factorResult) {

    interaction=unique (factorResult[grep ("genotype:", factorResult$Feature),]$Feature)
    featuresToKeep=c(interaction, "FULL_MODEL_LRT")
    colsToKeep=c("gene", "snp", "Feature", "Estimate", "zvalue", "pval")
    result=factorResult[factorResult$Feature %in% featuresToKeep, colsToKeep]
    return (result)
}

#BTotal is the genotype effect + the interaction score * latent factor score for a cell
#If no gene names are specified, return a matrix of cells (rows) by genes (columns) of
#this score.
#if a gene name is specified, only compute and return the matrix for the relevant data
getBTotal<-function (geneNames=NULL, factorName, summarizedResult, usageDF) {
    getOne<-function (geneName, factorName, summarizedResult, usageDF) {
        ge=summarizedResult[summarizedResult$gene==geneName,]$genotype_estimate
        ie=summarizedResult[summarizedResult$gene==geneName,]$interaction_estimate
        bTotal=ge+(usageDF[,factorName]*ie)
        return (bTotal)
    }

    #if the geneNames aren't set, use all gene names
    if (is.null(geneNames)) {
        geneNames=summarizedResult$gene
    }

    r=lapply(geneNames, getOne, factorName, summarizedResult, usageDF)
    rr=do.call(cbind, r)
    colnames (rr)=geneNames
    rr=data.frame(rr, stringsAsFactors = F, row.names = rownames(usageDF), check.names = F)
    return (rr)
}

#formula=fullFixedEffectsModel;data=df
#check to see if the fixed independent variables in the model are sufficiently non correlated
validateModel<-function (formula, df) {

    f=as.formula(formula)
    design <- model.matrix(f,data=df)

    if (qr(design)$rank!=ncol(design))
        stop ("design matrix may have redundant features")

    if (!limmaFullRank(design))
        stop ("Is full rank test failed!")

}

#factorName="FACTOR_2"; gene="RERE"
#order cells by latent factor score, compute pseudobulk expression for bins of the factor score.
getBinnedExpression<-function (gene, factorName, expData, metricsDF, usageDF, nBins=3,
                                         method=c("pseudobulk", "avgLog2PlusOne")) {
    #get the ordering of cells.
    i=order(usageDF[[factorName]])
    chunks=chunk2(i, nBins)
    r=lapply(chunks, getSummarizedExpression, gene, expData, metricsDF, method)
    return (r)
}

#compute the Btotal score for a gene, order the cells by the factor score and partition into bins
#and compute the average Btotal of the bins.
#returns a vector of length nBins with the average Btotal scores.
getBinnedBtotal<-function (gene, factorName, summarizedResult, usageDF, nBins=3) {
    i=order(usageDF[[factorName]])
    chunks=chunk2(i, nBins)
    btotal=getBTotal (geneNames=gene, factorName, summarizedResult, usageDF)
    getOne<-function (idx) mean (btotal[idx,1])
    sapply(chunks, getOne)
}

getSummarizedExpression<-function (idxCells=NULL, gene, expData, metricsDF, method=c("pseudobulk", "avgLog2PlusOne")) {
    if (method=="pseudobulk")
        return(createPseudoBulkExpression(idxCells, gene, expData, metricsDF))
    if (method=="avgLog2PlusOne")
        return(getAverageLog2Expression(idxCells, gene, expData, metricsDF))
    stop ("summarization method [", method, "] not recognized")
}

#Each point represents the average log2(UMI counts + 1) across all cells in the indicated
#subset of cells in a donor (n = 259), grouped by genotype.
getAverageLog2Expression<-function (idxCells=NULL, gene, expData, metricsDF) {

    if (is.null(idxCells))
        idxCells=1:dim(metricsDF)[1]

    #average log2(UMI counts + 1)
    avgLog2PlusOne<-function (x)
        mean (log2(x+1))

    x=tapply(expData[gene,idxCells], INDEX=metricsDF[idxCells,]$DONOR, FUN=avgLog2PlusOne)
    return (x)

}


#idx=1:50000
createPseudoBulkExpression<-function(idxCells=NULL, gene, expData, metricsDF) {
    if (is.null(idxCells))
        idxCells=1:dim(metricsDF)[1]

    x=tapply(expData[gene,idxCells], INDEX=metricsDF[idxCells,]$DONOR, FUN=sum)
    libTotal=tapply(metricsDF[idxCells,]$num_retained_transcripts, INDEX=metricsDF[idxCells,]$DONOR, FUN=sum)
    r=x/libTotal*100000
    return (r)
}


fixedEffectsRegression<-function (fullFixedEffectsModel, data) {

    summary(m1 <- glm(fullFixedEffectsModel, family="poisson", data=data))


}

limmaFullRank<-function (design) {
    x<- as.matrix(design)
    e <- eigen(crossprod(x), symmetric = TRUE, only.values = TRUE)$values
    e[1] > 0 && abs(e[length(e)]/e[1]) > 1e-13

}

scalePeerFactors<-function (df) {
    peerFactors=colnames(df)[grep ("PEER", colnames (df))]
    df[, (peerFactors) := lapply(.SD, scale), .SDcols = peerFactors]
    return (df)
}


###################################
# PLOTTING
###################################

plotEgeneSummary<-function (gene, snp, factorName, expData, metricsDF, usageDF, genotypeDF,
                            method=c("pseudobulk", "avgLog2PlusOne"), uMapCoordinates, summarizedResult=NULL, compressedImages=T) {

    #TODO: make the uMAP plot a PNG, not the entire plot.
    #that will reduce storage space and leave the PDF searchable.

    base_size=11
    legendKeySize=0.5
    if (compressedImages) {
        base_size=32
        legendKeySize=1
    }

    #uMapPlot=plotUmapFactor (factorName, usageDF, uMapCoordinates)
    uMapPlot=plotUmapBtotal (gene, snp, factorName, summarizedResult, usageDF, uMapCoordinates, titleSize=0, legendKeySize=legendKeySize, legendTitleStr="BTotal", base_size=base_size)

    eQTLPlots=plotEqtlBinnedByFactorScore (gene, snp, factorName, expData, metricsDF, usageDF, genotypeDF, nBins=3,
                                           method=method, summarizedResult, base_size)

    emptyPlot=ggplot() + theme_void()
    bottomPlots=cowplot::plot_grid(emptyPlot, uMapPlot, emptyPlot, nrow=1, ncol=3, rel_widths=c(0.25, 0.5, 0.25))
    result=cowplot::plot_grid(eQTLPlots, bottomPlots, nrow = 2, ncol=1, rel_heights=c(0.6, 0.4))

    if (compressedImages) {
        tempFile=paste(tempfile(), ".png", sep="")
        #ggplot2::ggsave(tempFile, bg = "white", width = 7, height = 7)
        png(tempFile, width=1600, height=1600, pointsize=10)
        print (result)
        dev.off()
        plotPNG(tempFile, xleft=0, ybottom=0, xright=1, ytop=1)
        file.remove(tempFile)
    } else {
        print (result)
    }

}

#gene="RERE";snp="chr1:8388317:C:CTT";
#if summarizedResult is not null, it is used to calculate the Btotal scores for the bins.
plotEqtlBinnedByFactorScore<-function (gene, snp, factorName, expData, metricsDF, usageDF, genotypeDF, nBins=3,
                                       method=c("pseudobulk", "avgLog2PlusOne"), summarizedResult=NULL, base_size=11) {
    r=getExpressionBinnedByFactorScore(gene, snp, factorName, expData, metricsDF, usageDF, genotypeDF, nBins, method, summarizedResult)
    data=r$data
    bTotal=r$summary$bTotal
    qvalue=NULL
    if (!is.null(summarizedResult)) {
        qvalue=summarizedResult[summarizedResult$gene==gene & summarizedResult$snp==snp,]$qvalue
    }
    r=plotManyEqtls(data, snp, gene, bTotal, qvalue, base_size)
    return (r)
}



#' Plot a group of eQTLs
#'
#' @param data Data as generated by getExpressionBinnedByFactorScore
#' @param gene The name of the gene being plotted
#' @param snp The name of the SNP being plotted
#' @param bTotal The BTotal score per cell as generated by getExpressionBinnedByFactorScore
#' @param qvalue Qvalue of the eQTL (optional)
#' @param base_size The ggplot2 base size for the plot
#'
#' @return A plot_grid result
#' @import cowplot ggplot2 viridis
plotManyEqtls<-function (data, gene, snp, bTotal, qvalue=NULL, base_size=11) {

    cell.types = names(data)
    max.value=max(sapply(data, function (x) max(x$expr)))
    max.value=max.value*1.1
    plots = lapply(cell.types, function(cell.type) {
        plotSingleEqtl(data[[cell.type]], cell.type, max.value, bTotal[[cell.type]], base_size)
    })
    strTitle=paste0(snp, ' / ', gene)
    if (!is.null(qvalue)) strTitle=paste(strTitle, "qvalue [", qvalue, "]")
    r=cowplot::plot_grid(
        #ggdraw() + draw_label(paste0(snp, ' / ', gene), size = 14, fontfamily = "Courier"),
        ggdraw() + draw_label(paste0(snp, ' / ', gene), size = 12*(base_size/11)),
        cowplot::plot_grid(plotlist = plots, nrow = 1),
        ncol=1, rel_heights=c(0.1, 1)
    )
    return (r)
}

#From Seva's code: /humgen/cnp04/sandbox/skashin/projects/CHDI/scripts/plot_eQTLs.R
#cell.type=cell.types[1];df=data[[cell.type]];bTotalBin=bTotal[[cell.type]]
plotSingleEqtl<-function (df, strTitle="", max.value, bTotalBin=NULL, base_size=11) {
    s=summary(lm(df$expr~ df$genotypes))
    intercept=round (coef(s)[1,1],2)
    beta=round (coef(s)[2,1],3)

    scalingFactor=base_size/11

    p=ggplot(df, aes(x = ALLELE, y = expr)) +
        #geom_point(size=1.5*(base_size/11)) +
        geom_violin(color = NA, fill = '#45A5F5') +
        geom_jitter(width = 0.2, size=1.5*scalingFactor) +
        theme_bw(base_size=base_size) +
        xlab('Allele') +
        ylab('Normalized expression') +
        ggtitle(strTitle) +
        coord_cartesian(ylim = c(0, max.value), xlim=c(1,3)) +
        theme(plot.title = element_text(hjust = 0.5, size = 10*scalingFactor)) +
        geom_smooth(aes(group = 1), color = 'red', formula = 'y ~ x', method = 'lm', se = FALSE) +
        #stat_cor(aes(group = 1, label = paste(..r.label.., ..p.label.., sep = "~`,`~")), label.y = max.value - 0.3, size = 3) +
        geom_text(x = 2, y = max.value*1.02, label = paste("int [", intercept, "] beta [", beta, "]", sep=""), parse = FALSE, size=2.5*scalingFactor)

    if (!is.null(bTotalBin)) {
        bTotalBin=round(bTotalBin,3)
        #p= p + geom_text(x = 2, y = max.value*0.98, label = paste("\u03B2total [",bTotalBin,"]", sep=""), parse = F, size=3)
        p= p + geom_text(x = 2, y = max.value*0.98, label = paste("Btotal [",bTotalBin,"]", sep=""), parse = F, size=2.5*scalingFactor)
    }
    return (p)
}



getExpressionBinnedByFactorScore<-function (gene, snp, factorName, expData, metricsDF, usageDF, genotypeDF,
                                            nBins=3, method=c("pseudobulk", "avgLog2PlusOne"), summarizedResult=NULL) {

    allCells=getSummarizedExpression (idxCells=NULL, gene, expData, metricsDF, method)
    b=getBinnedExpression (gene, factorName, expData, metricsDF, usageDF, nBins=3, method)
    g=genotypeDF[genotypeDF$id==snp,]

    #need to build a DF for each bin
    #expression=b[[1]]
    buildData<-function (expression) {
        ref=strsplit (snp, ":")[[1]][3]
        alt=strsplit (snp, ":")[[1]][4]
        alleleList=c(paste(ref, ref, sep="/"),paste(ref, alt, sep="/"), paste(alt, alt, sep="/"))
        donors=names(expression)
        genotypes=unlist(g[,donors,with=F])
        alleles=factor(alleleList[genotypes+1], levels=alleleList)
        df=data.frame(donor=donors, genotypes=genotypes, ALLELE=alleles, expr=expression)
        return (df)
    }

    binnedData=lapply(b, buildData)
    allCellsData=list(buildData(allCells))
    data=c(allCellsData, binnedData)

    if (nBins==3) {
        names(data)=c("All cells", "Bottom third", "Middle third","Top third")
    } else {
        names(data)=c("all cells", paste("Bin", 1:nBins))
    }

    #calculate bTotal and other metadata here:
    bTotal=NULL
    if (!is.null(summarizedResult)) {
        bTotalAll=getBinnedBtotal(gene, factorName, summarizedResult, usageDF, nBins=1)
        bTotal=getBinnedBtotal(gene, factorName, summarizedResult, usageDF, nBins)
        bTotal=c(bTotalAll, bTotal)
        names(bTotal)=names(data)
    }
    result=list(data=data, summary=list(bTotal=bTotal))
    return (result)
}

plotUmapFactor<-function (factorName, usageDF, uMapCoordinates) {
    plotLatentFactorOnUmap(gene=factorName, snp="", factorName="", usageDF, uMapCoordinates)
}

plotUmapBtotal<-function (gene, snp, factorName, summarizedResult, usageDF, uMapCoordinates, ...) {
    btotal=getBTotal (geneNames=gene, factorName, summarizedResult, usageDF)
    plotLatentFactorOnUmap(gene, snp, factorName, btotal, uMapCoordinates, ...)

}

plotLatentFactorOnUmap<-function (gene, snp, factorName="FACTOR_2", btotal, uMapCoordinates, titleSize=24, legendKeySize=0.5, legendTitleStr="", base_size=11) {

    score=btotal[,gene]
    strTitle=paste(gene, snp, factorName)
    p=ggplot2::ggplot(uMapCoordinates, aes(UMAP1, UMAP2)) +
        geom_point(size=0.002, aes(colour = score)) +
        xlab("UMAP1") + ylab("UMAP2") + labs(color=legendTitleStr) +
        viridis::scale_color_viridis() +
        #theme_classic(base_family = 'Courier') +
        theme_classic(base_size=base_size) +
        theme(plot.title = element_text(size = titleSize, face = "bold", hjust = 0.5), legend.key.size = unit(legendKeySize, 'cm'), axis.ticks = element_blank()) +
        ggtitle(strTitle) +
        scale_y_continuous(breaks=NULL) +
        scale_x_continuous(breaks=NULL)
    return (p)
}


#QQ plot of the tested
plotEmpiricPvalueDistribution<-function (r, qvalueThreshold=0.05,strPrefix="") {
    #how many hits?
    numHits=length(which(r$qvalue<=qvalueThreshold))
    qqman::qq(r$lrt_pval)
    title (paste(strPrefix, numHits, "of", dim(r)[1], "eGenes with qvalue <=", qvalueThreshold))
}

plotPNG<-function (pngFile, xleft=0, ybottom=0, xright=1, ytop=1) {
    op <- par(no.readonly = TRUE)
    r=png::readPNG(pngFile)
    par(mai=c(0,0,0,0))
    plot(c(0,1),c(0,1),type="n", axes=F)
    rasterImage(r,xleft, ybottom, xright, ytop)
    par(op)
}




###########################################################
# BASELINE eQTL plots - compare single cell to pseudobulk
###########################################################

plotBaseLineZscore<-function (df, qvalueThreshold=1, titleStrPrefix="", xyRange=NULL, ...) {
    df2=df[df$bulkQval<=qvalueThreshold,]

    b=df2$bulkZscore
    s=df2$singleCellZscore

    if (is.null(xyRange)) xyRange=range(c(b, s))
    correlation=cor(b,s)

    plot (b, s, xlab="pseudobulk eQTL discovery Z-score", ylab="single cell eQTL discovery Z-score",
          main=paste(titleStrPrefix, " Comparison of Z-scores correlation [ ", round (correlation,3), "]", sep=""),
          cex=0.25, xlim=xyRange, ylim=xyRange, ...)

    #points(b[df2$bulkQval<0.05], s[df2$bulkQval<0.05], col="blue", cex=0.25)
    fit=lm(s ~ b)
    abline (fit, col="red", lty=2)
    abline (0,1, col='black')
    legend("topleft", legend=c("best linear fit ", "identity"), fill=c("red", "black"))
}

plotBaselineEqtlPvalues<-function (df, qvalueThreshold=1, titleStrPrefix="", xyMax=NULL, ...) {
    df2=df[df$bulkQval<=qvalueThreshold,]

    b=-log10(df2$bulkPval)
    s=-log10(df2$singleCellPval)

    if (is.null(xyMax)) xyMax=max(c(b, s))

    corLog10=cor(b,s)
    plot (b, s, xlab="pseudobulk eQTL discovery empiric pvalue", ylab="single cell eQTL discovery lrt pvalue",
          main=paste(titleStrPrefix, " Comparison of empiric pvalues correlation [ ", round (corLog10,3), "]", sep=""), cex=0.25, xlim=c(0, xyMax), ylim=c(0, xyMax), ...)

    fit=lm(s ~ b)
    abline (fit, col="red", lty=2)
    abline (0,1, col='black')
    legend("topleft", legend=c("best linear fit ", "identity"), fill=c("red", "black"))

}

#df is a data frame that has the following columns: bulkPval, bulkQval, singleCellPval, bulkEffectSize, singleCellEffectSize
plotBaselineEqtlSignTest<-function (df, qvalueThreshold=1, titleStrPrefix="", ...) {

    df2=df[df$bulkQval<=qvalueThreshold,]
    numSignAgree=length(which(sign(df2$bulkEffectSize)==sign(df2$singleCellEffectSize)))
    totalCount=dim(df2)[1]
    fracAgree=numSignAgree/totalCount
    strTitle=paste(titleStrPrefix, " Fraction of Genes with concordant effect direction [", round (fracAgree, 3), "]", sep="")
    plot (df2$bulkEffectSize, df2$singleCellEffectSize, cex=0.25, main=strTitle, xlab="bulk eQTL beta", ylab="single cell beta", ...)
    abline (h=0, col="red", lty=2)
    abline (v=0, col="red", lty=2)
}




#####################################
# DATA PROCESSING
#####################################

prepareDataForSingleCellEqtl<-function (snpFile, h5SeuratFile,
                                        cellMetaDataFile, donorCovariatesFile, permutedPseudobulkEqtlResultsFile,
                                        usageFile=NULL, qvalueThreshold=0.05, eQTLGeneList=NULL, chromosome=NULL,
                                        factorList=NULL, quantizationErrorThreshold=NULL) {

    log_info("Reading in data for single cell eQTL analysis")

    #validate inputs exist.
    r=sapply(c(snpFile, h5SeuratFile,
               cellMetaDataFile, donorCovariatesFile, permutedPseudobulkEqtlResultsFile, usageFile),
             validateSingleFileExists)

    #read in the results first, so the various inputs can be filtered to the subset of critical genes/SNPs.
    r=getBulkeQTLs(permutedPseudobulkEqtlResultsFile, qvalueThreshold, eQTLGeneList, chromosome)
    bulkResults=r$bulkResults
    geneList=r$geneList
    snpList=r$snpList

    #parse genotypes
    genotypeDF=parseGenotypes(snpFile, snpList)

    #get the list of donors to work with
    donorList = colnames (genotypeDF)[-1]

    #read in the single cell metrics data.
    metricsDF=getSingleCellMetrics (cellMetaDataFile, donorList, quantizationErrorThreshold)
    cellBarcodeList=metricsDF$CELL_BARCODE_FINAL

    #generate a new feature - scaled log of the UMIs per cell
    metricsDF$SCALED_LOG_UMIS=scale(log(metricsDF$num_retained_transcripts))

    #read in the expression data
    expData=getExpressionMatrix(h5SeuratFile, geneList, cellBarcodeList, !is.null(quantizationErrorThreshold))

    #order the metrics cells in the same ordering as the expression data.
    #slim down to the mapping from cell barcode to donor, and a few cell level features.
    #drop the PREFIX if the data has been collapsed.
    featuresToRetain=c("CELL_BARCODE", "CELL_BARCODE_FINAL", "DONOR", "SCALED_LOG_UMIS", "pct_mt", "num_retained_transcripts")
    if ("QUANT_ERROR" %in% colnames (metricsDF)) {
        featuresToRetain=c(featuresToRetain, c("QUANT_ERROR", "NUM_CELLS"))
    } else {
        featuresToRetain=c(featuresToRetain, "PREFIX")
    }
    #featuresToRetain=c(featuresToRetain, "PREFIX")
    metricsDF=metricsDF[match(colnames(expData), metricsDF$CELL_BARCODE_FINAL),featuresToRetain,with=F]
    cellBarcodeList=metricsDF$CELL_BARCODE_FINAL

    #get the latent factor scores matrix for the cells retained
    usage=parseCnmfFactor(usageFile=usageFile, cellBarcodeList)

    #parse donor covariates and add to cell level metrics
    metricsDF=addCovariatesToMetricsDf(metricsDF, donorCovariatesFile)

    #assert line up the genotypeDF and metricsDF by cell
    if (!all (colnames (expData)==metricsDF$CELL_BARCODE_FINAL))
        stop ("Genotype DF and metrics DF not in the same order, something went wrong!")

    #add the latent factors to the metrics data frame if provided
    if (!is.null(usage))
        metricsDF=addLatentFactorsToMetricsDF(metricsDF, usage, factorList)

    log_info("Data reading complete")

    return (list(metricsDF=metricsDF, bulkResults=bulkResults, genotypeDF=genotypeDF, expData=expData, usage=usage))
}
#
addLatentFactorsToMetricsDF<-function (metricsDF, usageDF, factorList) {
    #if no factor is selected use all factors.
    #TODO: experimental, subtract the means from the factors to remove structural correlation?
    if (is.null(factorList)) {
        #This doesn't help.  Design matrix still not full rank
        #log_info("Center usageDF")
        #usageDF=data.frame(scale(usageDF, center=T, scale=F))
        factorList=colnames(usageDF)
    }

    df=data.table::copy(metricsDF)
    df[,(factorList):=usageDF[df$CELL_BARCODE_FINAL ,factorList]]
    return (df)

}

parseGenotypes<-function (snpFile, snpList=NULL) {
    a=data.table::fread(snpFile, key="id")
    if (!is.null(snpList)) {
        missingSNPs=setdiff(snpList, a$id)
        if (length(missingSNPs)>0)
            stop("Missing SNPs in geneotype file [", paste(missingSNPs, collapse=","), "]")
        a=a[match(snpList, a$id),]
    }
    return (a)
}


#' Add donor level covariates to the single cell metrics data frame
#'
#' @param metricsDF A data frame containing cell level metrics.  Must contain a DONOR column.
#' @param donorCovariatesFile A donor covariates file in the eQTL format
#'
#' @return The metrics file with the additonal donor level metrics
#' @export
addCovariatesToMetricsDf<-function (metricsDF, donorCovariatesFile) {
    covs=parseDonorCovariates(donorCovariatesFile)
    missingCovs=setdiff(metricsDF$DONOR, rownames(covs))
    if (length(missingCovs)>0)
        stop("Missing donor entries in covariates file [", paste(missingCovs, collapse=","), "]")
    idx=match(metricsDF$DONOR, rownames(covs))
    metricsDF=cbind(metricsDF, covs[idx,])
    return (metricsDF)
}

#some conversion of the existing eQTL format
parseDonorCovariates<-function (donorCovariatesFile) {
    donorCovars=read.table(donorCovariatesFile, header=T, stringsAsFactors = F, sep="\t")
    donorCovars=t(donorCovars)
    colnames(donorCovars)=donorCovars[1,]
    donorCovars=donorCovars[-1,]
    result=apply(donorCovars, 2, as.numeric)
    rownames(result)=rownames(donorCovars)
    return (result)
}

correlateDonorCovariates<-function (donorCovariatesFile) {
    donorCovars=parseDonorCovariates(donorCovariatesFile)
    z=cor(donorCovars)
    heatmap(z)
}

#' Extracts the sparse expression data matrix from an h5 seurat file.
#'
#' @param h5SeuratFile The .h5seurat file to load
#' @param geneList An optional list of genes.  Data will be subset to that list of genes
#' @param cellBarcodeList A list of cell barcodes.  If filterExpression is false, check that these are in the expected data.
#' If filterExpression is true, limit the expression data to the intersect of the cell barcodes in the list and expression data.
#' @return a sparse matrix (dgCMatrix) of expression data.  Genes in rows, cells in columns.
#' @import SeuratDisk Seurat
#' @export
getExpressionMatrix<-function (h5SeuratFile, geneList=NULL, cellBarcodeList=NULL, filterExpression=F) {
    data <- SeuratDisk::LoadH5Seurat(h5SeuratFile, assays = "RNA", verbose =TRUE, mode="r")

    if (!is.null(geneList)) {
        #Seurat doesn't support "_" in gene names.
        idxUnderscoreGenes=grep ("_", geneList)
        if (length(idxUnderscoreGenes)>0)
            log_warn("Some genes in the analysis had their names altered by Seurat. [",paste(geneList[idxUnderscoreGenes], collapse=","), "]")
        data=data[geneList,]
    }

    if (!is.null(cellBarcodeList)) {
        if (filterExpression) {
            both=intersect (colnames (data), cellBarcodeList)
            data=data[,both]
        }
        missingMetrics=setdiff(colnames (data), cellBarcodeList)
        if (length(missingMetrics)>0)
            log_warn("Some cell barcodes in the expression data had no metadata. [",paste(missingMetrics, collapse=","), "]")
    }

    expData=Seurat::GetAssayData(data)
    rm (data)

    return (expData)
}

#' Get cell level metrics, optionally filtered
#'
#' @param cellMetaDataFile A matrix of cell level metadata
#' @param donorList An optional list of donors to subset the cell level metadata.
#' @param quantizationErrorThreshold If the cell level metadata contains the column QUANT_ERROR, the data has been previously compressed, and will be filtered to cells with values less than this score.
#' @return A dataframe containing cell level metadata.
#' @export
getSingleCellMetrics<-function (cellMetaDataFile, donorList=NULL, quantizationErrorThreshold=NULL) {
    metricsDF=data.table::fread(cellMetaDataFile)
    if (!is.null(donorList)) {
        missingDonors=setdiff(donorList, metricsDF$DONOR)
        if (length(missingDonors)>0) {
            stop("Donors requested but not found in cell level metrics [", paste(missingDonors, collapse=","), "]")
        }
        idx=metricsDF$DONOR %in% donorList
        metricsDF=metricsDF[idx,]
    }

    if (!is.null(quantizationErrorThreshold)) {
        if (!"QUANT_ERROR" %in% colnames (metricsDF))
            stop("Requested filtering on QUANT_ERROR, but column not found in cellMetaDataFile [", cellMetaDataFile, "]")
        log_info("Filtering Cell Metrics to cells with a quantization error < ", quantizationErrorThreshold)
        metricsDF=metricsDF[metricsDF$QUANT_ERROR<=quantizationErrorThreshold,]
        log_info("Median quantitative error [", round (median (metricsDF$QUANT_ERROR), 3), "]")
    }

    metricsDF$CELL_BARCODE_FINAL=paste(metricsDF$PREFIX, metricsDF$CELL_BARCODE, sep="_")
    return (metricsDF)
}

#an eQTL gene list will override the qvalue threshold.
getBulkeQTLs<-function (permutedPseudobulkEqtlResultsFile=NULL, qvalueThreshold=0.05, eQTLGeneList=NULL, chromosome=NULL) {
    #short circuit
    if (is.null(permutedPseudobulkEqtlResultsFile)) {
        result=list(bulkResults=NULL, geneList=NULL, snpList=NULL)
        return (result)
    }
    bulkResults=data.table::fread(permutedPseudobulkEqtlResultsFile)
    #phenotype_id is tensorQTL data.
    if (colnames (bulkResults)[1]=="phenotype_id")
        bulkResults=remapTensorQTLColumnNames(bulkResults)

    if (!is.null(eQTLGeneList)) {
        missing=setdiff(eQTLGeneList, bulkResults$gene)
        if (length(missing)>0)
            stop ("Requested genes from bulk results for analysis that were not tested [", paste(missing, collapse=","), "]")
        bulkResults=bulkResults[match(eQTLGeneList, bulkResults$gene)]
        log_info("Loaded [", dim (bulkResults)[1], "] significant eQTLs from specified gene list")
    }
    if (!is.null(chromosome)) {
        chr=sapply(strsplit (bulkResults$SNP, split=":", fixed=T), function (x) x[1])
        if (!chromosome %in% chr)
            stop (paste("Requested chromosome [", chromosome, "] not found in chromosomes [", paste(unique(sort(chr)), collapse=","), "]"))
        bulkResults=bulkResults[chr %in% chromosome,]
    }

    if (!is.null(qvalueThreshold)) {
        bulkResults=bulkResults[bulkResults$qvalue<= qvalueThreshold,]
        log_info("Loaded [", dim (bulkResults)[1], "] significant eQTLs")
    } else {
        log_info("Loaded [", dim (bulkResults)[1], "] eQTLs")
    }



    result=list(bulkResults=bulkResults, geneList=sort(unique(bulkResults$gene)), snpList=sort(unique(bulkResults$SNP)))
    return (result)
}

#remap column names for SNP, snp, beta, p-value, qvalue
remapTensorQTLColumnNames<-function (bulkResults) {
    colnames (bulkResults) [which(colnames(bulkResults)=="phenotype_id")]<-"gene"
    colnames (bulkResults) [which(colnames(bulkResults)=="variant_id")]<-"SNP"
    colnames (bulkResults) [which(colnames(bulkResults)=="slope")]<-"beta"
    colnames (bulkResults) [which(colnames(bulkResults)=="pval_nominal")]<-"p-value"
    colnames (bulkResults) [which(colnames(bulkResults)=="qval")]<-"qvalue"
    return (bulkResults)
}

parseGeneLocations<-function (geneLocationFile, geneList=NULL) {
    a=data.table::fread(geneLocationFile)
    if (!is.null(geneList)) {
        missingGenes=setdiff(geneList, a$geneid)
        if (length(missingGenes)>0)
            stop("Missing genes in gene location file [", paste(missingGenes, collapse=","), "]")
        a=a[match(geneList, a$geneid),]
    }
    colnames (a)=c("gene", "chr", "start", "end")
    return (a)
}

#All cells in the barcode list are expected in this data set if a list is provided.
#A cell barcode list will sort the usage in the same order.
parseCnmfFactor<-function (usageFile=NULL, cellBarcodeList=NULL) {
    if (is.null(usageFile)) return (NULL)
    u=read.table(usageFile, header=T, stringsAsFactors=F, sep="\t", row.names=1)
    #if there are multiple factors, rescale the factors to sum to 1.
    #TODO: experiment with removing scaling to sum to 1
    # if (dim (u)[2]>1)
    #      u=u/rowSums(u)
    #TODO: experiment with scale/center factor scores
    #maybe not the best, but all the linear scalings are have the same result.
    u=scale(u)
    u=data.frame(u, stringsAsFactors = F)
    if (!is.null(cellBarcodeList)) {
        missingBarcodes=setdiff(cellBarcodeList, rownames (u))
        if (length(missingBarcodes)>0)
            stop("cNMF cell level scores missing some barcodes. Examples: [", paste(head(missingBarcodes), collapse=","), "]")
        u=u[cellBarcodeList,,drop=F]
    }
    #if any columns don't match "FACTOR" then rename
    if (!all(grepl("FACTOR", colnames (u))))
        colnames(u)=paste ("FACTOR", 1:dim(u)[2], sep="_")
    return (u)
}

parseUmapCoordinates<-function (uMAPCoordinatesFile, cellBarcodeList=NULL) {
    uMapCoordinates=read.table(uMAPCoordinatesFile, header=T, stringsAsFactors = F, sep="\t")
    cellsFixed=stringr::str_replace(string = uMapCoordinates$cells, pattern = "-[:digit:]+$", replacement = "")
    uMapCoordinates$cells=cellsFixed

    #optionally filter and order on an existing cell barcode list
    if (!is.null(cellBarcodeList)) {
        missing=setdiff(cellBarcodeList, uMapCoordinates$cells)
        if (length(missing)>0)
            stop(paste("uMAP coordinates missing requested cell barcodes [", paste(head(missing), collapse=","), "]"))
        uMapCoordinates=uMapCoordinates[match(cellBarcodeList, uMapCoordinates$cells),]
    }
    return (uMapCoordinates)
}


#IE: make metacells.
aggregateCountsByDonor2<-function (h5SeuratFile, cellMetaDataFile) {
    metricsDF=getSingleCellMetrics (cellMetaDataFile)
    data <- SeuratDisk::LoadH5Seurat(h5SeuratFile, assays = "RNA", verbose =TRUE, mode="r")

    #subset the data.
    cellsBoth=intersect(colnames(data), metricsDF$CELL_BARCODE_FINAL)
    data <- subset(data, cells = cellsBoth)

    data$donor=metricsDF[match(colnames(data), metricsDF$CELL_BARCODE_FINAL),]$DONOR
    metacell=AggregateExpression(data, group.by="donor", slot="counts")
    metacell=metacell$RNA
    metacell=data.frame(cbind(GENE=rownames(metacell), metacell), stringsAsFactors = F, row.names = NULL)
    #write.table(metacell, "/broad/mccarroll/nemesh/BA46.n191.raw.dge.astrocyte.All.noOutliers.metacells.txt", row.names = F, col.names = T, quote=F, sep="\t")

}


#https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks
chunk2 <- function(x,n) {
    if (n==1) return (list(x))
    split(x, cut(seq_along(x), n, labels = FALSE))
}


###############################################
# DATA REDUCTION
################################################


#' Reduce the single cell representation to a less granular data set to reduce computational burden
#'
#' This is the command line version of reduceData that reads in inputs from files and emits flat files
#' containing the new expression data, cell features, and factor usage generated by data collapse.
#'
#' These new files can then be 1-1 replacements for the original inputs and used in analysis.
#'
#' The permutedPseudobulkEqtlResultsFile provides a pre-filtered set of genes to operate on.  This
#' can additionally be filtered to a specific chromosome, gene list, or qvalue threshold.
#'
#' @param outCellFeaturesFile A cell level metrics file of the compressed cells
#' @param outExpressionFile A seurat h5 file of expression for the compressed cells
#' @param outUsageFile A latent factor matrix containing the selected factor for the compressed cells
#' @export
#' @inheritParams reduceData
reduceDataCmdLine<-function (snpFile, h5SeuratFile, cellMetaDataFile, donorCovariatesFile,
                             usageFile, factorName, useAllLatentFactors=F, usePeerFactors=F,
                             numClustersPerDonor=20, permutedPseudobulkEqtlResultsFile,
                             qvalueThreshold=1, eQTLGeneList=NULL, chromosome=NULL,
                             reduceDataPDF=NULL, outCellFeaturesFile, outExpressionFile, outUsageFile) {

    #allow the factor list be null to load all factors.
    factorList=factorName
    if (useAllLatentFactors)
        factorList=NULL

    #if usageFile is not null and factorList is null, all factors are loaded
    d=prepareDataForSingleCellEqtl(snpFile, h5SeuratFile,
                                   cellMetaDataFile, donorCovariatesFile, permutedPseudobulkEqtlResultsFile,
                                   usageFile, qvalueThreshold, eQTLGeneList, chromosome,
                                   factorList=factorList)

    metricsDF=d$metricsDF
    expData=d$expData
    usageDF=d$usage

    #TODO: remove hardcoded fixed features?
    featureList=getFeatureList (fixedFeatures=c("SCALED_LOG_UMIS", "pct_mt"), factorName, usePeerFactors, metricsDF, useAllLatentFactors, usageDF)
    log_info(paste("Feature List [", paste(featureList, collapse=" "), "]"))
    r=reduceData(numClustersPerDonor, factorName, featureList, metricsDF, expData, usageDF, reduceDataPDF)
    data.table::setDF(r$metricsDF)
    #strip the cluster from the data table, and make a new usage matrix to write out.
    #TODO: this is a problem when useAllLatentFactors is TRUE.

    if (useAllLatentFactors) {
        usageDFNew=r$metricsDF[,colnames(usageDF),drop=F]
    } else {
        usageDFNew=r$metricsDF[,c(factorName),drop=F]
    }

    rownames(usageDFNew)=r$metricsDF$CELL_BARCODE_FINAL

    data.table::fwrite(usageDFNew, outUsageFile, row.names=T, col.names = T, quote=F, sep="\t")
    #z=parseCnmfFactor(usageFile=outUsageFile, cellBarcodeList=NULL)
    #all.equal(usageDF, z)

    data.table::fwrite(r$metricsDF, outCellFeaturesFile, row.names = F, col.names = T, quote=F, sep="\t")
    #written as a dense matrix - matrix is ~ 63% sparse, but that will change by K.
    #data.table::fwrite(as.matrix(r$expData), outExpressionFile, row.names = T, col.names = T, quote=F, sep="\t")

    seu <- CreateSeuratObject(counts=r$expData)
    SaveH5Seurat(seu, outExpressionFile, overwrite = TRUE)

}

#' Reduce the single cell representation to a less granular data set to reduce computational burden
#'
#' This leverages vector quantization to split each donor's cells into N groups where cells in a group have
#' similar scores for known covariate information.  In effect, this bins each donor's total cells into N bins of data representing
#' the cluster center for the group
#'
#' # The expression data is modified to reflect the median of the cells in the group.  This preserves the counts nature of the data.
#' The expression data is modified to reflect the sum of cells in the group.  The number of total transcripts is recalculated to the sum of the cells
#' in the group, and the SCALED_LOG_UMIS feature is recalculated for the new exemplar cells.
#'
#' The cell barcode names for the metricsDF and expression data are changed to reflect the "bins" of donor expression.
#'
#' @param numClustersPerDonor The number of cells to emit for each donor.  These cells are the "average" cells for each cluster.
#' @param featureName A column in the metricsDF.  This feature will be plotted in the output PDF,
#' and a measurement of the error of each cell to it's median will be calculated.  This only affects plotting, not the data outputs.
#' This is typically a latent factor of interest.  This parameter only affects the error calculation.
#' This feature should also be in the featureList (and will be added if it is not.)
#' @param featureList A vector of features that will used to calculate the distances between cells.
#' This list typically includes latent factor(s) of interest and known covariates.  All features
#' must be numeric.
#' @param metricsDF A data frame containing the features for each cell barcode.
#' @param expData A matrix containing expression data for each cell - cells in rows, genes in columns.
#' @param usageDF A data frame containing the latent factor usage per cell.  Each row repesents a single cell,
#' columns represent each factor score.  Cell IDs are stored in the rownames of the dataframe.
#' @param reduceDataPDF An optional PDF to emit plots to.  Each plot contains a summary of how the cells of a donor
#' were reduced to a set of exemplar cells.
#'
#' @return A list containing the modified expression and metrics dataframes.  These are drop-in replacements for the
#' original expData and metricsDF data that was passed in.
reduceData<-function (numClustersPerDonor=NULL, featureName, featureList, metricsDF, expData, usageDF, reduceDataPDF=NULL) {

    #make sure the featureName is in the feature list.
    if (!featureName %in% featureList)
        featureList=c(featureName, featureList)

    #validate all features are already in the metricsDF
    missingFeatures=setdiff(featureList, colnames(metricsDF))
    if (length(missingFeatures)>0)
        stop ("Requested features missing from cell level metrics [", paste(missingFeatures, collapse=""), "]", sep="")

    donorList=sort(unique(metricsDF$DONOR))
    showPlots=F
    if (!is.null(reduceDataPDF)) {
        showPlots=T
        pdf(reduceDataPDF)
    }

    log_info("Generating [", numClustersPerDonor, "] prototype cells per donor")

    r<-lapply(donorList, reduceDataDonor, metricsDF, featureList, nclust=numClustersPerDonor, normalize=T, featureName=featureName, showPlots=showPlots)
    if (!is.null(reduceDataPDF)) dev.off()

    #rbind the results
    summary=do.call(rbind, lapply(r, function (x) x$summary))
    exemplarCells=do.call(rbind, lapply(r, function (x) x$exemplarCells))
    clusterMap=do.call(rbind, lapply(r, function (x) x$clusterMap))

    #Sum the UMIs across cells that are grouped with the donor, then recalculate the total transcripts and scaling factor.
    log_info("Calculating expression for exemplar cells")
    expDataNew<-calculateExemplarExpression(exemplarCells, expData, clusterMap)
    log_info("Finished calculating expression for exemplar cells")

    #validate total expression of each gene remains the same, but it's distributed differently.
    #if (rowSums(expData)!=rowSums(expDataNew))
    #    stop ("Something went wrong")

    #recalculate the num_retained_transcripts and SCALED_LOG_UMIS
    exemplarCells=recalculateTotalTranscripts(exemplarCells, clusterMap, metricsDF)
    exemplarCells$SCALED_LOG_UMIS=scale(log(exemplarCells$num_retained_transcripts))

    #the metrics DF should not have a CELL_BARCODE_FINAL column.  It should have a CELL_BARCODE column where the prefix isn't appended.
    #can recover the cell barcode from the original metricsDF
    exemplarCells$CELL_BARCODE=metricsDF[match(exemplarCells$CELL_BARCODE_FINAL, metricsDF$CELL_BARCODE_FINAL),]$CELL_BARCODE

    #reorder columns to look like the original data, PLUS extra columns.
    colsOrder=c(intersect (colnames(metricsDF), colnames(exemplarCells)), setdiff(colnames(exemplarCells), colnames(metricsDF)))
    exemplarCells=exemplarCells[,colsOrder,with=F]

    result=list(metricsDF=exemplarCells, expData=expDataNew)
    return (result)
}

#After a lot of testing of multiple ways to calculate this, two things are clear:
#1) Summarizing by the sum of expression is the way to go - if you summarize by the mean/median, the sparsity of the data
# means that many genes have an average expression of 0 across all donor meta cells.  This isn't very good.
#2) Of the ways I've thought of to aggregate the data (rowSums, multiplying the sparse matrix by a design matrix), Seurat's aggregate
#is orders of magnitude faster and doesn't have a large impact on memory, so stick with that!
calculateExemplarExpression<-function (exemplarCells, expData, clusterMap, sparseMatrixOutput=T) {

    #for each donor:cluster in the map, add the exemplar cell the data will map back to.
    mapExemplarBarcode<-function (clusterMap, exemplarCells) {
        clusterMap2=clusterMap[match(colnames(expData), clusterMap$barcode),]
        clusterMap2$key=factor(paste(clusterMap2$donor, clusterMap2$cluster, sep=":"))
        exemplarCellsKey=paste(exemplarCells$DONOR, exemplarCells$cluster, sep=":")
        idx=match(clusterMap2$key, exemplarCellsKey)
        clusterMap2$id=exemplarCells[idx,]$CELL_BARCODE_FINAL
        return (clusterMap2)
    }

    clusterMap2=mapExemplarBarcode (clusterMap, exemplarCells)
    clusterMap2=clusterMap2[match(colnames (expData), clusterMap2$barcode),]

    d <- CreateSeuratObject(counts=expData)
    #use the id of the exemplar cell to group cells via seurat.
    d$ident=clusterMap2$id
    #this is already transposed!
    system.time (exp<-Seurat::AggregateExpression(d, group.by="ident", slot="counts"))
    exp=Matrix::Matrix(exp$RNA, sparse = sparseMatrixOutput)
    return (exp)

}

recalculateTotalTranscripts<-function (exemplarCells, clusterMap, metricsDF) {
    #exemplarCell=exemplarCells[1,]$CELL_BARCODE_FINAL
    recalculateTotalTranscripts<-function (exemplarCell, clusterMap, metricsDF) {
        x=clusterMap[clusterMap$barcode==exemplarCell,]
        #log_info("Cell [", exemplarCell, "] Donor [", x$donor, "] ", "cluster [", x$cluster, "]")
        cells=clusterMap[clusterMap$cluster==x$cluster & clusterMap$donor==x$donor,]$barcode
        total=sum (metricsDF[match(cells, metricsDF$CELL_BARCODE_FINAL),]$num_retained_transcripts)
        return (total)
    }
    total=sapply(exemplarCells$CELL_BARCODE_FINAL, recalculateTotalTranscripts, clusterMap, metricsDF)
    exemplarCells$num_retained_transcripts=total
    return (exemplarCells)
}


#donor="AN00216"; nclust=20
reduceDataDonor<-function (donor, metricsDF, featureList, nclust=10, normalize=T, featureName="FACTOR_2", showPlots=F, quant_method="kmeans") {
    #slim the data to a single donor, restrict to the features for the vector quantization + identifier.
    log_info("Finding exemplar cells for donor [", donor, "]")
    idxDonor=which(metricsDF$DONOR==donor)
    df=metricsDF[idxDonor,]
    data.table::setDF(df, rownames = df$CELL_BARCODE_FINAL)
    df=df[,featureList]

    #if nclust >= number of cells, then return the original data.
    if (nclust>=dim(df)[1]) {
        log_info("Not enough cells to generate [", nclust, "] exemplars")
        #hard code all cells to belong to cluster 1.
        summary=data.frame(donor=donor, avgError=NA, stringsAsFactors = F)
        z=metricsDF[idxDonor,]
        z$cluster=1:dim(z)[1]
        z$QUANT_ERROR=NA
        z$NUM_CELLS=1
        clusterMap=data.frame(donor=donor, cluster=1:dim(z)[1], barcode=z$CELL_BARCODE_FINAL)
        result=list(summary=summary, exemplarCells=z, clusterMap=clusterMap)
        return (result)
    }

    if (quant_method=="kmedoids") {
        suppressMessages(hvt.results<-muHVT::HVT(dataset=df, nclust=nclust, depth=1, quant.err = 0.2, projection.scale = 10,
                                                 normalize = normalize, distance_metric = "L1_Norm", error_metric = "mean", quant_method="kmedoids"))
        codeBook=hvt.results[[3]][['summary']]

        #use predict to get back for each data point what cluster it belongs to.
        suppressMessages(r<-muHVT::predictHVT(df, hvt.results, error_metric="mean",  distance_metric = "L1_Norm")$scoredPredictedData)
        rownames (r)=rownames (df)

        df$assigned_cluster=factor(r$Segment.Child)
        #get the cluster center cell barcodes
        medoids=r[r$Quant.Error==0,]
        medoids=medoids[order(medoids$Segment.Child),]
        medoidCellBarcodes=rownames (medoids)

        #these are the final "exemplar" cells with their original values.
        exemplarCells=metricsDF[match(medoidCellBarcodes, metricsDF$CELL_BARCODE_FINAL),]
        exemplarCells$cluster=medoids$Segment.Child

    }

    if (quant_method=="kmeans") {
        suppressMessages(hvt.results<-muHVT::HVT(dataset=df, nclust=nclust, depth=1, quant.err = 0.2, projection.scale = 10,
                                                 normalize = normalize, distance_metric = "L1_Norm", error_metric = "mean", quant_method="kmeans"))
        codeBook=hvt.results[[3]][['summary']]

        #use predict to get back for each data point what cluster it belongs to.
        suppressMessages(r<-muHVT::predictHVT(df, hvt.results, error_metric="mean",  distance_metric = "L1_Norm")$scoredPredictedData)
        rownames (r)=rownames (df)

        #reverse the scale/center on the codebook back to the original values
        #fn="SCALED_LOG_UMIS"
        reverseScale<-function (fn, codeBook, dfScaled) {
            (codeBook[,fn] * attr(dfScaled, 'scaled:scale')[fn]) + attr(dfScaled, 'scaled:center')[fn]
        }
        dfScaled=scale(df)
        codeBook[,featureList]=do.call(cbind, lapply(featureList, reverseScale, codeBook, dfScaled))

        # get a representitive cell barcode from the donor/cluster
        # so that other data can look up genotypes/etc for this cell/donor
        getStandInCellBarcodes<-function (clusterID, r) {
            rr=r[r$Segment.Child==clusterID,]
            rownames (rr[which.min(rr$Quant.Error),])
        }

        #get the cell that is closest to the minimum.
        clusterIDs=sort(unique(r$Segment.Child))
        kMeansCellBarcodes=sapply(clusterIDs, getStandInCellBarcodes, r)

        #these are the final "exemplar" cells with their original values.
        exemplarCells=metricsDF[match(kMeansCellBarcodes, metricsDF$CELL_BARCODE_FINAL),]
        exemplarCells$cluster=clusterIDs

        exemplarCells[,featureList]=codeBook[,featureList]

        #also copy the quantization error
        exemplarCells$QUANT_ERROR=codeBook$Quant.Error
        exemplarCells$NUM_CELLS=codeBook$n

        #add the cluster label to the df
        df$assigned_cluster=factor(r$Segment.Child)

    }

    #calculate the sum of the distances from the center point to the cells for each cluster [codeword]
    avgError=calculateAvgError(r, df, exemplarCells, featureName)
    summary=data.frame(donor=donor, avgError=avgError, stringsAsFactors = F)

    if (showPlots) reduceDataDonorPlots(nclust, hvt.results, df, exemplarCells, featureName, featureList, donor)

    #need a mapping from the prototype cell to all cell barcodes the cell is related to.
    clusterMap=data.frame(donor=donor, cluster=r$Segment.Child, barcode=rownames (r))
    clusterMap=clusterMap[order(clusterMap$cluster),]

    #return a list containing the summary and the downsampled cell barcodes.
    result=list(summary=summary, exemplarCells=exemplarCells, clusterMap=clusterMap)
    return (result)
}

reduceDataDonorPlots<-function (nclust, hvt.results, df, exemplarCells, featureName="FACTOR_2", featureList, donor) {

    #bind up the cluster ID to the exemplar cells.
    exemplarDF=df[exemplarCells$CELL_BARCODE_FINAL,]

    #plot the distribution by reaction
    # df$assigned_cluster=r$Segment.Child
    # df$PREFIX=unlist(metricsDF[idxDonor, "PREFIX"])
    # par(mar=c(10,4,4,4))
    # vioplot::vioplot (df$assigned_cluster ~ df$PREFIX, las=2, xlab="", ylab="assigned cluster", main=donor, col="light blue")
    # par(mar=c(5, 4, 4, 2) + 0.1)

    plotQuantErrorByCluster<-function (codeBook) {
        ggplot(codeBook, aes(n, Quant.Error)) +
            geom_text(aes(label = Segment.Child)) +
            xlab("Number of cells in code word") +
            ylab("Quantization Error") +
            theme_classic(base_size = 14)
    }

    #plot the distribution by codeword (cluster center)
    #https://stackoverflow.com/questions/5677885/ignore-outliers-in-ggplot2-boxplot
    plotFactorByCluster<-function (df, featureName, donor) {
        z=data.frame(cluster=factor(df$assigned_cluster), score=df[[featureName]])
        sortedClusterOrder=sort(tapply(z$score, INDEX=z$cluster, FUN = median))
        z$cluster=factor(df$assigned_cluster, levels=c(names(sortedClusterOrder)))

        #exemplarScores=data.frame(cluster=as.numeric(names(sortedClusterOrder)))
        #exemplarScores$score=exemplarDF[exemplarScores$cluster,][[featureName]]

        #x=reorder(cluster, score, FUN=median)
        ggplot(z, aes(x=cluster, y=score)) +
            #geom_violin(fill="light blue") +
            #geom_point(data = exemplarScores, aes(x = cluster, y = score)) +
            geom_boxplot(fill="light blue", varwidth = TRUE, outlier.shape = NA) +
            #stat_summary(fun.y=median, geom="point", size=2, color="red") +
            theme_classic(base_size = 10) +
            #geom_boxplot(outlier.shape = NA) +
            #theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
            ylab(featureName) +
            xlab("cluster") +
            ggtitle(donor)
            #scale_y_continuous(limits = quantile(z$score, c(0, 0.9)))
    }

    z=muHVT::hvtHmap(hvt.results, df, child.level = 1, hmap.cols = "Quant.Error", line.width = c(0.2),
               color.vec = c("#141B41"), palette.color = 6, centroid.size = 1.5, show.points = T, quant.error.hmap = 0.2, nclust.hmap = nclust)

    # z=hvtHmap(hvt.results, df, child.level = 1, hmap.cols = featureName, line.width = c(0.2),
    #          color.vec = c("#141B41"), palette.color = 6, centroid.size = 1.5, show.points = T, quant.error.hmap = 0.2, nclust.hmap = nclust)

    p1<-plotFactorByCluster(df, featureName, donor)
    p2<-plotJointFeatureDistributions2(featureList, df, exemplarDF)

    topPlot=cowplot::plot_grid(p1, z, ncol=2)
    finalPlot=cowplot::plot_grid(topPlot, p2, nrow=2)
    print (finalPlot)

    #TODO: plot quantization error by number of points in cluster.  Label by cluster ID?
    #plotQuantErrorByCluster(codeBook)


}
#compute the average error for the featureName
#that is, the average absolute distance from the cluster center to the point.
calculateAvgError<-function (r, df, exemplarCells, featureName="FACTOR_2") {

    #cellBarcodeMediod=medoidDF[1,]$CELL_BARCODE_FINAL
    getTotalDistanceCluster<-function (cellBarcodeMediod, r, df, featureName) {
        #get the cluster id for the cluster center, find all the cells in the cluster
        #and compute the distance from the cluster center to each cell for the given field.
        clusterID=r[cellBarcodeMediod,]$Segment.Child
        allCellsInCluster=rownames (r[r$Segment.Child==clusterID,])
        allButExemplarCB=setdiff(allCellsInCluster, cellBarcodeMediod)
        #short circuit if this cell barcode clustered by itself.
        if (length(allCellsInCluster)==0)
            return (NA)
        centerScore=df[cellBarcodeMediod, featureName]
        otherScores=df[allButExemplarCB, featureName]
        return (mean (abs(centerScore-otherScores)))
    }

    mean(sapply(exemplarCells$CELL_BARCODE_FINAL, getTotalDistanceCluster, r, df, featureName), na.rm=T)

}

plotFeatureDistributions<-function (featureList, df, strTitle="") {
    plotOne<-function (featureName, df) {
        plot(density(df[[featureName]]), xlab="", ylab="", main="")
        title(main=featureName, cex.main = 1, line=0.5)
    }
    nrows=ceiling(length(featureList)/2)
    par(mfrow=c(nrows, 2), mar=c(2,2.5,2.5,1), mpg=c(2,0,0))
    lapply(featureList, plotOne, df)
    mtext(strTitle, side = 3, line = - 1.25, outer = TRUE, cex=1.25)
    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1, mpg=c(3,1,0))
}

#df1=df; df2=exemplarDF
plotJointFeatureDistributions2<-function (featureList, df1, df2) {
    #TODO: plot both distributions at the same time using ggplot2
    #This will keep the density on the same scale on the X.
    #Layout as length(featureList) columns, 1 row.
    df1$source="all"
    df2$source="codeword"
    both=rbind(df1, df2)

    plotOne<-function (f, both) {
        ggplot(data=both, aes(x=.data[[f]], fill=source)) +
            geom_density(alpha=0.5) +
            ggtitle(f) +
            theme(legend.position="bottom", legend.title=element_blank(), axis.title.y = element_blank(),
                  axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
            scale_fill_manual( values = c("red","blue"))

    }
    l=lapply(featureList, plotOne, both)
    z=cowplot::plot_grid(plotlist=l, nrow=1)
    return (z)
}


###########################################################
#
#gene="RERE";snp="chr1:8388317:C:CTT";

# use the null model to get the fitted residuals
getResiduals<-function (metricsDF, bulkResults, expData, usePeerFactors=T) {
    #set up a dataframe with slots for all the analysis.
    #this may be faster to run many times if we use data table to change column values per snp/gene tested.
    df=data.table::copy(metricsDF)
    df[,expression:=as.numeric (expData[gene,])]

    nullModelStr=paste("expression ~ ", "(1 | DONOR) + (1 | PREFIX)+ SCALED_LOG_UMIS + pct_mt", sep=" ")

    if (usePeerFactors) {
        peerFactors=colnames(df)[grep ("PEER", colnames (df))]
        log_info("Using [", length(peerFactors), "] PEER factors")
        nullModelStr=paste(nullModelStr, paste(peerFactors, collapse=" + "), sep=" + ")
    }

    nullModelFormula=as.formula(nullModelStr)
    null_model <- lme4::glmer(formula=nullModelFormula, data=df, control = lme4::glmerControl(optimizer = "nloptwrap"), family = "poisson", nAGQ=0)

    #residuals are normally distributed around 0.
    res=residuals(null_model)
    names (res)=colnames (expData)
    return (res)
}



# scatterPlotExpressionBinnedByFactorScore<-function (gene, snp, factorName, expData, metricsDF, usageDF, genotypeDF, nBins=3, method=c("pseudobulk", "avgLog2PlusOne")) {
#     r=getExpressionBinnedByFactorScore(gene, snp, factorName, expData, metricsDF, usageDF, genotypeDF, nBins, method)
#     data=r$data
#     #extract all to get the min/max info, etc.
#     all=data[names (data) %in% "all cells"][[1]]
#     maxY=max (all$expr)
#
#     #remove the "all cells" data set.
#     data=data[!names (data) %in% "all cells"]
#     cols=brewer.pal("Set2", n=3)
#
#     addSeries<-function (index, data) {
#         d=data[[index]]
#         points (jitter(d$genotypes), d$expr, col=cols[index], cex=0.75)
#         abline (lm(d$expr~ d$genotypes), col=cols[index], lty=2, lwd=2)
#     }
#
#     plot (c(-0.3,2.3), c(0, maxY), xlab="", ylab="Normalized expression", type='n', axes=F)
#     axis(side=1, at=0:2, labels=levels(all$ALLELE), col = NA, col.ticks = NA)
#     axis(side=2)
#     sapply(1:3, addSeries, data)
#     legend("topleft", legend=c("low", "medium", "high"), fill=cols, title=paste(factorName, "score"), ncol=nBins)
#     title(paste0(gene, ' / ', snp))
# }


# Despite claims to the contrary, this is slow AF.
# https://github.com/fabsig/GPBoost/issues/21
# runManyGenotypeRegressionsGPBoost<-function (metricsDF, bulkResults, genotypeDF, expData, usePeerFactors=T) {
#     #use GPBoost package instead.
#
#     #group_data is a matrix or vector with categorical grouping variable(s) specifying the random effects
#     #structure. If there are multiple (crossed or nested) random effects, the corresponding grouping
#     #variables should be in the columns of group_data
#     #y is a vector with response variable data
#     #X is a matrix with fixed effects covariate data
#
#     df=data.table::copy(metricsDF)
#     df[,genotype:=NA]
#     df[,expression:=NA]
#
#     idxToGenotype=match(metricsDF$DONOR, colnames(genotypeDF))
#
#     #group data is fixed across all tests.
#     group_data=as.matrix (df[,c("DONOR", "PREFIX")])
#
#     #X has all columns except the genotype.
#     fixedEffectColumns=c("SCALED_LOG_UMIS", "pct_mt")
#     if (usePeerFactors) {
#         peerFactors=colnames(df)[grep ("PEER", colnames (df))]
#         log_info("Using [", length(peerFactors), "] PEER factors")
#         fixedEffectColumns=c(fixedEffectColumns, peerFactors)
#     }
#
#     X=as.matrix(df[,fixedEffectColumns,with=F])
#     X=cbind(intercept=1, X)
#
#     #index=23;
#     runOne<-function (index, bulkResults, X, groupData, genotypeDF) {
#         gene=bulkResults[index,]$gene
#         snp=bulkResults[index,]$SNP
#         log_info("Gaussian process genotype fit gene [", gene, "] SNP [", snp, "]")
#         y=as.numeric (expData[gene,])
#         g=as.numeric (genotypeDF[genotypeDF$id==snp,idxToGenotype, with=F])
#         XX=cbind(genotype=g, X)
#         #the null model doesn't include genotype
#         system.time(null_gp_model <- gpboost::fitGPModel(group_data=group_data, likelihood="poisson", y=y, X=X, params=list(optimizer_cov = "nelder_mead")))
#         system.time(gp_model <- gpboost::fitGPModel(group_data=group_data, likelihood="poisson", y=y, X=XX, params=list(optimizer_cov = "nelder_mead")))
#     }
#
# }

#Lots of testing code for the original version, with many different ways to do this.
#Sum works, median doesn't, and  Seurat is the speed winner for aggregation by orders of magnitude.
# reduceDataTesting<-function (numClustersPerDonor=NULL, featureName, featureList, metricsDF, expData, usageDF, reduceDataPDF=NULL) {
#     #Used to shortcut testing.
#     # d=prepareDataForSingleCellEqtl(snpFile, h5SeuratFile,
#     #                                cellMetaDataFile, donorCovariatesFile, permutedPseudobulkEqtlResultsFile=NULL,
#     #                                usageFile, qvalueThreshold=NULL, eQTLGeneList=NULL, chromosome=NULL,
#     #                                factorList=factorName)
#     #
#     # metricsDF=d$metricsDF
#     # expData=d$expData
#     # usageDF=d$usage
#     # featureList=getFeatureList (fixedFeatures=c("SCALED_LOG_UMIS", "pct_mt"), featureName, usePeerFactors, metricsDF)
#
#     #make sure the featureName is in the feature list.
#     if (!featureName %in% featureList)
#         featureList=c(featureName, featureList)
#
#     #validate all features are already in the metricsDF
#     missingFeatures=setdiff(featureList, colnames(metricsDF))
#     if (length(missingFeatures)>0)
#         stop ("Requested features missing from cell level metrics [", paste(missingFeatures, collapse=""), "]", sep="")
#
#     donorList=sort(unique(metricsDF$DONOR))
#     showPlots=F
#     if (!is.null(reduceDataPDF)) {
#         showPlots=T
#         pdf(reduceDataPDF)
#     }
#
#     log_info("Generating [", numClustersPerDonor, "] prototype cells per donor")
#
#     r<-lapply(donorList, reduceDataDonor, metricsDF, featureList, nclust=numClustersPerDonor, normalize=T, featureName=featureName, showPlots=showPlots)
#     if (!is.null(reduceDataPDF)) dev.off()
#
#     #rbind the results
#     summary=do.call(rbind, lapply(r, function (x) x$summary))
#     exemplarCells=do.call(rbind, lapply(r, function (x) x$exemplarCells))
#     clusterMap=do.call(rbind, lapply(r, function (x) x$clusterMap))
#
#     #expression is the median of the cells expression in the exemplar group's expression
#     #that value needs to be rounded if the median is between two values
#     #Median summarization doesn't work - many genes have a median expression of 0 for all donor prototype cells.
#     #What about a sum of the UMIs, recalculate the num_retained transcripts to be the sum,
#     #and regenerate the SCALED_LOG_UMIS.
#
#     # log_info("Calculating expression for exemplar cells")
#     # system.time (expDataNew<-calculateExemplarExpression(exemplarCells, expData, clusterMap, summaryFunction="sum"))
#     # log_info("Finished calculating expression for exemplar cells")
#
#     # log_info("Calculating expression for exemplar cells")
#     # system.time (expDataNew2<-calculateExemplarExpression(exemplarCells, expData, clusterMap, summaryFunction="sum_faster"))
#     # log_info("Finished calculating expression for exemplar cells")
#
#     # log_info("Calculating expression for exemplar cells")
#     # expDataNew3<-calculateExemplarExpression(exemplarCells, expData, clusterMap, summaryFunction="sum_vectorized")
#     # log_info("Finished calculating expression for exemplar cells")
#
#     log_info("Calculating expression for exemplar cells")
#     expDataNew4<-calculateExemplarExpression(exemplarCells, expData, clusterMap, summaryFunction="Seurat")
#     log_info("Finished calculating expression for exemplar cells")
#
#     # expDataNew4=expDataNew4[,match(colnames(expDataNew2), colnames(expDataNew4))]
#     #all (expDataNew2==expDataNew4)
#
#     #recalculate the num_retained_transcripts and SCALED_LOG_UMIS
#     if (summaryFunction=="sum") {
#         #validate total expression of each gene remains the same, but it's distributed differently.
#         if (rowSums(expData)!=rowSums(expDataNew))
#             stop ("Something went wrong")
#
#         exemplarCells=recalculateTotalTranscripts(exemplarCells, clusterMap, metricsDF)
#         exemplarCells$SCALED_LOG_UMIS=scale(log(exemplarCells$num_retained_transcripts))
#     }
#
#     result=list(metricsDF=exemplarCells, expData=expDataNew)
#     return (result)
# }


compare<-function () {

    d=prepareDataForSingleCellEqtl(snpFile, h5SeuratFile,
                                   cellMetaDataFile, donorCovariatesFile, permutedPseudobulkEqtlResultsFile,
                                   usageFile, qvalueThreshold=qvalueThreshold, eQTLGeneList=eQTLGeneList, chromosome=chromosome,
                                   factorList=factorName)

    metricsDF=d$metricsDF
    bulkResults=d$bulkResults
    genotypeDF=d$genotypeDF
    expData=d$expData
    usageDF=d$usage

    #compare metrics.
    oldMetricsDF=read.table("/downloads/single_cell_eQTL/old_code/metricsDF.txt", header=T, stringsAsFactors = F, sep="\t")
    oldMetricsDF=oldMetricsDF[match(metricsDF$CELL_BARCODE_FINAL, oldMetricsDF$CELL_BARCODE_FINAL),]
    cols=intersect (colnames(oldMetricsDF), colnames(metricsDF))
    data.table::setDF(metricsDF)

    #FACTOR_2 is clearly screwed up!
    cor (as.numeric (metricsDF$FACTOR_2), as.numeric (oldMetricsDF$FACTOR_2))


    oldExpression=read.table("/downloads/single_cell_eQTL/old_code/expression.txt", header=T, stringsAsFactors=F, sep="\t", check.names=F)
    oldExpression=Matrix(as.matrix(oldExpression), sparse = T)



}


#calculate the permuted pvalue in two ways
# 1- use only the empiric pvalues
# 2 - fit the permuted distribution to a beta distribution, use the CDF to calculate the probabiltiy of
# drawing the empiric result.
# getPermutedPvalue<-function (empiricP, permutedPDist) {
#
#     mean=mean(permutedPDist)
#     var=var(permutedPDist)
#     beta_shape1=mean*(mean*(1-mean) /var-1)
#     beta_shape2=beta_shape1*(1/mean-1)
#     p = seq(0, 1, length=100)
#     plot(p, dbeta(p, beta_shape1, beta_shape2), type='l')
#
#     f=fitdistrplus::fitdist(data = permutedPDist, "beta", start = list(shape1 = beta_shape1, shape2 = beta_shape2))
#
#     f=fitdistrplus::fitdist(data = permutedPDist, "beta", control=list(trace=1, REPORT=1))
#     pval_perm=(length(which(permutedPDist<=empiricP))+1)/(length(permutedPDist)+1)
#     pval_beta=pbeta(empiricP, f$estimate[1], f$estimate[2], lower.tail = TRUE)
#     data.frame(pval_perm=pval_perm, pval_beta=pval_beta)
# }
