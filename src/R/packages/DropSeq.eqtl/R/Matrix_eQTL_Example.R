# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2017 by the
# Broad Institute/Massachusetts Institute of Technology. All rights are
# reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for its
# use, misuse, or functionality.

# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# This example take from base example provided on website, then enhanced.

# library(MatrixEQTL);library (data.table);library (tictoc)
# source ("/Users/nemesh/dropseqrna/transcriptome/R/packages/DropSeq.eqtl/R/Matrix_eQTL_CommonFunctions.R")
# source ("/Users/nemesh/dropseqrna/transcriptome/R/packages/DropSeq.eqtl/R/AddGeneSNPPositions.R")
## Location of the package with the data files.
# base.dir="/downloads/eQTL/villageC/out_nextseq_novaseq"
# base.name="DropulationVillageC_all"

#' eQTL discovery.
#'
#' Runs Matrix eQTL package on a data set.
#'
#' This validates files exist, loads up data, runs the analysis, and writes out the results.
#' Optionally you can filter variants on MAF, or quantile normalize the expression data.
#' @seealso \code{\link{Matrix_eQTL_main}}.
#' @param params_file_name Matrix_eQTL_main param elements (except multi-valued ones) written here in tabular form.
#' @param cis_eqtl_file_name Results of Matrix_eQTL_main written here.
#' @param SNP_file_name A matrix of genotypes, 1 row per SNP, 1 column per sample.  Encoded as 0,1,2 copies of the alternate allele.
#' @param snps_location_file_name The location of each SNP.  Same number of lines as the SNP_file_name, encodes the snp name, chromosome, position.
#' @param expression_file_name A matrix of expression data, 1 row per gene, 1 column per sample.
#' @param gene_location_file_name The location of each gene.  Same number of lines as the expression_file_name, encodes the gene name, chromosome, start, end.
#' @param covariates_file_name (optional) A matrix of 1 or more covariates.  Each row is a covariate, each column a sample.
#' @param useQuantileNormalization Should the data be quantile normalized?
#' @param minorAlleleFreqThreshold Should the SNPs be filtered on minor allele frequency?
#' @param pvOutputThreshold_cis cis eQTLs must have a pvalue less than or equal to this threshold to be reported.
#' @param useModel The model MatrixEQTL uses.
#' @param cisDist defines the window around each gene to find associated SNPs.
#' @param verbose output detailed information about the eQTL process as it runs (passed to MatrixEQTL::Matrix_eQTL_main)
#' @param noFDRsaveMemory passed through to MatrixEQTL::Matrix_eQTL_main.
#' @param snp_gene_map_file_name optional file with columns 'snp' and 'gene' that represents a many-to-many
#'                               map.  If present, only the list of SNPs mapped to a gene are eQTL candidates.
#' @param errorCovarianceMatrixFile If supplied, a matrix of the relationships between donors, which is calculated as (kinship *2).
#' This can be a superset of the donors in the other files.
#' @return the eQTL result data frame.
#' @import data.table utils
#' @seealso MatrixEQTL::Matrix_eQTL_main
#' @export
runEQTL<-function (params_file_name, cis_eqtl_file_name, SNP_file_name, snps_location_file_name,
    expression_file_name, gene_location_file_name,
    covariates_file_name=NULL, useQuantileNormalization=F, minorAlleleFreqThreshold=NULL,
    pvOutputThreshold_cis=1, useModel=modelLINEAR, cisDist=100000, verbose=TRUE, noFDRsaveMemory = TRUE,
    snp_gene_map_file_name=NULL, errorCovarianceMatrixFile=NULL) {

    validateFilesExist(snps_location_file_name, SNP_file_name, expression_file_name, gene_location_file_name, covariates_file_name)
    validateSingleFileExists(snp_gene_map_file_name)
    validateSingleFileExists(errorCovarianceMatrixFile)
    if (is.null(snp_gene_map_file_name)) {
        snp_gene_map=NULL
    } else {
        snp_gene_map = DropSeq.utilities::fastRead(snp_gene_map_file_name, comment_regexp='^#', header = TRUE)
    }

    #load up the genotype data and optionally filter by MAF.
    r=loadGenotypeData(SNP_file_name, snps_location_file_name, minorAlleleFreqThreshold = minorAlleleFreqThreshold,
                       snp_gene_map=snp_gene_map)
    snps=r$snps
    snpspos=r$snpspos

    gene=loadExpressionData(expression_file_name, useQuantileNormalization = useQuantileNormalization)
    cvrt=loadCovariates(covariates_file_name)
    ## Run the analysis
    genepos = fread(gene_location_file_name, data.table=F);
    if (!is.null(snp_gene_map)) {
        # make the contig be the same as the gene
        genepos$chr = genepos$geneid
    }

    errorCovariance=getErrorCovariance(errorCovarianceMatrixFile)
    errorCovariance=validateAndFilterErrorCovarianceMatrix (errorCovariance, gene)

    system.time(me<- MatrixEQTL::Matrix_eQTL_main(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        pvOutputThreshold     = 0,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = verbose,
        output_file_name.cis = cis_eqtl_file_name,
        pvOutputThreshold.cis = pvOutputThreshold_cis,
        snpspos = filterSNPsPosForMatrix_eQTL(snpspos),
        genepos = genepos,
        cisDist = cisDist,
        pvalue.hist = FALSE,
        min.pv.by.genesnp = T,
        noFDRsaveMemory = noFDRsaveMemory))

    ## Results:
    cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
    params_df = do.call(cbind.data.frame, Filter(function(x) length(x) == 1, me$param))
    write.table(params_df, params_file_name, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
    #z=me$cis$eqtls
    #cat("Number of genes with FDR<=0.05", length(unique(z[z$FDR<0.05,]$gene)), "\n")
    #plot(me, main=basename(output_file_name_cis))
    invisible(me)
}

#' eQTL discovery.
#'
#' Annotates eQTLs produced by Matrix eQTL package.
#'
#' @seealso \code{\link{Matrix_eQTL_main}}.
#' @param output_file_name_cis Result written here.
#' @param cis_eqtl_file_name output produced by Matrix_eQTL_main(output_file_name.cis)
#' @param params_file_name Most Matrix_eQTL_main returned params, in tabular format.
#' @param SNP_file_name A matrix of genotypes, 1 row per SNP, 1 column per sample.  Encoded as 0,1,2 copies of the alternate allele.
#' @param snps_location_file_name The location of each SNP.  Same number of lines as the SNP_file_name, encodes the snp name, chromosome, position.
#' @param expression_file_name A matrix of expression data, 1 row per gene, 1 column per sample.
#' @param gene_location_file_name The location of each gene.  Same number of lines as the expression_file_name, encodes the gene name, chromosome, start, end.
#' @param useQuantileNormalization Should the data be quantile normalized?
#' @param minorAlleleFreqThreshold Should the SNPs be filtered on minor allele frequency?
#' @param snp_gene_map_file_name optional file with columns 'snp' and 'gene' that represents a many-to-many
#'                               map.  Note that currently this is only used to indicate that SNP names
#'                               have gene suffix that needs to be removed
#'
#' @import data.table utils
#' @seealso MatrixEQTL::Matrix_eQTL_main
#' @export
annotateEQTL<-function (output_file_name_cis, cis_eqtl_file_name, params_file_name, SNP_file_name,
    snps_location_file_name, expression_file_name, gene_location_file_name, useQuantileNormalization=F,
    minorAlleleFreqThreshold=NULL, snp_gene_map_file_name=NULL) {

    result=DropSeq.utilities::fastRead(cis_eqtl_file_name)

    if (!is.null(snp_gene_map_file_name)) {
        # SNP looks like this: chr22:17084998:C:A:IL17RA ; i.e. contig:position:A-allele:B-allele:gene
        # This function chops off the :gene
        result$SNP = undoPerGeneSnpName(result$SNP)
    }

    params = DropSeq.utilities::fastRead(params_file_name)
    result[,r2 := (result$`t-stat` / sqrt( params$dfFull + result$`t-stat`^2 ))^2]
    result[,beta_se :=result$beta/result$`t-stat`]

    #add on some annotations to help later analysis.
    #add snp and gene location, and distance from snp to gene.
    result=getDistanceToGene (result, snps_location_file_name, gene_location_file_name)

    #load up the genotype data and optionally filter by MAF.
    r=loadGenotypeData(SNP_file_name, snps_location_file_name, minorAlleleFreqThreshold = minorAlleleFreqThreshold)
    snps=r$snps

    #add the SNP AF.
    afAll=getAlleleFrequency(snps)
    idx=match(result$SNP, names(afAll))
    result[,MAF := as.numeric (round(as.numeric (afAll[idx]),3))]

    #add the average expression and the effect size normalized by expression.
    result=addMedianExpression(result, loadExpressionData(expression_file_name, useQuantileNormalization = useQuantileNormalization))

    #SNP	gene	beta	t-stat	p-value

    #format the output before we write it out to not have silly numbers of significant digits.
    #for R CMD CHECK
    "t-stat"=NULL;"p-value"=NULL;beta=NULL;r2=NULL;beta_se=NULL;MAF=NULL;median_expression=NULL;effect_size=NULL;
    result[,`t-stat`:=round(result$`t-stat`,5)]
    result[,`p-value`:=format(result$`p-value`, scientific = T, digits=5)]
    result[,beta:=round(result$beta,5)]
    result[,r2:=round(result$r2,5)]
    result[,beta_se:=round(result$beta_se,5)]
    result[,MAF:=round(result$MAF,5)]
    result[,median_expression:=round(result$median_expression,5)]
    result[,effect_size:=round(result$effect_size,5)]

    conn = DropSeq.utilities::open_conn(output_file_name_cis, open="w")
    write.table(result, conn, row.names=F, col.names=T, quote=F, sep="\t")
    close(conn)
}


#' eQTL discovery with assumed input file names.
#'
#' run eQTL analysis, using a convenience method to define the file name inputs/outputs to a standard format.
#' @seealso \code{\link{runEQTL}}.
#' @param base.dir The directory data resides in.
#' @param base.name The prefix of all files for this data set.
#' @inheritParams runEQTL
#' @return the eQTL result data frame.
#' @export
runEQTLFromBaseFileNames<-function (base.dir, base.name, useQuantileNormalization=F, minorAlleleFreqThreshold=NULL, pvOutputThreshold_cis=1, useModel=modelLINEAR, cisDist=100000) {
    # Genotype file name
    snps_location_file_name = getSNPLocationFileName(base.dir, base.name)
    SNP_file_name= getSNPFileName(base.dir, base.name)

    # Gene expression file name
    expression_file_name = getExpressionFileName(base.dir, base.name)
    gene_location_file_name = getGeneLocationFileName(base.dir, base.name)

    # Covariates file name
    covariates_file_name = getCovariateFileName(base.dir, base.name)
    # Output file name
    output_file_name_cis=getOutCiseQTLFileName(base.dir, base.name, useQuantileNormalization, cisDist)
    #useQuantileNormalization=F, minorAlleleFreqThreshold=NULL, pvOutputThreshold_cis=1, useModel=modelLINEAR, cisDist=100000
    runEQTL(SNP_file_name, snps_location_file_name, expression_file_name, gene_location_file_name, covariates_file_name, output_file_name_cis,
            useQuantileNormalization, minorAlleleFreqThreshold, pvOutputThreshold_cis, useModel, cisDist)

}

# Remove columns from snps_location data frame that are not expected by Matrix_eQTL
filterSNPsPosForMatrix_eQTL<-function(snpspos) {
    return(snpspos[,c('snp', 'chr', 'pos')])
}

#' Parse in the error covariance matrix if set, or return an empty numeric vector as expected by eQTL Matrix.
#' @param errorCovarianceMatrixFile A tab seperated square matrix of donor kinship scores (multipled by 2), with row and column names.
#' @return A square matrix of donor kinship scores if errorCovarianceMatrixFile is not null, otherwise a length 0 numeric vector
getErrorCovariance<-function (errorCovarianceMatrixFile=NULL) {
    if (is.null(errorCovarianceMatrixFile)) return (numeric())
    ecm=read.table(errorCovarianceMatrixFile, header=T, stringsAsFactors = F, sep="\t", check.names = F)
    ecm=as.matrix (ecm)
    return (ecm)
}

#' Validate that if the errorCovariance is not null, the donors are in the same order as the expression data.
#' If the errorCovariance donors are a superset of the donors in the genes data set, restrict the error covariance to be
#' the same set of donors.
#' If the errorCovariance does not include all donors in the gene object, but is supplied, this will generate an error
#' and analysis will not continue.
#' @param errorCovariance The error covariance, which is either length 0, or a matrix with column names that should match the expression data
#' @param gene A SlicedData object containing the expression data
#' @return The error errorCovariance matrix if it is valid.  If invalid, return NULL.
validateAndFilterErrorCovarianceMatrix<-function (errorCovariance, gene) {
    if (length(errorCovariance)==0) return (errorCovariance)
    if (all(colnames(gene)==colnames(errorCovariance)))
        return (errorCovariance)

    #otherwise, you need to do more work.
    donorsNeeded=colnames(gene)
    idx=match(donorsNeeded, colnames(errorCovariance))
    if (any(is.na(idx)))
        stop("Error Covariance Matrix present, but does not contain all donors in the expression data.  Can not proceed with analysis!")

    #otherwise, subset the error covariance matrix and return
    return (errorCovariance[idx,idx])
}


