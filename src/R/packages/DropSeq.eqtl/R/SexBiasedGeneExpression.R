# library (data.table); library(qvalue); library (MatrixEQTL);library (yaml); library (parallel); library (DropSeq.eqtl); library (qqman)
# source ("~/dropseqrna/transcriptome/R/packages/DropSeq.eqtl/R/Matrix_eQTL_CommonFunctions.R")

# base.dir="/downloads/eQTL/genderBias"
# base.name="d14-42_NGN2.maf_0.20_cisDist_10kb"
# contigGroupsFile="/downloads/h37_m38.contig_groups.yaml"

# base.dir="/downloads/eQTL/genderBias"
# base.name="d5QueenBasic.maf_0.20_cisDist_10kb"
# contigGroupsFile="/downloads/GRCh38_maskedAlt.contig_groups.yaml"

# base.dir="/downloads/eQTL/genderBias"
# base.name="d21STORMI.maf_0.20_cisDist_10kb"
# contigGroupsFile="/downloads/GRCh38_maskedAlt.contig_groups.yaml"

# covar=c("SEX", "FRACTION_X"); numPermutations=0; random.seed=1
# out.dir="/downloads/eQTL/genderBias/out"
# filterCovarName="SEX"; filterCovarValue=NULL;renormAutosomeExpression=T

#input files
# expression_file_name = getExpressionFileName(base.dir, base.name)
# gene_location_file_name = getGeneLocationFileName(base.dir, base.name)
# covariates_file_name = getCovariateFileName(base.dir, base.name)

#outfiles
# outPDF=paste(out.dir, "/", base.name, ".xx_biased_genes.pdf", sep="")
# outFile=paste(out.dir, "/", base.name, ".xx_biased_genes.txt", sep="")
# permutedResultsFile=outFile
# pValueThreshold=1e-5

#' Find genes with biased expression based SEX or FRACTION_X
#'
#' Converts FRACTION_X into X reactivation.  This normalizes fraction_X by the mean FRACTION_X of individuals with SEX=1.
#' This makes all SEX=1 individuals close to 1.
#'
#' Analysis is broken down into a few steps:
#' (Optionally) Normalize autosome expression to expression/sum(expression) [this removes the X chromosome influence on normalization]
#' 1) For the subset of XX individuals, regression : expression ~ X reactivation
#' 1a) Run permutation testing by permuting expression labels
#' 2) Normalize XX individuals expression to account for X reactivation covariate.
#' 3) Test the effects of SEX on normalized expression (expression [conditional on X reactivation] ~ SEX)
#' 3a) Run permutations on SEX regression
#'
#' @param expression_file_name A matrix of expression data, 1 row per gene, 1 column per sample.
#' @param gene_location_file_name The location of each gene.  Same number of lines as the expression_file_name, encodes the gene name, chromosome, start, end.
#' @param covariates_file_name A matrix of 1 or more covariates.  Each row is a covariate, each column a sample.  Contains a row with the donors' gender.
#' @param numPermutations How many permutations should be run to generate permuted p-values
#' @param outFile The report of each gene tested and the p-values and effect sizes.
#' @param contigGroupsFile A YAML file containing each contig name and groups it belongs to.  In this case the contigs should belong to the groups autosome and chrX.  If provided, renormalizes autosome expression data to not include the x chromosomne.
#' @param renormAutosomeExpression If a contigGroupsFile and this is true, then renormalize expression of the autosome to the sum of the autosome genes (instead of all genes)
#' @param random.seed A random seed to ensure reproducibility across runs of the same data.
#' @import data.table qvalue utils
#' @export
sexBiasedGeneExpression<-function (expression_file_name, covariates_file_name, gene_location_file_name=NULL,
                                   contigGroupsFile=NULL, outFile=NULL,
                                   numPermutations=10000, renormAutosomeExpression=T,
                                   random.seed=1234) {

    set.seed(random.seed)

    validateSingleFileExists(expression_file_name)
    validateSingleFileExists(covariates_file_name)

    covariatesToProcess=c("SEX", "FRACTION_X")
    expData=read.table(expression_file_name, header=T, stringsAsFactors = F, sep="\t", check.names = F)
    covars=read.table(covariates_file_name, header=T, stringsAsFactors = F, sep="\t", check.names = F)
    geneLoc=read.table(gene_location_file_name, header=T, stringsAsFactors = F, sep="\t", check.names = F)

    missingCovars=setdiff(covariatesToProcess, covars$id)
    if (length(missingCovars)>0) {
        warning("Missing required covariates [", paste(missingCovars, collapse=","), "]", immediate.=TRUE)
        return()
    }
    covars=convertToXReactivationPhenotype(covars)

    if (!is.null(contigGroupsFile) & renormAutosomeExpression)
        expData=DropSeq.eqtl::renormalizeAutosomeExpression(expData, geneLoc, contigGroupsFile)

    #normalize sex expression so that X reactivation = 1.
    expDataXReactivationNormalized=removeXReactivationEffectOnExpressionAllDonors(expData, covars)

    #get the upper and lower bounds of the X reactivation regression
    fitDF=calculateConfIntervalsXReactivationEffectOnExpressionAllDonors(expData, covars, confInterval=0.95)

    #the lower and uppber bounds to the data.
    expLow=applyXReactivationToExpression(expData, covars, fitDF$min_slope)
    expHigh=applyXReactivationToExpression(expData, covars, fitDF$max_slope)

    #run the covariate, and reorder results to match the expData ordering.
    testSex<-function (covarName="SEX", covars, expData, numPermutations=numPermutations) {
        sexResult=DropSeq.eqtl::runOneCovarMatrixeQTL (covarName="SEX", covars, expData=expData, numPermutations=numPermutations)
        idx=match(expData$id, sexResult$gene)
        return (sexResult[idx,])
    }

    sexResult=testSex (covarName="SEX", covars, expData=expDataXReactivationNormalized, numPermutations=numPermutations)
    sexResultNoNorm=testSex (covarName="SEX", covars, expData=expData, numPermutations=numPermutations)
    sexResultLow=testSex (covarName="SEX", covars, expData=expLow, numPermutations=numPermutations)
    sexResultHigh=testSex (covarName="SEX", covars, expData=expHigh, numPermutations=numPermutations)

    #in the same order as the expData.
    sexResultNonParametric = testSexNonParametric(expDataXReactivationNormalized, covars, covarName="SEX")
    sexResultNonParametricNoNorm = testSexNonParametric(expData, covars, covarName="SEX")


    #x reactivation and permutation
    #note that the BETA from X reactivation is the same as the beta from the fitDF!
    #all.equal(fitDF$slope, xReactivationResult$beta) #this is TRUE.
    xReactivationResult=testXReactivation(expData, covars, covarName="FRACTION_X", numPermutations=numPermutations)

    #order X reactivation results
    xReactivationResult=xReactivationResult[match(sexResult$gene, xReactivationResult$gene),]

    if (length(which(sexResult$gene!=xReactivationResult$gene))>0)
        warning("Problem with gene names matching up.", immediate.=TRUE)


    #bind it all together.
    #chr	start	end	gene	SEX_pvalue	FRACTION_X_pvalue	FULL_MODEL	SEX_beta	FRACTION_X_beta SEX_permuted_p	FRACTION_X_permuted_p	FULL_MODEL_permuted_p	SEX_qval	FRACTION_X_qval	FULL_MODEL_qval	median_expression
    result=data.frame(gene=sexResult$gene, SEX_pvalue=sexResult$`p-value`, SEX_beta=sexResult$beta,
                  SEX_raw_pvalue=sexResultNoNorm$`p-value`, SEX_raw_beta=sexResultNoNorm$beta,
                  SEX_low_pvalue=sexResultLow$`p-value`, SEX_low_beta=sexResultLow$beta,
                  SEX_high_pvalue=sexResultHigh$`p-value`, SEX_high_beta=sexResultHigh$beta,
                  SEX_raw_nonParametric=sexResultNonParametric, SEX_nonParametric=sexResultNonParametricNoNorm,
                  FRACTION_X_pvalue=xReactivationResult$`p-value`, FRACTION_X_beta=xReactivationResult$beta,
                  FRACTION_X_beta_min=fitDF$min_slope, FRACTION_X_beta_max=fitDF$max_slope,
                  stringsAsFactors = F)

    #this is a bit of a hack until we figure out the best way to solve the lack of variance in X reactivation
    #in some experiments leading to a poor estimate of the regression slope.
    idxRawBest=which(result$SEX_pvalue<result$SEX_raw_pvalue)
    result$Sex_conservative_pvalue=result$SEX_pvalue
    result[idxRawBest,]$Sex_conservative_pvalue=result[idxRawBest,]$SEX_raw_pvalue
    result$Sex_conservative_beta=result$SEX_beta
    result[idxRawBest,]$Sex_conservative_beta=result[idxRawBest,]$SEX_raw_beta

    #same with nonparametric result
    idxRawBest=which(result$SEX_nonParametric<result$SEX_raw_nonParametric)
    result$Sex_conservative_pvalue_nonparametric=result$SEX_nonParametric
    result[idxRawBest,]$Sex_conservative_pvalue_nonparametric=result[idxRawBest,]$SEX_raw_nonParametric

    #result$Sex_conservative_pvalue=pmax(result$SEX_pvalue, result$SEX_raw_pvalue)

    result=addMedianExpression(result, expData)

    #sort by pvalue.
    #result=result[order(result$FULL_MODEL_permuted_p),]

    idx=match(result$gene, geneLoc$geneid)
    loc=data.frame(chr=geneLoc[idx,]$chr, start=geneLoc[idx,]$s1, end=geneLoc[idx,]$s2)

    result=data.frame(loc, result)

    #formatting
    colsSci=c("SEX_pvalue", "SEX_raw_pvalue", "SEX_low_pvalue",  "SEX_high_pvalue",
           "FRACTION_X_pvalue")

    colsOther=c("SEX_beta", "SEX_raw_beta", "SEX_low_beta", "SEX_high_beta",
                "FRACTION_X_beta", "FRACTION_X_beta_min", "FRACTION_X_beta_max", "median_expression")

    for (col in colsSci)
        result[[col]]=format(result[[col]], scientific=T, digits=5)

    for (col in colsOther)
        result[[col]]=round(result[[col]], 5)

    write.table(result, outFile, row.names=F, col.names = T, quote=F, sep="\t")

}



#' Plots the SEX/FRACTION_X regression results.
#'
#' @param expression_file_name A matrix of expression data, 1 row per gene, 1 column per sample.
#' @param covariates_file_name A matrix of 1 or more covariates.  Each row is a covariate, each column a sample.  Contains a row with the donors' gender.
#' @param gene_location_file_name The location of each gene.  Same number of lines as the expression_file_name, encodes the gene name, chromosome, start, end.
#' @param permutedResultsFile The report of each gene tested and the p-values and effect sizes for the covariates SEX and FRACTION_X
#' @param contigGroupsFile A YAML file containing each contig name and groups it belongs to.  In this case the contigs should belong to the groups autosome and chrX.  If provided, renormalizes autosome expression data to not include the x chromosomne.
#' @param renormAutosomeExpression If a contigGroupsFile and this is true, then renormalize expression of the autosome to the sum of the autosome genes (instead of all genes)
#' @param pValueThreshold Plot all genes that have at a pvalue < this score in one category.
#' @param outPDF Output file for plots.
#' @export
plotGenesXXDosage<-function (expression_file_name, covariates_file_name, gene_location_file_name, permutedResultsFile, contigGroupsFile, renormAutosomeExpression=T, pValueThreshold=5e-08, outPDF=NULL) {

    validateSingleFileExists(covariates_file_name)
    #finish the covariate validation first - this can short circuit the analysis.
    covars=read.table(covariates_file_name, header=T, stringsAsFactors = F, sep="\t", check.names = F)
    covariatesToProcess=c("SEX", "FRACTION_X")

    missingCovars=setdiff(covariatesToProcess, covars$id)
    if (length(missingCovars)>0) {
        warning("Missing required covariates [", paste(missingCovars, collapse=","), "]", immediate.=TRUE)
        return()
    }

    #if the analysis was valid, then validate the rest of the input files.
    validateSingleFileExists(expression_file_name)
    validateSingleFileExists(gene_location_file_name)
    validateSingleFileExists(permutedResultsFile)

    expData=read.table(expression_file_name, header=T, stringsAsFactors = F, sep="\t", check.names = F)
    result=read.table(permutedResultsFile, header=T, stringsAsFactors=F, check.names = F)
    geneLoc=read.table(gene_location_file_name, header=T, stringsAsFactors = F, sep="\t", check.names = F)

    covars=convertToXReactivationPhenotype(covars)

    if (!is.null(contigGroupsFile) & renormAutosomeExpression) {
        expData=DropSeq.eqtl::renormalizeAutosomeExpression(expData, geneLoc, contigGroupsFile)
    }


    #expDataXReactivationNormalized=removeXReactivationEffectOnExpressionXXDonors(expData, covars)
    expDataXReactivationNormalized=removeXReactivationEffectOnExpressionAllDonors(expData, covars)

    if (!is.null(outPDF)) pdf(outPDF, width=10, height = 7)
    #plot summary stats
    plotXReactivationVsSex(covars)
    resultGeneList=plotGenesXXDosageSummaryMetrics(result, contigGroupsFile, pValueThreshold, outPDF=NULL)

    z=sapply(resultGeneList$resultAutosomeSex, plotExampleGeneXXDosage, expData, expDataXReactivationNormalized, covars, result, "X Dosage")
    z=sapply(resultGeneList$resultXSex, plotExampleGeneXXDosage, expData, expDataXReactivationNormalized, covars, result, "X Dosage")
    z=sapply(resultGeneList$resultAutosomeBoth, plotExampleGeneXXDosage, expData, expDataXReactivationNormalized, covars, result, "Both")
    z=sapply(resultGeneList$resultXBoth, plotExampleGeneXXDosage, expData, expDataXReactivationNormalized, covars, result, "Both")
    z=sapply(resultGeneList$resultAutosomeFrac, plotExampleGeneXXDosage, expData, expDataXReactivationNormalized, covars, result, "Reactivation")
    z=sapply(resultGeneList$resultXFrac, plotExampleGeneXXDosage, expData, expDataXReactivationNormalized, covars, result, "Reactivation")
    #plotExampleGeneXXDosage("GPX3", expData, expDataXReactivationNormalized, covars, result, "")
    #z=sapply(sample (result[result$SEX_pvalue>0.1 & result$FRACTION_X_pvalue>0.1 & result$median_expression>5,]$gene, 20), plotExampleGeneXXDosage, expData, expDataXReactivationNormalized, covars, result, "Reactivation")
    if (!is.null(outPDF)) dev.off()

}

#For autosomal genes (1 for all autosomal genes, 1 for genes correlated with X reactivation)
# Plot variance of gene expression of XX vs XY (scatter plot with XX on the X axis)
# Plot variance of gene expression of XX vs XY after correcting for the X reactivation (scatter plot with XX on the X axis)
# Is there more variance in the XX or the XY lines in each plot?  (Maybe plot the densities of variance for each population XX/XY?)
# After the correction, is there more variance in the XX XY line?

plotGeneExpressionVarianceBySex<-function (expression_file_name, covariates_file_name, gene_location_file_name, permutedResultsFile) {
    validateSingleFileExists(expression_file_name)
    validateSingleFileExists(covariates_file_name)
    validateSingleFileExists(gene_location_file_name)
    validateSingleFileExists(permutedResultsFile)

    expData=read.table(expression_file_name, header=T, stringsAsFactors = F, sep="\t", check.names = F)
    covars=read.table(covariates_file_name, header=T, stringsAsFactors = F, sep="\t", check.names = F)
    result=read.table(permutedResultsFile, header=T, stringsAsFactors=F, check.names = F)
    geneLocs=read.table(gene_location_file_name, header=T, stringsAsFactors = F, sep="\t", check.names = F)

    covarNames=sub("_pvalue", "", colnames(result)[grep("_pvalue", colnames(result))])
    expectedCovarNames=c("SEX", "FRACTION_X")
    missingCovars=setdiff(expectedCovarNames, covarNames)
    if (length(missingCovars)>0)
        stop(paste("Must have the expected covariants to use this plot", paste(expectedCovarNames, collapse=" ")))

    covars=convertToXReactivationPhenotype(covars)
    #expDataXReactivationNormalized=removeXReactivationEffectOnExpressionXXDonors(expData, covars)
    expDataXReactivationNormalized=removeXReactivationEffectOnExpressionAllDonors(expData, covars)

    sex=as.numeric (covars[covars$id=="SEX",-1])
    x_reactivation=as.numeric (covars[covars$id=="FRACTION_X",-1])

    getExpression<-function (s=1, sex, expData, geneLocs, geneList, xChromosomeOnly=F) {
        idxX=match(geneLocs[geneLocs$chr=="X",]$geneid, expData$id)
        if (xChromosomeOnly) {
            z=expData[idxX,]
        } else {
            z=expData[-idxX,]
        }

        if (!is.null(geneList)) {
            idx=match(intersect (geneList, z$id), z$id)
            z=z[idx,]
        }
        rownames(z)=z$id
        z=z[,-1]
        z=z[,which(sex==s)]
        return (z)
    }

    getExpressionVariance<-function (s=1, sex, expData, geneLocs, geneList=NULL, xChromosomeOnly=F, calculateCOV=F) {
        z=getExpression(s, sex, expData, geneLocs, geneList, xChromosomeOnly)
        if (calculateCOV) {
            v=apply(z, 1, var)
            m=apply(z, 1, mean)
            r=v/m
        } else {
            r=apply(z, 1, var)
        }

        return (r)
    }

    plotOne<-function (xxVar, xyVar, strTitle="", xyLim=NULL, ...) {
        if (is.null(xyLim)) xyLim=c(floor (min (log10(c(xxVar, xyVar)))), ceiling (max (log10(c(xxVar, xyVar)))))
        strTitle=paste(strTitle, "median var XX=", round(median(xxVar),3), "median var XY=", round(median (xyVar),3))
        smoothScatter(log10(xxVar), log10(xyVar), xlim=xyLim, ylim=xyLim, main=strTitle, cex.axis=1.25, cex.lab=1.25, ...)
        abline (0,1, lty=2, lwd=3, col='red')
    }

    #what about using the genes highly correlated to X reactivation?
    corToReactivation=getCorrelationToXReactivation(sex, expData, geneLocs)
    lowerBound=quantile (corToReactivation, probs=seq(0, 1, 0.05))[["10%"]]
    upperBound=quantile (corToReactivation, probs=seq(0, 1, 0.05))[["90%"]]
    hist (corToReactivation, breaks=100, main="Gene expression correlation to X reactivation [XX Donors]", xlab="expression correlation to X reactivation", cex.lab=1.25, cex.axis=1.25)
    abline (v=c(lowerBound, upperBound), col="red")


    xyLim=c(-6, 6)
    xxVarN=getExpressionVariance(s=2, sex, expData, geneLocs, NULL, T)
    xyVarN=getExpressionVariance(s=1, sex, expData, geneLocs, NULL, T)
    plotOne(xxVarN, xyVarN, "X chromosome genes only", xlab="XX expression variance [log10]", ylab="XY expression variance [log10]", xyLim=xyLim)

    xxVarN=getExpressionVariance(s=2, sex, expDataXReactivationNormalized, geneLocs, NULL, T)
    xyVarN=getExpressionVariance(s=1, sex, expDataXReactivationNormalized, geneLocs, NULL, T)
    plotOne(xxVarN, xyVarN, "X chromosome genes only", xlab="XX expression variance [normalized] [log10]", ylab="XY expression variance [log10]", xyLim=xyLim)


    xxVarN=getExpressionVariance(s=2, sex, expData, geneLocs, NULL, T, calculateCOV=T)
    xyVarN=getExpressionVariance(s=1, sex, expData, geneLocs, NULL, T, calculateCOV=T)
    plotOne(xxVarN, xyVarN, "X chromosome genes only", xlab="XX expression COV [log10]", ylab="XY expression variance [log10]", xyLim=xyLim)

    xxVarN=getExpressionVariance(s=2, sex, expDataXReactivationNormalized, geneLocs, NULL, calculateCOV=T)
    xyVarN=getExpressionVariance(s=1, sex, expDataXReactivationNormalized, geneLocs, NULL, calculateCOV=T)
    plotOne(xxVarN, xyVarN, "X chromosome genes only", xlab="XX expression COV [normalized] [log10]", ylab="XY expression variance [log10]", xyLim=xyLim)

    xxVar=getExpressionVariance(s=2, sex, expData, geneLocs)
    xyVar=getExpressionVariance(s=1, sex, expData, geneLocs)
    plotOne(xxVar, xyVar, "All autosome genes", xlab="XX expression variance [log10]", ylab="XY expression variance [log10]", xyLim=xyLim)

    xxVarN=getExpressionVariance(s=2, sex, expDataXReactivationNormalized, geneLocs)
    xyVarN=getExpressionVariance(s=1, sex, expDataXReactivationNormalized, geneLocs)
    plotOne(xxVarN, xyVarN, "All autosome genes", xlab="XX expression variance [normalized] [log10]", ylab="XY expression variance [log10]", xyLim=xyLim)

    xxVar=getExpressionVariance(s=2, sex, expData, geneLocs, calculateCOV=T)
    xyVar=getExpressionVariance(s=1, sex, expData, geneLocs, calculateCOV=T)
    plotOne(xxVar, xyVar, "All autosome genes", xlab="XX expression COV [log10]", ylab="XY expression variance [log10]", xyLim=xyLim)

    xxVarN=getExpressionVariance(s=2, sex, expDataXReactivationNormalized, geneLocs, calculateCOV=T)
    xyVarN=getExpressionVariance(s=1, sex, expDataXReactivationNormalized, geneLocs, calculateCOV=T)
    plotOne(xxVarN, xyVarN, xyLim=xyLim, "All autosome genes", xlab="XX expression COV [normalized] [log10]", ylab="XY expression COV [log10]")



    #BOTTOM 10%
    xxVar=getExpressionVariance(s=2, sex, expData, geneLocs, geneList=names(corToReactivation[corToReactivation<=lowerBound]))
    xyVar=getExpressionVariance(s=1, sex, expData, geneLocs, geneList=names(corToReactivation[corToReactivation<=lowerBound]))
    plotOne(xxVar, xyVar, xyLim=xyLim, "Genes correlated to X reactivation [bottom 10% negative correlation]\n", xlab="XX expression variance [log10]", ylab="XY expression variance [log10]")

    xxVar=getExpressionVariance(s=2, sex, expDataXReactivationNormalized, geneLocs, geneList=names(corToReactivation[corToReactivation<=lowerBound]))
    xyVar=getExpressionVariance(s=1, sex, expDataXReactivationNormalized, geneLocs, geneList=names(corToReactivation[corToReactivation<=lowerBound]))
    plotOne(xxVar, xyVar, xyLim=xyLim, "Genes correlated to X reactivation [bottom 10% negative correlation]\n", xlab="XX expression variance X reactivation normalized [log10]", ylab="XY expression variance [log10]")

    #BOTTOM 10% COV
    xxVar=getExpressionVariance(s=2, sex, expData, geneLocs, geneList=names(corToReactivation[corToReactivation<=lowerBound]), calculateCOV=T)
    xyVar=getExpressionVariance(s=1, sex, expData, geneLocs, geneList=names(corToReactivation[corToReactivation<=lowerBound]), calculateCOV=T)
    plotOne(xxVar, xyVar, xyLim=xyLim, "Genes correlated to X reactivation [bottom 10% negative correlation]\n", xlab="XX expression COV [log10]", ylab="XY expression variance [log10]")

    xxVar=getExpressionVariance(s=2, sex, expDataXReactivationNormalized, geneLocs, geneList=names(corToReactivation[corToReactivation<=lowerBound]), calculateCOV=T)
    xyVar=getExpressionVariance(s=1, sex, expDataXReactivationNormalized, geneLocs, geneList=names(corToReactivation[corToReactivation<=lowerBound]), calculateCOV=T)
    plotOne(xxVar, xyVar, xyLim=xyLim, "Genes correlated to X reactivation [bottom 10% negative correlation]\n", xlab="XX expression COV X reactivation normalized [log10]", ylab="XY expression variance [log10]")


    #TOP 10%
    xxVar=getExpressionVariance(s=2, sex, expData, geneLocs, geneList=names(corToReactivation[corToReactivation>=upperBound]))
    xyVar=getExpressionVariance(s=1, sex, expData, geneLocs, geneList=names(corToReactivation[corToReactivation>=upperBound]))
    plotOne(xxVar, xyVar, xyLim=xyLim, "Genes correlated to X reactivation [top 10% positive correlation]\n", xlab="XX expression variance [log10]", ylab="XY expression variance [log10]")

    xxVar=getExpressionVariance(s=2, sex, expDataXReactivationNormalized, geneLocs, geneList=names(corToReactivation[corToReactivation>=upperBound]))
    xyVar=getExpressionVariance(s=1, sex, expDataXReactivationNormalized, geneLocs, geneList=names(corToReactivation[corToReactivation>=upperBound]))
    plotOne(xxVar, xyVar, xyLim=xyLim, "Genes correlated to X reactivation [top 10% positive correlation]\n", xlab="XX expression variance normalized [log10]", ylab="XY expression variance [log10]")

    #TOP 10% COV
    xxVar=getExpressionVariance(s=2, sex, expData, geneLocs, geneList=names(corToReactivation[corToReactivation>=upperBound]),calculateCOV=T)
    xyVar=getExpressionVariance(s=1, sex, expData, geneLocs, geneList=names(corToReactivation[corToReactivation>=upperBound]), calculateCOV=T)
    plotOne(xxVar, xyVar, xyLim=xyLim, "Genes correlated to X reactivation [top 10% positive correlation]\n", xlab="XX expression COV [log10]", ylab="XY expression variance [log10]")

    xxVar=getExpressionVariance(s=2, sex, expDataXReactivationNormalized, geneLocs, geneList=names(corToReactivation[corToReactivation>=upperBound]))
    xyVar=getExpressionVariance(s=1, sex, expDataXReactivationNormalized, geneLocs, geneList=names(corToReactivation[corToReactivation>=upperBound]))
    plotOne(xxVar, xyVar, xyLim=xyLim, "Genes correlated to X reactivation [top 10% positive correlation]\n", xlab="XX expression COV normalized [log10]", ylab="XY expression variance [log10]")

    #take the genes highly correlated to XX reactivation, and plot their corr
    geneList=names(corToReactivation[corToReactivation>=upperBound])
    XYExpressionUpper=getExpression(s=1, sex, expData, geneLocs, geneList=geneList, xChromosomeOnly=F)
    XXExpressionUpper=getExpression(s=2, sex, expData, geneLocs, geneList=geneList, xChromosomeOnly=F)
    #XYExpression=getExpression(s=1, sex, expData, geneLocs, geneList=NULL, xChromosomeOnly=F)

    #control plot genes correlated to X reactivation on XX and XY individuals.

    #lets try plotting the expression of the genes highly correlated to X reactivation in XY donors.
    #heatmap uses the distance function dist and hclust : hclust(dist(1-cor))
    #library (heatmap3)
    heatmap3::heatmap3 (cor(t(XYExpressionUpper), use = "pa"), method="median", symm=T, main="XY individuals positively correlated to X reactivation", showColDendro=F, labRow=FALSE, labCol=FALSE, margins=c(3,3))
    heatmap3::heatmap3 (cor(t(XXExpressionUpper), use = "pa"), method="median", symm=T, main="XX individuals positively correlated to X reactivation", showColDendro=F, labRow=FALSE, labCol=FALSE, margins=c(3,3))

    xxMedianCor=median (cor(t(XXExpressionUpper)))
    xyMedianCor=median (cor(t(XYExpressionUpper)))

    #sample same number of genes as XXExpressionUpper from the autosomes and calculate correlation.
    xy=getExpression(s=1, sex, expData, geneLocs, geneList=NULL, xChromosomeOnly=F)
    getPerm<-function (xy, numGenes=dim (XXExpressionUpper)[1], plot=F) {
        idx=sort(sample(1:dim(xy)[1], numGenes))
        xyp=xy[idx,]
        if (plot) {
            heatmap3::heatmap3 (cor(t(xyp), use = "pa"), method="median", symm=T, main="XY individuals random genes", showColDendro=F, labRow=FALSE, labCol=FALSE, margins=c(3,3))
        }
        return (median(cor(t(xyp))))
    }

    permMedian=replicate(1000, getPerm(xy, numGenes=dim (XXExpressionUpper)[1], plot=F))
    hist (permMedian, breaks=100, main="Median correlation of random gene set", xlab="median correlation", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
    abline (v=xyMedianCor, col='red', lwd=3)
    legend("top", legend=c("random genes", "top 10% X genes correlated to X reactivation"), fill=c("black", "red"), cex=1.25)

}


plotXReactivationVsSex<-function (covars, strTitle=NULL) {
    s=as.numeric (covars[covars$id=="SEX",-1])
    xr=as.numeric (covars[covars$id=="FRACTION_X",-1])
    ylim=c(min(xr)*0.9, 2)
    cols=c("blue", "red")
    plot (jitter(s, factor=0.2), xr, xlim=c(0.75, 2.25), xlab="XX/XY", ylab="Fraction of mRNAs from chrX genes (normalized)", axes=F, ylim=ylim, col=cols[s], cex.lab=1.25, cex.axis=1.25, pch=16, cex=1.25)
    axis(1, at=c(1,2), labels=c("XY", 'XX'))
    axis(2)
    if (is.null(strTitle)) strTitle= "X reactivation vs XX/XY status"
    title(strTitle)

}

#get the correlation between expression data and x reactivation for XX donors on autosomal genes
getCorrelationToXReactivation<-function (sex, expData, x_reactivation, geneLocs, xChromosomeString="X") {
    idx=match(geneLocs[geneLocs$chr!=xChromosomeString,]$geneid, expData$id)
    z=expData[idx,-1]
    z=z[,which(sex==2)]
    zz=as.vector(cor (x_reactivation[sex==2], t(z)))
    names (zz)=expData[idx,]$id
    return (zz)
}

plotGenesXXDosageSummaryMetrics<-function (result, contigGroupsFile, pValueThreshold=5e-08, outPDF=NULL) {
    if (!is.null(outPDF)) pdf(outPDF)

    z=splitResultsByContigGroup (result, contigColumn="chr", contigGroupsFile=contigGroupsFile)
    resultAutosome=z$autosomes
    resultX=z$chrX

    #sexPvalueColumn="Sex_bounded_pvalue"
    sexPvalueColumn="Sex_conservative_pvalue"

    #XX/XY biased but not x reactivation
    resultAutosomeSex=resultAutosome[which(resultAutosome[[sexPvalueColumn]]<=pValueThreshold & resultAutosome$FRACTION_X_pvalue>pValueThreshold),]
    resultXSex=resultX[which(resultX[[sexPvalueColumn]]<=pValueThreshold & resultX$FRACTION_X_pvalue>pValueThreshold),]

    #x reactivation but not XX/XY biased
    resultAutosomeFrac=resultAutosome[which(resultAutosome[[sexPvalueColumn]]>pValueThreshold & resultAutosome$FRACTION_X_pvalue<pValueThreshold),]
    resultXFrac=resultX[which(resultX[[sexPvalueColumn]]>pValueThreshold & resultX$FRACTION_X_pvalue<=pValueThreshold),]

    #BOTH
    resultAutosomeBoth=resultAutosome[which(resultAutosome[[sexPvalueColumn]]<=pValueThreshold & resultAutosome$FRACTION_X_pvalue<=pValueThreshold),]
    resultXBoth=resultX[which(resultX[[sexPvalueColumn]]<=pValueThreshold & resultX$FRACTION_X_pvalue<=pValueThreshold),]

    resultGeneList=list(resultAutosomeSex=resultAutosomeSex$gene, resultXSex=resultXSex$gene,
                        resultAutosomeBoth=resultAutosomeBoth$gene, resultXBoth=resultXBoth$gene,
                        resultAutosomeFrac=resultAutosomeFrac$gene, resultXFrac=resultXFrac$gene)

    summaryCounts=sapply(resultGeneList, length)

    ylim=c(0, max(ceiling(log10(summaryCounts)))+1)
    #fix for all summary counts being 0
    ylim[is.infinite(ylim)]<-1
    z=barplot (log10(as.numeric(summaryCounts+1)), names.arg=c("X Dosage", "X Dosage", "Both", "Both", "X Reactivation", "X Reactivation"), ylab="counts [log10]", col=c("light blue", "light green"), main="Summary of genes with FDR <=0.05", ylim=ylim)
    avg=mean(range(log10(summaryCounts+1)))
    #fix for all summary counts being 0
    if (avg<1) avg=0.5
    text(z[,1], rep(avg, dim(z)[1]), labels=as.numeric(summaryCounts), cex=1.5)
    legend("topleft", legend=c("Autosome", "X"), fill=c("light blue", "light green"), title="Chromosome")

    ylim=c(0, ceiling(max(-log10(c(result$SEX_pvalue,result$FRACTION_X_pvalue)))))
    manhattenPlotResults (genes=result$gene, pvalColumn="SEX_pvalue", result, contigGroupsFile, strTitle="Donor Sex [XX/XY]", ylim=ylim)
    manhattenPlotResults (genes=result$gene, pvalColumn="SEX_raw_pvalue", result, contigGroupsFile, strTitle="Donor Sex [not normalized] [XX/XY]", ylim=ylim)
    manhattenPlotResults (genes=result$gene, pvalColumn="Sex_conservative_pvalue", result, contigGroupsFile, strTitle="Donor Sex [conservative] [XX/XY]", ylim=ylim)

    manhattenPlotResults (genes=result$gene, pvalColumn="FRACTION_X_pvalue", result, contigGroupsFile, strTitle="X Reactivation", ylim=ylim)

    qqPlot(result$SEX_pvalue, main="Donor Sex [XX/XY] qqplot")
    qqPlot(resultAutosome$SEX_pvalue, main="Donor Sex [XX/XY] Autosomes qqplot")
    qqPlot(result$SEX_raw_pvalue, main="Donor Sex [not normalized] [XX/XY] qqplot")
    qqPlot(resultAutosome$SEX_pvalue, main="Donor Sex [not normalized] [XX/XY] Autosomes qqplot")
    qqPlot(result[[sexPvalueColumn]], main="Donor Sex [conservative] [XX/XY] qqplot")
    qqPlot(resultAutosome[[sexPvalueColumn]], main="Donor Sex [conservative] [XX/XY] Autosomes qqplot")
    qqPlot(resultAutosome$Sex_conservative_pvalue_nonparametric, main="Donor Sex [conservative nonparametric] [XX/XY] Autosomes qqplot")

    xylim=c(0, ceiling(max (-log10(resultAutosome$Sex_conservative_pvalue), -log10(resultAutosome$Sex_conservative_pvalue_nonparametric))))
    smoothScatter(-log10(resultAutosome$Sex_conservative_pvalue), -log10(resultAutosome$Sex_conservative_pvalue_nonparametric), xlab="Sex linear regression pvalue [-log10]", ylab="Sex Wilcox Test pvalue [-log10]", xlim=xylim, ylim=xylim, cex.lab=1.5, cex.axis=1.5)
    abline (0,1, col='red', lty=2)

    qqPlot(result$FRACTION_X_pvalue, main="X Reactivation qqplot")
    qqPlot(resultAutosome$FRACTION_X_pvalue, main="X Reactivation Autosomes qqplot")

    foldChange=resultAutosome$Sex_conservative_beta/resultAutosome$median_expression
    hist (foldChange, main="Fold change [beta/median expression] [conservative] autosome genes", xlab="Fold Change of sex biased expression", breaks=1000, xlim=c(-5,5))
    smoothScatter(foldChange, -log10(resultAutosome$Sex_conservative_pvalue), xlim=c(-5, 5), ylab="Sex p-value [-log10, conservative]", xlab="Sex fold change [conservative]", nbin=256, main="Autosome genes")
    abline (h=-log10(pValueThreshold), lty=2)
    # xyLim=c(0, ceiling(max (c(-log10(result$SEX_raw_pvalue), -log10(result$SEX_pvalue)))))
    # smoothScatter (-log10(result$SEX_raw_pvalue), -log10(result$SEX_pvalue), xlab="Donor Sex [not normalized] pvalue",
    #       ylab="Donor Sex pvalue", xlim=xyLim, ylim=xyLim)
    # abline (h=-log10(pValueThreshold), v=-log10(pValueThreshold))

    #permutation plots, no useful at this point.
    # numPerms=(1/min (c(result$SEX_permuted_p, result$FRACTION_X_permuted_p)))-1
    # maxXY=c(0, ceiling (log10(numPerms)))
    # smoothScatter (-log10(result$SEX_pvalue), -log10(result$SEX_permuted_p), xlim=maxXY, ylim=maxXY, xlab="Sex empiric p-value", ylab="Sex permuted p-value",
    #                main="Comparison of empiric and permuted p-values")
    # abline (0,1, col="red")
    # smoothScatter (-log10(result$FRACTION_X_pvalue), -log10(result$FRACTION_X_permuted_p), xlim=maxXY, ylim=maxXY, xlab="Sex empiric p-value", ylab="Sex permuted p-value",
    #                main="Comparison of empiric and permuted p-values")
    # abline (0,1, col="red")

    if (!is.null(outPDF)) dev.off()
    return (resultGeneList)

}

#for a list of genes, plot the results manhatten plot style.
#genes=c(resultGeneList$resultAutosomeSex, resultGeneList$resultXSex);pvalColumn="SEX_pvalue"
#genes=result$gene
manhattenPlotResults<-function (genes, pvalColumn="SEX_pvalue", result, contigGroupsFile, strTitle="", ylim=NULL) {
    replaceChr<-function (x, start, end) {
        idx=which(x$chr==start)
        if (length(idx)>0) x[idx,]$chr=end
        return (x)
    }

    contigGroups=DropSeq.eqtl::getContigGroups(contigGroupsFile)
    x=result[match(genes, result$gene),]

    x=replaceChr(x, contigGroups$chrX, 23)
    x=replaceChr(x, contigGroups$chrY, 24)
    x=replaceChr(x, contigGroups$chrMT, 25)
    x$chr=sub("chr", "", x$chr)
    idx=which(!is.na(match(x$chr, 1:25)))
    x=x[idx,]
    x$chr=as.numeric(x$chr)
    if (is.null(ylim)) ylim=c(0, ceiling(max(-log10(result[[pvalColumn]]))))
    #setting the SNP name as gene to get rid of the warning
    #x=x[order(x[[pvalColumn]]),]
    #top10AutosomeSNPs=head (x[x$chr!=23,], n=10)$gene
    qqman::manhattan(x, chr = "chr", bp = "start", snp="gene", p=pvalColumn, main=strTitle, ylim=ylim, col=c("darkslateblue","cadetblue"))
}

formatPvalue<-function (x, threshold=0.001) {
    if (x>=threshold) return (round(x, digits=round(log10(1/threshold))))
    format(x, digits=2, scientific=T)
}


#geneName="ELAVL4"
plotExampleGeneXXDosage<-function (geneName, expData, expDataXReactivationNormalized, covars, result, plotTitlePrefix=NULL) {
    df=data.frame(id=colnames(covars)[-1], sex=as.numeric (covars[covars$id=="SEX",-1]), x_reactivation=as.numeric (covars[covars$id=="FRACTION_X",-1]),
                  expression=as.numeric (expData[expData$id==geneName,-1]),
                  fit=as.numeric(expDataXReactivationNormalized[expDataXReactivationNormalized$id==geneName,-1]),  stringsAsFactors = F)

    x=result[result$gene==geneName,]
    #last iteration of analysis.
    old.par <- par(no.readonly = TRUE)
    par(mfrow=c(2,2), mar=c(5,5,2,1), oma=c(0,0,3,0))
    plotSEX(df$sex, df$expression, x, removeFracXEffects=F, showPvalue=F, cex.axis=1.25, cex.lab=1.25)
    df=plotFracXDonorXXNormalized(df, x, cex.axis=1.25, cex.lab=1.25)

    #normalized gene expression chromosome x

    #plot (df$fit, as.numeric(expDataXReactivationNormalized[expDataXReactivationNormalized$id==geneName, match(df$id, colnames(expDataXReactivationNormalized)),])
    plotSEX(df$sex, df$fit, x, removeFracXEffects=T, showPvalue=T, cex.axis=1.25, cex.lab=1.25)
    #legend("bottom",ncol=2, legend=c("XX", "XY"), fill=c("red", "blue"))
    # plot (df$sex, df$fit, xlab="genomic copies X chromosome", ylab="gene expression [fitted to x_reactivation]", col=cols[df$sex], main=strTitle)
    # abline (fit, lty=2)
    r=list(XY=df[df$sex==1,]$expression, XX=df[df$sex==2,]$expression, XX_FIT=df[df$sex==2,]$fit)
    boxplot (r, names=c("XY", "XX raw", "XX fitted"), main=paste(geneName, "X reactivation fitting results"), cex.axis=1.25, cex.lab=1.25, ylab="gene expression")
    titleStr=paste(geneName, paste(x$chr, ":", x$start, "-", x$end, sep=""))
    if (!is.null(plotTitlePrefix)) titleStr=paste(plotTitlePrefix, titleStr)
    mtext(titleStr, line=1, side=3, outer=T, cex=1.5)
    par(old.par)
}

# 1) For the subset of XX individuals, regression : expression ~ X reactivation
# 1a) Run permutation testing by permuting expression labels
testXReactivation<-function (expData, covars, covarName="FRACTION_X", numPermutations=1000) {
    s=as.numeric(covars[covars$id=="SEX",-1])

    idxDonors=which(s==2) #without the ID column, this is off by one.
    expDataXX=expData[,c(1,idxDonors+1)]
    phenotype=covars[,c(1,idxDonors+1)]
    r=DropSeq.eqtl::runOneCovarMatrixeQTL (covarName=covarName, phenotype, expDataXX, numPermutations=numPermutations)
}

testSexNonParametric<-function (expData, covars, covarName="SEX") {
    sex=as.numeric(covars[covars$id==covarName,-1])
    #what is the index of the XX donors?  Need this for grouping.
    idxXX=which(sex==2)
    expDataT=t(expData[,-1]); colnames(expDataT)=expData[,1]
    result=apply(expDataT, 2, function (x, idxXX) wilcox.test(x=jitter(x[idxXX]),y=jitter(x[-idxXX]))$p.value, idxXX)

}

#this removes the effects of X reactivation on the XX donors.  XY donors remain the same.
#this is a special case where we normalize the XX donors to have expression equivilent to an X reactivation of 1.
#this would probably be best to also apply to XY individuals, though it'll have little effect, it will be more consistent to apply it.
removeXReactivationEffectOnExpressionXXDonors<-function (expData, covars) {
    s=as.numeric(covars[covars$id=="SEX",-1])

    idxDonors=which(s==2) #without the ID column, this is off by one.
    expDataXX=expData[,c(1,idxDonors+1)]
    covarsXX=covars[,c(1,idxDonors+1)]

    #BESPOKE VERSION.  We're going to train on XX donors
    expDataT=t(expDataXX[,-1]); colnames(expDataT)=expDataXX[,1]
    xReactivation=unlist (covarsXX[which(covarsXX$id=="FRACTION_X"),-1])
    z=lm(expDataT ~ xReactivation)
    coef=coef(z)

    #apply to all the data.
    for (i in 1:dim (expDataT)[2]) {
        expDataT[,i]=expDataT[,i]-(coef[2,i]*(xReactivation-1))
    }
    expDataXXNorm=t(expDataT)

    idxDonor=match(colnames(expDataXXNorm), colnames(expData))
    expData[,idxDonor]=expDataXXNorm
    return (expData)
}

#this learns the effects of X reactivation on the XX donors.
#this is then applied to all donors.  The effect is relatively large on XX donors, and very modest on XY donors.
#Maybe more consistent than removeXReactivationEffectOnExpressionXXDonors

removeXReactivationEffectOnExpressionAllDonors<-function (expData, covars) {
    s=as.numeric(covars[covars$id=="SEX",-1])

    idxDonors=which(s==2) #without the ID column, this is off by one.
    expDataXX=expData[,c(1,idxDonors+1)]
    covarsXX=covars[,c(1,idxDonors+1)]

    #BESPOKE VERSION.  We're going to train on XX donors
    expDataT=t(expDataXX[,-1]); colnames(expDataT)=expDataXX[,1]
    xReactivation=unlist (covarsXX[which(covarsXX$id=="FRACTION_X"),-1])
    z=lm(expDataT ~ xReactivation)
    coef=coef(z)

    #apply to all the data.
    #abstracted method
    expDataFixed=applyXReactivationToExpression (expData, covars, slopes=coef[2,])
    return (expDataFixed)
}

calculateConfIntervalsXReactivationEffectOnExpressionAllDonors<-function (expData, covars, confInterval=0.99) {
    s=as.numeric(covars[covars$id=="SEX",-1])

    idxDonors=which(s==2) #without the ID column, this is off by one.
    expDataXX=expData[,c(1,idxDonors+1)]
    covarsXX=covars[,c(1,idxDonors+1)]

    #BESPOKE VERSION.  We're going to train on XX donors
    expDataT=t(expDataXX[,-1]); colnames(expDataT)=expDataXX[,1]
    xReactivation=unlist (covarsXX[which(covarsXX$id=="FRACTION_X"),-1])
    z=lm(expDataT ~ xReactivation)
    coef=coef(z)

    #V2
    z3=apply(expDataXX[,-1], 1, function(x) lm(x ~ xReactivation))
    coef2=do.call(rbind, lapply(z3, coef))

    #which(t(coef)!=coef2) # this is length 0, we're good.

    #what about lower/upper confidence intervals?
    #conf=confint (z, "xReactivation", level=0.95)
    conf=do.call(rbind, lapply(z3, function (x) confint(x, "xReactivation", level=confInterval)))

    fitDF=data.frame(gene=expData$id, intercept=coef2[,1], slope=coef2[,2], min_slope=apply(conf, 1, min),
                  max_slope=apply(conf, 1, max), stringsAsFactors = F)

    return (fitDF)

}

applyXReactivationToExpression<-function (expData, covars, slopes) {
    #rebuild expression and X reactivation for all donors.
    expDataT=t(expData[,-1]); colnames(expDataT)=expData[,1]
    xReactivation=unlist (covars[which(covars$id=="FRACTION_X"),-1])

    for (i in 1:dim (expDataT)[2]) {
        expDataT[,i]=expDataT[,i]-(slopes[i]*(xReactivation-1))
    }
    expDataNorm=t(expDataT)

    idxDonor=match(colnames(expDataNorm), colnames(expData))
    expData[,idxDonor]=expDataNorm
    return (expData)
}



convertToXReactivationPhenotype<-function (covars) {
    #fracX/(mean (fracX[s==1]))
    fracX=as.numeric(covars[covars$id=="FRACTION_X",-1])
    s=as.numeric (covars[covars$id=="SEX",-1])
    reactivation=fracX/(mean (fracX[s==1]))
    covars[covars$id=="FRACTION_X",-1]=reactivation
    return (covars)

}






#df has columns for sex, x_reactivation, expression
plotFracXDonorXXNormalized<-function (df, x, ...) {

    #need this for the slope of the fit.
    fitXReactivation=lm (expression ~ x_reactivation, data=df[df$sex==2,])
    pval=x$FRACTION_X_pvalue
    beta=x$FRACTION_X_beta
    cols=c("blue", "red")
    strTitle=paste("Fit to XX donors \n pval [", formatPvalue(pval), "], beta [", round(beta,3), "]")
    # if (!is.null(x)) {
    #     strTitle=paste(strTitle, "\n", "permuted P [", formatPvalue(x$FRACTION_X_permuted_p), "] FDR [", formatPvalue(x$FRACTION_X_FDR), "]")
    # }
    plot (df$x_reactivation, df$expression, xlab="normalized gene expression chromosome X", ylab="gene expression", main=strTitle, xlim=c(1,2), col=cols[df$sex], ...)
    abline (fitXReactivation)
    return (df)
}

plotSEX<-function (s, m, x, removeFracXEffects=F, showPvalue=F, ...) {
    ylab="gene expression"
    #if true, regress the expression against FRACTION_X, and replace the expression with the residuals.
    if (removeFracXEffects)
        ylab="gene expression [X reactivation removed]"

    z=lm(m ~ s)

    plot (s, as.numeric(m), type='n', axes=F, xlab="", ylab=ylab, xlim=c(0.75,2.25), ...)
    #axis (1, at=1:2, labels=1:2, ...)
    axis (1, at=1:2, labels=c("XY", "XX"), ...)
    axis(2, ...)
    idx=which(s==1)
    points(s[idx], m[idx], col="blue")
    points(s[-idx], m[-idx], col="red")
    yStart=z$coefficients[[1]]+(min(s)*z$coefficients[[2]])
    yEnd=z$coefficients[[1]]+(max(s)*z$coefficients[[2]])
    lines (x=c(1,2), y=c(yStart, yEnd), lty=2, col="black")
    if (showPvalue) {
        #title(paste("p-value [", formatPvalue(x$SEX_pvalue), "]\npermuted P [", formatPvalue(x$SEX_permuted_p), "] FDR [", formatPvalue(x$SEX_FDR), "]"))
        title(paste("p-value [", formatPvalue(x$SEX_pvalue), "] beta [", round (x$SEX_beta,3), "]"))
    } else {
        title(paste("p-value [", formatPvalue(summary(z)$coefficients[2,4]), "]"))
    }
    #legend("topright",ncol=2, legend=c("XX", "XY"), fill=c("red", "blue"))
}

#wanted to look at DZIP1, KIAA0753
testOne<-function (geneName="DZIP1", covars, expDataXReactivationNormalized) {

    e=as.numeric(expDataXReactivationNormalized[expDataXReactivationNormalized$id==geneName,-1])
    p=as.numeric(covars[covars$id=="SEX",-1])
    fracX=as.numeric(covars[covars$id=="FRACTION_X",-1])
    # idx=which.max(e)
    # e=e[-idx]
    # p=p[-idx]
    fit=lm(e ~ p)
    pEmpiric=coef(summary(fit))[2,4]
    pTest=t.test(x=e[p==1], y=e[p==2])$p.value

    plot (p, e, xlab="SEX", ylab="expression", main=paste(geneName, "linear pval", formatPvalue(pEmpiric), " t-test pval", formatPvalue(pTest)), col=p)
    abline (fit)

    summary(lm(e ~ p + fracX))
    summary(lm(e ~ p * fracX))

    runOne<-function (e, p) {
        pp=sample(p, size=length(p))
        coef(summary(lm(e ~ pp)))[2,4]
    }
    runTwo<-function (e, p) {
        ee=sample (e, size=length(e))
        coef(summary(lm(ee ~ p)))[2,4]
    }
    runThree<-function (e, p) {
        ee=sample (e, size=length(e))
        t.test(x=e[p==1], y=ee[p==2])$p.value
    }

    numPerm=10000
    pPermuted=unlist(mclapply(1:numPerm, function (x, e, p) runOne(e,p), e,p))
    pPermuted2=unlist(mclapply(1:numPerm, function (x, e, p) runTwo(e,p), e,p))
    pPermuted3=unlist(mclapply(1:numPerm, function (x, e, p) runThree(e,p), e,p))

    hist (pPermuted, breaks=100, main=paste(geneName, "distribution permuted values [permuted phenotype]\n", "num exeed p [", length(which(pPermuted<pEmpiric)), "/", length(pPermuted), "]"), xlim=range(c(pEmpiric, pPermuted)))
    abline (v=pEmpiric, col="red")
    hist (pPermuted2, breaks=100, main=paste(geneName, "distribution permuted values [permuted expression]\n", "num exeed p [", length(which(pPermuted2<pEmpiric)), "/", length(pPermuted2), "]"), xlim=range(c(pEmpiric, pPermuted2)))
    abline (v=pEmpiric, col="red")
    hist (pPermuted3, breaks=100, main=paste(geneName, "distribution t-test permuted values [permuted expression]\n", "num exeed p [", length(which(pPermuted3<pTest)), "/", length(pPermuted3), "]"), xlim=range(c(pTest, pPermuted3)))
    abline (v=pTest, col="red")
    hist (log10(pPermuted3), breaks=100, xlim=c(log10(pTest), 0), main=paste(geneName, "distribution t-test permuted values [permuted expression]\n", "num exeed p [", length(which(pPermuted3<pTest)), "/", length(pPermuted3), "]"))
    abline (v=log10(pTest), col="red")

    #what about the "best case scenario"?
    ee=sort(e)
    pp=c(rep(1, length(which(p==1))), rep(2, length(which(p==2))))
    pTest=t.test(x=ee[pp==1], y=ee[pp==2])$p.value
    fit=lm(ee ~ pp)
    pEmpiric2=coef(summary(fit))[2,4]
    plot (pp, ee, xlab="SEX", ylab="expression", main=paste(geneName, "linear pval", formatPvalue(pEmpiric2), " t-test pval", formatPvalue(pTest)), col=pp)
    abline (fit)

}

plotXIST<-function () {
    # metaCellFile="/downloads/eQTL/genderBias/d14-42_NGN2.meta_cell.expression.txt"
    # covarsFile="/downloads/eQTL/genderBias/d14-42_NGN2.maf_0.20_cisDist_10kb.covariates.txt"

    metaCellFile="/downloads/eQTL/genderBias/d5QueenBasic.meta_cell.expression.txt"
    covarsFile="/downloads/eQTL/genderBias/d5QueenBasic.maf_0.20_cisDist_10kb.covariates.txt"

    # metaCellFile="/downloads/eQTL/genderBias/d21STORMI.meta_cell.expression.txt"
    # covarsFile="/downloads/eQTL/genderBias/d21STORMI.maf_0.20_cisDist_10kb.covariates.txt"

    # metaCellFile="/downloads/eQTL/genderBias/DropulationVillageC_all.meta_cell.expression.txt"
    # covarsFile="/downloads/eQTL/genderBias/DropulationVillageC_all.meta_cell.covars.txt"

    m=read.table(metaCellFile, header=T, stringsAsFactors = F, check.names = F)
    cov=read.table(covarsFile, header=T, stringsAsFactors = F, check.names = F)
    cov=convertToXReactivationPhenotype(cov)
    x_reactivation=as.numeric (cov[cov$id=="FRACTION_X",-1])
    s=as.numeric (cov[cov$id=="SEX",-1])

    plot(x_reactivation, log10(as.numeric (m[m$GENE=="XIST",-1])+1), col=s, main="Raw XIST expression [log10+1]", xlab="X reactivation", ylab="XIST expression [raw] [log10]")
    ylim=c(-2, 2)
    totalExpression=colSums(m[,-1])
    normalizedExpression=((m[m$GENE=="XIST",-1]+1)/totalExpression)*100000
    plot (x_reactivation, log10(normalizedExpression), col=s, main="XIST Fractional expression", xlab="X reactivation", ylab="XIST expression [normalized] [log10]", ylim=ylim)

}

plotGenesCorrelatedToXreactivation<-function () {

    getCorToXReactivation<-function (expFile, covarsFile, geneLocsFile, xChromosomeString="X") {
        geneLocs=read.table(geneLocsFile, header=T, stringsAsFactors = F, sep="\t")
        expData=read.table(expFile, header=T, stringsAsFactors = F, check.names = F)
        cov=read.table(covarsFile, header=T, stringsAsFactors = F, check.names = F)
        cov=convertToXReactivationPhenotype(cov)
        x_reactivation=as.numeric (cov[cov$id=="FRACTION_X",-1])
        s=as.numeric (cov[cov$id=="SEX",-1])
        result=getCorrelationToXReactivation(s, expData, x_reactivation, geneLocs, xChromosomeString)
        return (result)
    }

    getMedianExpressionXXDonors<-function (expFile, covarsFile, geneLocsFile, xChromosomeString="X", useXChromosomeOnly=T) {
        geneLocs=read.table(geneLocsFile, header=T, stringsAsFactors = F, sep="\t")
        expData=read.table(expFile, header=T, stringsAsFactors = F, check.names = F)
        cov=read.table(covarsFile, header=T, stringsAsFactors = F, check.names = F)
        cov=convertToXReactivationPhenotype(cov)
        s=as.numeric (cov[cov$id=="SEX",-1])
        if (useXChromosomeOnly)
            idx=match(geneLocs[geneLocs$chr==xChromosomeString,]$geneid, expData$id)
        else
            idx=match(geneLocs[geneLocs$chr!=xChromosomeString,]$geneid, expData$id)
        z=expData[idx,]
        rownames (z)=z$id
        z=z[,-1]
        z=z[,which(s==2)]
        mExp=apply (z, 1, median)
        return (mExp)
    }

    expFile="/downloads/eQTL/genderBias/d14-42_NGN2.maf_0.20_cisDist_10kb.gene_expression.txt"
    covarsFile="/downloads/eQTL/genderBias/d14-42_NGN2.maf_0.20_cisDist_10kb.covariates.txt"
    geneLocsFile="/downloads/eQTL/genderBias/d14-42_NGN2.maf_0.20_cisDist_10kb.gene_locations.txt"
    expName="d14-42NGN2"
    #unique (read.table(geneLocsFile, header=T)$chr)
    r1=getCorToXReactivation(expFile, covarsFile, geneLocsFile, "X")
    mExp1X=getMedianExpressionXXDonors(expFile, covarsFile, geneLocsFile, xChromosomeString="X", useXChromosomeOnly=T)
    mExp1A=getMedianExpressionXXDonors(expFile, covarsFile, geneLocsFile, xChromosomeString="X", useXChromosomeOnly=F)
    hist (r1, breaks=100, main=expName, xlab="correlation to X reactivation phenotype")

    expFile="/downloads/eQTL/genderBias/d5QueenBasic.maf_0.20_cisDist_10kb.gene_expression.txt"
    covarsFile="/downloads/eQTL/genderBias/d5QueenBasic.maf_0.20_cisDist_10kb.covariates.txt"
    geneLocsFile="/downloads/eQTL/genderBias/d5QueenBasic.maf_0.20_cisDist_10kb.gene_locations.txt"
    expName="Queen Basic"
    mExp2X=getMedianExpressionXXDonors(expFile, covarsFile, geneLocsFile, "chrX", useXChromosomeOnly=T)
    mExp2A=getMedianExpressionXXDonors(expFile, covarsFile, geneLocsFile, "chrX", useXChromosomeOnly=F)
    r2=getCorToXReactivation(expFile, covarsFile, geneLocsFile, "chrX")
    hist (r2, breaks=100, main=expName, xlab="correlation to X reactivation phenotype")

    expFile="/downloads/eQTL/genderBias/d21STORMI.maf_0.20_cisDist_10kb.gene_expression.txt"
    covarsFile="/downloads/eQTL/genderBias/d21STORMI.maf_0.20_cisDist_10kb.covariates.txt"
    geneLocsFile="/downloads/eQTL/genderBias/d21STORMI.maf_0.20_cisDist_10kb.gene_locations.txt"
    expName="STORMI"
    mExp3X=getMedianExpressionXXDonors(expFile, covarsFile, geneLocsFile, "chrX", useXChromosomeOnly=T)
    mExp3A=getMedianExpressionXXDonors(expFile, covarsFile, geneLocsFile, "chrX", useXChromosomeOnly=F)
    r3=getCorToXReactivation(expFile, covarsFile, geneLocsFile, "chrX")
    hist (r3, breaks=100, main=expName, xlab="correlation to X reactivation phenotype")

    expFile="/downloads/eQTL/genderBias/DropulationVillageC_all.gene_expression.txt"
    covarsFile="/downloads/eQTL/genderBias/DropulationVillageC_all.meta_cell.covars.txt"
    geneLocsFile="/downloads/eQTL/genderBias/DropulationVillageC_all.gene_locations.txt"
    expName="Village C [HESC]"
    mExp4X=getMedianExpressionXXDonors(expFile, covarsFile, geneLocsFile, "X", useXChromosomeOnly=T)
    mExp4A=getMedianExpressionXXDonors(expFile, covarsFile, geneLocsFile, "X", useXChromosomeOnly=F)
    r4=getCorToXReactivation(expFile, covarsFile, geneLocsFile, "X")
    hist (r4, breaks=100, main=expName, xlab="correlation to X reactivation phenotype")

    cols=c("blue", "green", "orange", "purple")
    plot(density(r1), xlab="Correlation to X reactivation phenotype", main="", col=cols[1], lwd=3, ylim=c(0, 2.5), cex.lab=1.5, cex.axis=1.5)
    lines(density(r2), col=cols[2], lwd=3)
    lines(density(r3), col=cols[3], lwd=3)
    lines(density(r4), col=cols[4], lwd=3)
    legend("topleft", legend=c("d14-42NGN2", "Queen Basic", "STORMI", "Village C [HESC]"), fill=cols, cex=1.5)

    #plot QB vs STORMI gene by gene for correlation of autosomal genes correlation to X reactivation
    smoothScatter (r2[intersect (names (r2), names(r3))], r3[intersect (names (r2), names(r3))], xlim=c(-1, 1), ylim=c(-1,1),
                   xlab="Queen Basic", ylab="STORMI", main="Autosome gene correlation to X reactivation XX donors",
                   cex.lab=1.5, cex.axis=1.5, cex.main=1.5)

    xyLim=range (log10(c(mExp2X, mExp3X)))
    smoothScatter (log10(mExp2X[intersect (names(mExp2X), names(mExp3X))]), log10(mExp3X[intersect (names(mExp2X), names(mExp3X))]),
                   xlab="Queen Basic [log10]", ylab="STORMI [log10]", main="Average X chromosome gene expression XX donors", cex.lab=1.5, cex.axis=1.5,
                   xlim=xyLim, ylim=xyLim)
    abline (0, 1, col='red', lty=2)
    xyLim=range (log10(c(mExp2A, mExp3A)))
    smoothScatter (log10(mExp2A[intersect (names(mExp2A), names(mExp3A))]), log10(mExp3A[intersect (names(mExp2A), names(mExp3A))]),
                   xlab="Queen Basic [log10]", ylab="STORMI [log10]", main="Average autosome gene expression XX donors", cex.lab=1.5, cex.axis=1.5,
                   xlim=xyLim, ylim=xyLim)
    abline (0, 1, col='red', lty=2)

}

# prototypePlotSexRegression<-function (geneName, expData, covars, result) {
#     df=data.frame(sex=as.numeric (covars[covars$id=="SEX",-1]), x_reactivation=as.numeric (covars[covars$id=="FRACTION_X",-1]),
#                   expression=as.numeric (expData[expData$id==geneName,-1]))
#
#     x=result[result$gene==geneName,]
#     #last iteration of analysis.
#     old.par <- par(no.readonly = TRUE)
#     par(mfrow=c(2,2), mar=c(5,5,2,1), oma=c(0,0,3,0))
#     plotSEX(df$sex, df$expression, x, removeFracXEffects=F, showPvalue=T, cex.axis=1.25, cex.lab=1.25)
#     df=plotFracXDonorXXNormalized(df, x, cex.axis=1.25, cex.lab=1.25)
#
#
#     plotSEX(df$sex, df$fit, x, removeFracXEffects=T, showPvalue=F, cex.axis=1.25, cex.lab=1.25)
#     legend("bottom",ncol=2, legend=c("XX", "XY"), fill=c("red", "blue"))
#     # plot (df$sex, df$fit, xlab="genomic copies X chromosome", ylab="gene expression [fitted to x_reactivation]", col=cols[df$sex], main=strTitle)
#     # abline (fit, lty=2)
#     r=list(XY=df[df$sex==1,]$expression, XX=df[df$sex==2,]$expression, XX_FIT=df[df$sex==2,]$fit)
#     boxplot (r, names=c("XY", "XX raw", "XX fitted"), main=paste(geneName, "X reactivation fitting results"), cex.axis=1.25, cex.lab=1.25, ylab="gene expression")
#     titleStr=paste(geneName, paste(x$chr, ":", x$start, "-", x$end, sep=""))
#
#     mtext(titleStr, line=1, side=3, outer=T, cex=1.5)
#     par(old.par)
# }

plotRegressionCoefficients<-function (out.dir) {
    chrX="X"

    base.dir="/downloads/eQTL/genderBias"

    base.name="d14-42_NGN2.maf_0.20_cisDist_10kb"
    resultFile1=paste(out.dir, "/", base.name, ".xx_biased_genes.txt", sep="")

    base.name="d21STORMI.maf_0.20_cisDist_10kb"
    resultFile2=paste(out.dir, "/", base.name, ".xx_biased_genes.txt", sep="")

    a=read.table(resultFile1, header=T, stringsAsFactors = F)
    b=read.table(resultFile2, header=T, stringsAsFactors = F)
    both=intersect (a$gene, b$gene)
    a=a[match(both, a$gene),]
    b=b[match(both, b$gene),]

    pdf("/downloads/hesc_neurons_vs_ipsc_neurons.pdf")
    plot (log10(a$median_expression), log10(b$median_expression), xlim=c(-2, 4), ylim=c(-2,4), xlab="d14-d42 NGN2 expression [log10]", ylab="d21 Stormi expression [log10]", main="Median Expression")
    abline (0,1, col='red', lty=2)
    a$Sex_conservative_fc=a$Sex_conservative_beta/a$median_expression
    b$Sex_conservative_fc=b$Sex_conservative_beta/b$median_expression

    a$FRACTION_X_fc=a$FRACTION_X_beta/a$median_expression
    b$FRACTION_X_fc=b$FRACTION_X_beta/b$median_expression

    #Sex, autosomes
    xylim=range(c(a[a$chr!=chrX,]$Sex_conservative_beta, b[a$chr!=chrX,]$Sex_conservative_beta))
    plot (a[a$chr!=chrX,]$Sex_conservative_beta, b[a$chr!=chrX,]$Sex_conservative_beta, xlim=xylim, ylim=xylim, xlab="d14-d42 NGN2 Sex beta", ylab="d21 Stormi Sex betae", main="Sex beta autosomes")

    xylim=c(-10,10)
    plot (a[a$chr!=chrX,]$Sex_conservative_beta, b[a$chr!=chrX,]$Sex_conservative_beta, xlim=xylim, ylim=xylim, xlab="d14-d42 NGN2 Sex beta", ylab="d21 Stormi Sex betae", main="Sex beta autosomes")

    #Sex X chromosome
    xylim=range(c(a[a$chr==chrX,]$Sex_conservative_beta, b[a$chr==chrX,]$Sex_conservative_beta))
    plot (a[a$chr==chrX,]$Sex_conservative_beta, b[a$chr==chrX,]$Sex_conservative_beta, xlim=xylim, ylim=xylim, xlab="d14-d42 NGN2 Sex beta", ylab="d21 Stormi Sex betae", main="Sex beta X chromosome")

    xylim=c(-5,10)
    plot (a[a$chr==chrX,]$Sex_conservative_beta, b[a$chr==chrX,]$Sex_conservative_beta, xlim=xylim, ylim=xylim, xlab="d14-d42 NGN2 Sex beta", ylab="d21 Stormi Sex beta", main="Sex beta X chromosome")

    #X reactivation
    xylim=range(c(a[a$chr!=chrX,]$FRACTION_X_beta, b[a$chr!=chrX,]$FRACTION_X_beta))
    plot (a[a$chr!=chrX,]$FRACTION_X_beta, b[a$chr!=chrX,]$FRACTION_X_beta, xlim=xylim, ylim=xylim, xlab="d14-d42 NGN2 X reactivation beta", ylab="d21 Stormi X reactivation beta", main="X reactivation beta autosomes")

    xylim=range (-100, 100)
    plot (a[a$chr!=chrX,]$FRACTION_X_beta, b[a$chr!=chrX,]$FRACTION_X_beta, xlim=xylim, ylim=xylim, xlab="d14-d42 NGN2 X reactivation beta", ylab="d21 Stormi X reactivation beta", main="X reactivation beta autosomes")

    xylim=range(c(a[a$chr==chrX,]$FRACTION_X_beta, b[a$chr==chrX,]$FRACTION_X_beta))
    plot (a[a$chr==chrX,]$FRACTION_X_beta, b[a$chr==chrX,]$FRACTION_X_beta, xlim=xylim, ylim=xylim, xlab="d14-d42 NGN2 X reactivation beta", ylab="d21 Stormi X reactivation beta", main="X reactivation beta X chromosome")

    xylim=c(-20, 50)
    plot (a[a$chr==chrX,]$FRACTION_X_beta, b[a$chr==chrX,]$FRACTION_X_beta, xlim=xylim, ylim=xylim, xlab="d14-d42 NGN2 X reactivation beta", ylab="d21 Stormi X reactivation beta", main="X reactivation beta X chromosome")

    dev.off()

}

calculateGenomicInflation<-function (pvals) {
    asChiSq <- qchisq(1-pvals,1)
    result=median(asChiSq)/qchisq(0.5,1)
    return (result)
}

qqPlot<-function (pvals, main="") {
    lambda=calculateGenomicInflation(pvals)
    qqman::qq(pvals, main=paste(main, "\nLambda [", round(lambda,2),"]", sep=""))
}
#
