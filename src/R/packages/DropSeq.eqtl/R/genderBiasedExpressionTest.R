# are there any genes that are explained by sex alone?
# new: are there geness that are explauned by dosage of the X chromosome by genotype or phenotype?

#
#
#library (data.table);library (MatrixEQTL);library (qvalue)
#findGenderBiasedGenesFromBaseFiles(base.dir="/downloads/eQTL_Test/test_1_covar", base.name="small.maf_0.20_cisDist_10kb", numPermutations=1000)
#findGenderBiasedGenesFromBaseFiles(base.dir="/downloads/eQTL_Test/test_2_covar", base.name="small.maf_0.20_cisDist_10kb", numPermutations=1000)

#' A convienience method for findGenderBiasedGenes
#'
#' @param base.dir The directory where the input files are located.
#' @param base.name The prefix name of the input files (usually an experiment name)
#' @param out.dir Where should the output files be generated?
#' @inheritParams findGenderBiasedGenes
#' @export
findGenderBiasedGenesFromBaseFiles<-function (base.dir, base.name, out.dir=base.dir, covar="SEX", numPermutations=100) {
    # Gene expression file name
    expression_file_name = getExpressionFileName(base.dir, base.name)
    gene_location_file_name = getGeneLocationFileName(base.dir, base.name)

    # Covariates file name
    covariates_file_name = getCovariateFileName(base.dir, base.name)

    #outfiles
    outPDFExemplars=paste(out.dir, "/", base.name, ".xDosage_biased_genes.pdf", sep="")
    outPDFSummary=paste(out.dir, "/", base.name, ".xDosage_biased_genes.summary.pdf", sep="")
    outFile=paste(out.dir, "/", base.name, ".xDosage_biased_genes.txt", sep="")

    findGenderBiasedGenes(expression_file_name, covariates_file_name, covar=covar, gene_location_file_name, numPermutations=numPermutations,
                          outPDFExemplars, outPDFSummary, outFile)


}


#' Find genes that are biased by the gender of the donors.
#'
#' Runs linear regression of gender (encoded as 1 for male, 2 for female - ie: copies of the X chromosome) vs expression of the gene.
#' Also runs permutations of the donor gender to generate a null, and permuted results by q-value.
#'
#' @param expression_file_name A matrix of expression data, 1 row per gene, 1 column per sample.
#' @param gene_location_file_name The location of each gene.  Same number of lines as the expression_file_name, encodes the gene name, chromosome, start, end.  Optional.
#' @param covariates_file_name A matrix of 1 or more covariates.  Each row is a covariate, each column a sample.  Contains a row with the donors' gender.
#' @param covar The row name that contains the donor copies of the X chromosome in the covariates file.
#' @param numPermutations How many permutations should be run to generate permuted p-values
#' @param outPDFExemplars Output each gene as a plot
#' @param outPDFSummary Output the distribution of p-vals as a set of plots.
#' @param outFile The report of each gene tested and the p-values and effect sizes.
#' @import data.table qvalue utils
#' @export
findGenderBiasedGenes<-function (expression_file_name, covariates_file_name, covar="SEX", gene_location_file_name=NULL, numPermutations=10000, outPDFExemplars=NULL, outPDFSummary=NULL, outFile=NULL) {
    validateSingleFileExists(expression_file_name)
    validateSingleFileExists(covariates_file_name)

    expData=read.table(expression_file_name, header=T, stringsAsFactors = F, sep="\t")
    b=read.table(covariates_file_name, header=T, stringsAsFactors = F, sep="\t")
    geneLoc=read.table(gene_location_file_name, header=T, stringsAsFactors = F, sep="\t")
    #the primary covariate to test.
    covars=b[match(covar, b$id),]
    #covariates other than the primary one.
    otherCovars=b[match(setdiff(b$id, covars$id),b$id),]
    cvrt=SlicedData$new()

    if (dim (otherCovars)[1]>0)
        cvrt=convertToSlicedData(b[match(setdiff(b$id, covars$id),b$id),])

    if (length(which(colnames (expData)!=colnames(covars)))>0) warning ("expression and covariance matrixes out of frame")

    #the fast way using eQTL matrix.
    ee=convertToSlicedData(expData)
    s=convertToSlicedData(covars)

    #if there are multiple covars they are run independently by runOneFasst
    result=runOneFast(ee,s, cvrt)
    result=runPermutations (result, ee, s, cvrt, numPermutations=numPermutations)

    #result$`p-value`<-as.numeric (result$`p-value`)

    #add the median expression
    result=addMedianExpression(result, ee)

    #sort by pvalue.
    result=result[order(result$fdr, result$beta),]

    idx=match(result$gene, geneLoc$geneid)
    result$chr=geneLoc[idx,]$chr
    result$start=geneLoc[idx,]$s1
    result$end=geneLoc[idx,]$s2

    #make it look like the permuted results output.
    idx=c(2,1,3:dim(result)[2])
    result=result[,idx]


    resultAutosome=result[result$chr!="X" & result$chr!="Y",]
    resultAutosomeFDRFilteredGenes=resultAutosome[resultAutosome$qvalue<=0.05,]$gene
    paste ("autosome eQTLs found", length(resultAutosomeFDRFilteredGenes))

    if (!is.null(outPDFExemplars) & (length(resultAutosomeFDRFilteredGenes)>0)) pdf(outPDFExemplars)
    if (length(resultAutosomeFDRFilteredGenes)>0) sapply(resultAutosomeFDRFilteredGenes, plotExemplarGender, expData, result, covars)
    if (!is.null(outPDFExemplars) & (length(resultAutosomeFDRFilteredGenes)>0)) dev.off()

    #positive controls
    #xGenes=result[result$chr=="X" & result$fdr<=0.05,]$id
    #if (!is.null(outPDFControls)) pdf(outPDFControls)
    #sapply(xGenes, plotExemplarGender, expData, result, sex)
    #if (!is.null(outPDFControls)) dev.off()

    #plot distribution
    if (!is.null(outPDFSummary)) {
        pdf(outPDFSummary)
        strTitle=paste("Gene expression affected by XX copy number", sep="")
        #if (permuteDonors) strTitle=paste(strTitle, "[permuted]", sep=" ")
        z=hist(result$`p-value`, breaks=100, main=strTitle, xlab="empiric pvalue of gene", ylab="number of genes", col="grey", cex.main=1.5)
        #plotExemplarOther(result[1,]$id, expData, result, sex)
        z2=hist(result[result$chr=="X",]$`p-value`, breaks=z$breaks, col="blue", add=T)
        legend("topright", legend=c("all data", "x chromosome"), fill=c("grey", "blue"))
        z3=hist(result[result$chr!="X",]$`p-value`, breaks=z$breaks, col="grey", main="Autosome gene expression affected by donor XX copy number", xlab="empiric pvalue of gene", ylab="number of genes")
        dev.off()
    }

    #drop the fdr column, don't need that AND q-value...
    #format columnns for writing out.
    result=result[,-match("fdr", colnames(result))]
    #rename SNP to phenotype
    colnames(result)[which(colnames(result)=="SNP")]<-"covariate"
    result$`t-stat`=round(result$`t-stat`, 5)
    result$`p-value`=format(result$`p-value`, scientific=T, digits=5)
    result$FDR=format(result$FDR, scientific=T, digits=5)
    result$beta=round(result$beta, 5)
    result$permuted_p=format(result$permuted_p, scientific=T, digits=5)
    result$qvalue=format(result$qvalue, scientific=T, digits=5)
    result$median_expression=round(result$median_expression, 5)
    result$effect_size=round(result$effect_size, 5)
    if (!is.null(outFile)) write.table(result, outFile, row.names=F, col.names=T, quote=F, sep="\t")

}

getSkewedDonors<-function (sex, lowerBound=-2.2, upperBound=-1.4) {
    x=as.numeric (sex)
    idx=which(x<=lowerBound | x>= upperBound)
    return (names(sex)[idx])
}

#x=expData[1,-1]
runOneGene<-function (x, sex) {
    z=lm(as.numeric(x) ~ as.numeric(sex))
    zz=summary(z)
    df=data.frame(intercept=zz$coefficients[1,1], beta=zz$coefficients[2,1], pval=zz$coefficients[2,4], r.squared=zz$r.squared)
    return (df)
}


#PVALUES are the same as running the long slow linear regression via lm().
runOneFast<-function (ee, s, cvrt=SlicedData$new()) {
    invisible(capture.output(z<-Matrix_eQTL_engine(snps=s, gene=ee, cvrt=cvrt, output_file_name=NULL, pvOutputThreshold=1, verbose=F),type="message"))
    zz=z$all$eqtls
    colnames(zz)[match(c("snps", "statistic", "pvalue"), colnames(zz))]=c("SNP", "t-stat","p-value")
    return (zz)
}





#result=runPermutations (resultFast, ee, s, numPermutations=numPermutations)
runPermutations<-function (resultFast, ee, s, cvrt, numPermutations=1000, random.seed=1) {
    set.seed(random.seed)
    resultFast$exceeds_empiric_p=0
    for (i in 1:numPermutations) {
        if (i%%100==0)cat (paste("Permutation [", i, "] of [", numPermutations, "]\n"))
        #permute genders.
        s$ColumnSubsample(sample (1:s$nCols()))
        permResult=runOneFast(ee, s, cvrt)
        idx=match(resultFast$gene, permResult$gene)
        permResult=permResult[idx,]
        idxExeeds=which(permResult$`p-value`<=resultFast$`p-value`)
        resultFast[idxExeeds,]$exceeds_empiric_p=resultFast[idxExeeds,]$exceeds_empiric_p+1
    }

    resultFast$permuted_p=(resultFast$exceeds_empiric_p+1)/(numPermutations+1)
    q=qvalue(resultFast$permuted_p)
    resultFast$fdr=p.adjust(resultFast$permuted_p, method="BH")
    resultFast$qvalue=q$qvalues
    return (resultFast)
}





#geneName="HINT1"

#' Plot a single gene gender vs expression with fit.
#'
#' @param geneName The gene name
#' @param expData A SlicedData set of expression data
#' @param result The eQTL-type scan of all genes
#' @param sex A SlicedData set of the sex of the donors
#' @import grDevices
plotExemplarGender<-function (geneName, expData, result, sex) {
    s=as.numeric (sex[-1])
    m=as.numeric (expData[expData$id==geneName,-1])
    plot (s, as.numeric(m), type='n', axes=F, xlab="copies of X chromosome", ylab="gene expression", xlim=c(0.5,2.5))
    axis (1, at=1:2, labels=1:2)
    axis(2)
    idx=which(s==1)
    points(s[idx], m[idx], col="blue")
    points(s[-idx], m[-idx], col="red")

    z=lm (m ~ s)

    gene=result[result$gene==geneName,]
    #strTitle=paste(gene$gene, " ", gene$chr, ":", gene$start, "-", gene$end, "\n beta ", round(gene$beta,3), " effect size ", round(gene$effect_size,4), " pval [-log10] ", round(-log10(gene$pval),2), " q-value ", round (gene$q_value,4), sep="")
    strTitle=paste(gene$gene, " ", "\n beta ", round(gene$beta,3), " effect size ", round(gene$effect_size,4), " pval [-log10] ", round(-log10(gene$`p-value`),2), " q-value ", round (gene$qvalue,4), sep="")
    yStart=z$coefficients[[1]]+(1*z$coefficients[[2]])
    yEnd=z$coefficients[[1]]+(2*z$coefficients[[2]])
    #abline (r, lty=2)
    lines (x=c(1,2), y=c(yStart, yEnd), lty=2, col="black")
    title (strTitle)
}


# plotExemplarOther<-function (geneName, expData, result, sex) {
#     m=expData[expData$id==geneName,-1]
#     plot (as.numeric(sex), as.numeric(m), axes=T, xlab="covariate", ylab="gene expression")
#     gene=result[result$id==geneName,]
#     strTitle=paste(gene$id, "\n", gene$chr, ":", gene$start, "-", gene$end, "\n beta ", round(gene$beta,3), " pval[-log10] ", round(-log10(gene$pval),2), sep="")
#     abline(gene$intercept,gene$beta)
#     title (strTitle)
# }


# runOneManual<-function (expData, sex) {
#     dd=t(expData[,-1]); colnames(dd)=expData[,1]
#     z=lm(dd ~ as.numeric (sex[1,-1]))
#     s=summary(z)
#     #df=data.frame(gene=expData[,1], chrX_expression_dosage_phenotype_pval=as.numeric (sapply(s, function (x) x$coefficients[2,4]))
#     #              , gender_pval=as.numeric (sapply(s, function (x) x$coefficients[3,4])), stringsAsFactors = F)
#     df=data.frame(gene=expData[,1], pval=as.numeric (sapply(s, function (x) x$coefficients[2,4])) , stringsAsFactors = F)
#     rownames (df)=NULL
#     return (df)
# }



