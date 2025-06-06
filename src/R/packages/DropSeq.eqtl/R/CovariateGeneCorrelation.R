
#source ("/Users/nemesh/dropseqrna/transcriptome/R/QTL/prepareEQTLData/buildExpressionMatrix.R")
#source ("/Users/nemesh/dropseqrna/transcriptome/R/packages/DropSeq.eqtl/R/genderBiasedExpressionTest.R")
#source ("~/dropseqrna/transcriptome/R/packages/DropSeq.eqtl/R/Matrix_eQTL_CommonFunctions.R")

#library (data.table); library(qvalue); library (MatrixEQTL);library (yaml); library (parallel); library (DropSeq.eqtl);

# base.dir="/downloads/eQTL/genderBias"
# base.name="d14-42_NGN2.maf_0.20_cisDist_10kb"
# covar=c("SEX", "FRACTION_X"); numPermutations=100; random.seed=1
# contigGroupsFile="/downloads/h37_m38.contig_groups.yaml"
#
# out.dir="/downloads/eQTL/genderBias/out"
# filterCovarName="SEX"; filterCovarValue=NULL;renormAutosomeExpression=T

#input files
# expression_file_name = getExpressionFileName(base.dir, base.name)
# gene_location_file_name = getGeneLocationFileName(base.dir, base.name)
# covariates_file_name = getCovariateFileName(base.dir, base.name)
#
# #outfiles
# outPDFExemplars=paste(out.dir, "/", base.name, ".linear_regression.xDosage_biased_genes.pdf", sep="")
# outPDFSummary=paste(out.dir, "/", base.name, ".linear_regression.xDosage_biased_genes.summary.pdf", sep="")
# outFile=paste(out.dir, "/", base.name, ".linear_regression.xDosage_biased_genes.txt", sep="")
# permutedResultsFile=outFile
# numThreads=4





#' Find genes with biased expression based on one or more covariates.
#'
#' Runs linear regression on all covariate terms linearly combined to explain the expression of a gene.
#' Also runs permutations of the donor label to generate a null, and permuted results by q-value.
#'
#' @param expression_file_name A matrix of expression data, 1 row per gene, 1 column per sample.
#' @param gene_location_file_name The location of each gene.  Same number of lines as the expression_file_name, encodes the gene name, chromosome, start, end.  Optional.
#' @param covariates_file_name A matrix of 1 or more covariates.  Each row is a covariate, each column a sample.  Contains a row with the donors' gender.
#' @param numPermutations How many permutations should be run to generate permuted p-values
#' @param numThreads The number of threads to run processing on.  Defaults to 1.
#' @param outFile The report of each gene tested and the p-values and effect sizes.
#' @param filterCovarName Filter donors by this covariate name to the subset that have the filterCovarValue
#' @param filterCovarValue Filter donors by this covariate name to the subset that have the filterCovarValue
#' @param contigGroupsFile A YAML file containing each contig name and groups it belongs to.  In this case the contigs should belong to the groups autosome and chrX.  If provided, renormalizes autosome expression data to not include the x chromosomne.
#' @param renormAutosomeExpression If a contigGroupsFile and this is true, then renormalize expression of the autosome to the sum of the autosome genes (instead of all genes)
#' @param covariatesToProcess If not null, filter covariates input to this list of IDs.
#' @param random.seed A random seed to ensure reproducibility across runs of the same data.
#' @import data.table qvalue utils
#' @export
findGeneCovariateCorrelation<-function (expression_file_name, covariates_file_name, gene_location_file_name=NULL,
                                        numPermutations=10000, numThreads=1, outFile=NULL, filterCovarName="SEX",
                                        filterCovarValue=NULL, contigGroupsFile=NULL, renormAutosomeExpression=F,
                                        covariatesToProcess=NULL, random.seed=1234) {
    validateSingleFileExists(expression_file_name)
    validateSingleFileExists(covariates_file_name)

    expData=read.table(expression_file_name, header=T, stringsAsFactors = F, sep="\t")
    covars=read.table(covariates_file_name, header=T, stringsAsFactors = F, sep="\t")
    geneLoc=read.table(gene_location_file_name, header=T, stringsAsFactors = F, sep="\t")

    #filter to the covariates to process
    if (!is.null(covariatesToProcess)) {
        b=intersect(covariatesToProcess, covars$id)
        covars=covars[match(b, covars$id),]
    }

    #filter on covariate value if requested.
    if (!is.null(filterCovarName) & !is.null(filterCovarValue)) {
        #keep the first ID column plus matching donor names
        idx=c(1,which(covars[covars$id==filterCovarName,]==filterCovarValue))
        covars=covars[,idx]
        expData=expData[,idx]
    }


    if (!is.null(contigGroupsFile) & renormAutosomeExpression)
        expData=renormalizeAutosomeExpression(expData, geneLoc, contigGroupsFile)

    #remove invariant covars. index of variant covars.
    idx=which(apply (covars[,-1], 1, var)>0)
    covars=covars[idx,]

    if (length(which(colnames (expData)!=colnames(covars)))>0) warning ("expression and covariance matrixes out of frame")
    result=runAnalysisManyCovariates(expData, covars, numPermutations, random.seed, numThreads)

    #add the median expression
    result=addMedianExpression(result, expData)

    #sort by pvalue.
    result=result[order(result$FULL_MODEL_permuted_p),]

    idx=match(result$gene, geneLoc$geneid)
    loc=data.frame(chr=geneLoc[idx,]$chr, start=geneLoc[idx,]$s1, end=geneLoc[idx,]$s2)

    result=data.frame(loc, result)
    write.table(result, outFile, row.names=F, col.names = T, quote=F, sep="\t")

}



getGenesPassQvalue<-function (result, qvalThreshold=0.05) {
    idx=grep("qval", colnames(result))
    idxPass=which(apply (result[,idx], 1, function (x) any(x<=qvalThreshold)))
    return (result[idxPass,])
}

#these plots are specifically written to look at the SEX and FRACTION_X covariants.


plotSEX<-function (s, m, x, removeFracXEffects=F, showPvalue=F, ...) {
    ylab="gene expression"
    #if true, regress the expression against FRACTION_X, and replace the expression with the residuals.
    if (removeFracXEffects)
        ylab="gene expression [X Reactivation removed]"

    z=lm(m ~ s)

    plot (s, as.numeric(m), type='n', axes=F, xlab="genomic copies X chromosome", ylab=ylab, xlim=c(0.75,2.25), ...)
    axis (1, at=1:2, labels=1:2, ...)
    axis(2, ...)
    idx=which(s==1)
    points(s[idx], m[idx], col="blue")
    points(s[-idx], m[-idx], col="red")
    yStart=z$coefficients[[1]]+(min(s)*z$coefficients[[2]])
    yEnd=z$coefficients[[1]]+(max(s)*z$coefficients[[2]])
    lines (x=c(1,2), y=c(yStart, yEnd), lty=2, col="black")
    if (showPvalue) {
        title(paste("p-value [", formatPvalue(x$SEX_pvalue), "]\npermuted P [", formatPvalue(x$SEX_permuted_p), "] qvalue [", formatPvalue(x$SEX_qval), "]"))
    } else {
        title(paste("p-value [", formatPvalue(summary(z)$coefficients[2,4]), "]"))
    }
    #legend("topright",ncol=2, legend=c("XX", "XY"), fill=c("red", "blue"))
}


#plotFracX<-function (fracX, m, x, removeSexEffects=F, showPvalue=F, ...) {
#    ylab="gene expression"
#    if (removeSexEffects)
#        ylab="gene expression [X chr CN removed]"
#
#    z=lm(m ~ fracX)
#    plot (fracX, as.numeric(m), type='n', axes=T, xlab="X Reactivation", ylab=ylab, bty='n', ...)
#    idx=which(s==1)
#    points(fracX[idx], m[idx], col="blue")
#    points(fracX[-idx], m[-idx], col="red")
#
#    yStart=z$coefficients[[1]]+(min(fracX)*z$coefficients[[2]])
#    yEnd=z$coefficients[[1]]+(max(fracX)*z$coefficients[[2]])
#    lines (x=c(min(fracX),max(fracX)), y=c(yStart, yEnd), lty=2, col="black")
#    if (showPvalue) {
#        title (paste("p-value [", formatPvalue(x$FRACTION_X_pvalue), "]\npermuted P [", formatPvalue(x$FRACTION_X_permuted_p), "] q-value [", formatPvalue(x$FRACTION_X_qval), "]"))
#    } else {
#        title(paste("p-value [", formatPvalue(summary(z)$coefficients[2,4]), "]"))
#    }
#    #legend("topright",ncol=2, legend=c("XX", "XY"), fill=c("red", "blue"))
#}

formatPvalue<-function (x, threshold=0.001) {
    if (x>=threshold) return (round(x, digits=round(log10(1/threshold))))
    format(x, digits=2, scientific=T)
}


#geneName="GRIA2";eQTLResults=resultAutosomeSex;
#geneName="ST3GAL1";
#This set of plots is incredibly specific to SEX and FRACTION_X covariates.


#recalculate qvalues based on permuted pvalues, useful if subsetting data where distributions are
#different, like autosomes vs X.k
recalculateQValues<-function (result) {
    idx=grep("qval", colnames(result))
    for (i in idx) {
        cn=colnames(result)[i]
        featureName=sub("_qval", "", cn)
        pvalFeatureName=paste(featureName, "permuted_p", sep="_")
        result[,i]=qvalue(result[[pvalFeatureName]])$qvalue
    }
    return (result)

}



runAnalysisManyCovariates<-function (expData, covars, numPermutations=1000, random.seed=1, numThreads=4) {
    #pre-process data one time for speed.
    expDataT=t(expData[,-1]); colnames(expDataT)=expData[,1]
    covarsList=lapply(covars$id, function (covarName) {as.numeric(covars[covars$id==covarName,-1])})
    names(covarsList)=covars$id

    #calculate empiric results
    system.time(result<-runOneManyCovars(expDataT, covarsList))

    #run permutations
    #system.time(permResult<-runPermutationManyCovars (result, expDataT, covarsList, numPermutations, random.seed, numThreads=1))
    system.time(permResult<-runPermutationManyCovars (result, expDataT, covarsList, numPermutations, random.seed, numThreads=numThreads))

    return (permResult)
}

#run multiple covariates at the same time and capture their pvalues.
#' Run eQTL analysis one iteration to look at many covariates
#'
#' This uses data that's been pre-processed to save time.
#'
#' @param expDataT The expData matrix, transposed so each row is a donor,
#' and the gene ID's made into column names.
#' @param covarsList The covariates from the matrix transformed into a list of vectors, with
#' each list element the name of the covariate.
#'
#' @return A data frame containing each gene, and the pvalue for the each covariate as well as the full model.
runOneManyCovars<-function (expDataT, covarsList) {

    strTerms=paste(sapply(1:length(covarsList), function (x) paste("covarsList[[",x,"]]", sep="")), collapse=" + ")
    formStr=paste("expDataT ~ ", strTerms)
    z=lm(formula(formStr))
    s=summary(z)

    #extract out all the pvalues.
    idx=2:(length(covarsList)+1)
    getPvals<-function (idx, s) {
        as.numeric (sapply(s, function (x) x$coefficients[idx,4]))
    }

    getBeta<-function (idx, x) {
        as.numeric (sapply(s, function (x) x$coefficients[idx,1]))
    }


    pvals=as.data.frame(sapply(idx, getPvals, s))
    #colnames (pvals)=names(covarsList)
    colnames(pvals)= paste(names(covarsList), "_pvalue", sep="")
    betas=as.data.frame(sapply(idx, getBeta, s))
    colnames (betas)=paste(names(covarsList), "_beta", sep="")

    #full model
    getFullModelP<-function (x) {
        fstat=x$fstatistic
        as.numeric (pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE))
    }

    fullModel=sapply(s, getFullModelP)
    #this might be bugged and need expDataT[1,] ...
    df=data.frame(gene=colnames(expDataT), pvals , FULL_MODEL=as.numeric(fullModel), betas, stringsAsFactors = F)
    rownames (df)=NULL
    return (df)
}

#' Run permutations on regresssions using multiple covariates
#'
#' @param result The empiric result
#' @param expDataT The expData matrix, transposed so each row is a donor,
#' and the gene ID's made into column names.
#' @param covarsList The covariates from the matrix transformed into a list of vectors, with
#' each list element the name of the covariate.
#' @param numPermutations The number of permutations to run
#' @param random.seed The random seed
#' @param numThreads Number of threads used for processing.
#' @return A data frame of permutation results
#' @import parallel
#' @export
runPermutationManyCovars<-function (result, expDataT, covarsList, numPermutations=1000, random.seed=1, numThreads=1) {
    set.seed(random.seed)
    #alter the result dataframe to get rid of the gene column and make that a rowname
    rownames(result)=result$gene
    result=result[,-1]
    #have to ignore the beta columns until the end.
    idxBeta=grep ("beta", colnames(result))

    #track the number of times the pvalue exceeds the empiric P.
    exceedsPvalDF=data.frame(result[,-idxBeta])
    exceedsPvalDF[,1:dim(exceedsPvalDF)[2]]<-0

    #run all permutations up front and hand them out to threads.
    #listPermutations=lapply(1:numPermutations, function (x) permuteCovariates(covarsList))

    runOne<-function (covarsPermuted, expDataT, permResult, idxBeta) {
        permResult=runOneManyCovars(expDataT, covarsPermuted)
        permResult=permResult[,-1]
        idxExeeds=which(permResult[,-idxBeta] < result[,-idxBeta], arr.ind=T)
        return (idxExeeds)
    }

    numPermutationsRun=0
    logSize=100
    while (numPermutationsRun<numPermutations) {
        if (numPermutationsRun%%logSize==0)cat (paste(format(Sys.time(), "%a %b %d %X"), "Permutation [", numPermutationsRun, "] of [", numPermutations, "]\n"))
        thisBatchSize=min(numThreads, numPermutations-numPermutationsRun)
        # if the batch size is one do it the easy way.
        if (thisBatchSize==1) {
            #permute covariates
            covarsPermuted=permuteCovariates(covarsList)
            #run one.
            permResult=runOneManyCovars(expDataT, covarsPermuted)
            permResult=permResult[,-1]
            idxExeeds=which(permResult[,-idxBeta] < result[,-idxBeta], arr.ind=T)
            exceedsPvalDF[idxExeeds]=exceedsPvalDF[idxExeeds]+1
            numPermutationsRun=numPermutationsRun+1
        } else {
            listPermutations=lapply(1:thisBatchSize, function (x) permuteCovariates(covarsList))
            r=parallel::mclapply(listPermutations, runOne, expDataT, permResult, idxBeta, mc.cores =numThreads)
            for (index in 1:length(r)) {
                idxExeeds=r[[index]]
                exceedsPvalDF[idxExeeds]=exceedsPvalDF[idxExeeds]+1
            }
            numPermutationsRun=numPermutationsRun+length(r)
        }
    }

    if (numPermutationsRun==numPermutations & numPermutationsRun%%logSize==0)
        cat (paste(format(Sys.time(), "%a %b %d %X"), "Permutation [", numPermutationsRun, "] of [", numPermutations, "]\n"))

    permutedPvals=(exceedsPvalDF+1)/(numPermutations+1)
    qvals=apply(permutedPvals, 2, function (x) qvalue::qvalue(x)$qvalues)
    colnames(qvals)=paste(sub ("_pvalue", "", colnames (qvals)), "_qval", sep="")
    colnames(permutedPvals)=paste(sub ("_pvalue", "", colnames (permutedPvals)), "_permuted_p", sep="")

    finalResult=cbind(gene=rownames(result), result, permutedPvals, qvals, stringsAsFactors=F)
    rownames (finalResult)=NULL
    return (finalResult)
}

#permutes all covariates, but keeps covariates in the same order relative to each other to
#preserve the relationship of covariates within a donor.
permuteCovariates<-function (covarsList) {
    p<-function (x, idx) {
        xx=x[idx]
        return (xx)
    }
    #get a sequence from 1 to the number of donors.
    idx=1:length(covarsList[[1]])
    #generate one permutation order, and apply it to all covariates
    idx=sample(idx, length(idx), replace = F)
    covarsPermuted=lapply(covarsList, p, idx)
    return (covarsPermuted)
}


#split the result data by the contig group into autosomes and X.
splitResultsByContigGroup<-function (result, contigGroupsFile) {
    contigGroups=getContigGroups(contigGroupsFile)
    idxAutosome=sort(which(!is.na(match(result$chr, contigGroups$autosomes))))
    resultAutosome=result[idxAutosome,]
    idxX=sort(which(!is.na(match(result$chr, contigGroups$chrX))))
    resultX=result[idxX,]
    return (list(autosome=resultAutosome, X=resultX))
}

#returns a list containing a list of autosome contig names, and the X,Y,MT contig names.
#if the contigGroupsFile is not null, use that.  If it is null, then it's expected that the other variabless are
#explicitly set.  If not, then fail.

#' Get the chromosomes that are part of the autosomes/X/Y/MT groups
#'
#' @param contigGroupsFile A yaml file containing the contig groups
#' @param autosomes The autosome contigs if not using the contigGroupsFile
#' @param chrX The X contigs if not using the contigGroupsFile
#' @param chrY The Y contigs if not using the contigGroupsFile
#' @param chrMT The MT contigs if not using the contigGroupsFile
#'
#' @return A list of the contigs on the autosome/X/Y/MT
#' @import yaml
getContigGroups<-function (contigGroupsFile=NULL, autosomes=NULL, chrX=NULL, chrY=NULL, chrMT=NULL) {
    result=list()
    getContigsInGroup<-function (group, r) {
        idx=which(sapply(r, function (x, group) {group %in% x}, group))
        return (names (r)[idx])
    }

    if (!is.null(contigGroupsFile)) {
        r=yaml::read_yaml(contigGroupsFile)
        result$autosomes=getContigsInGroup("autosome", r)
        result$chrX=getContigsInGroup("X", r)
        result$chrY=getContigsInGroup("Y", r)
        result$chrMT=getContigsInGroup("MT", r)
        return (result)
    }

    if (is.null(autosomes) || is.null(chrX) || is.null (chrY) || is.null(chrMT)) {
        stop ("Must set either the contigGroupsFile OR all of autosomes/chrX/chrY/chrMT")
    }

    result$autosomes=autosomes
    result$chrX=chrX
    result$chrY=chrY
    result$chrMT=chrMT
    return (result)
}

#buildSummaryPlots(permutedResultsFile, contigGroupsFile, dataColumnType="permuted_p")
#buildSummaryPlots(permutedResultsFile, contigGroupsFile, dataColumnType="qval")

buildSummaryPlots<-function (permutedResultsFile, contigGroupsFile, dataColumnType="permuted_p") {
    result=read.table(permutedResultsFile, header=T, stringsAsFactors=F)

    z=splitResultsByContigGroup (result, contigGroupsFile)
    resultAutosome=z$autosome
    resultX=z$X

    #covarName=covarNames[2];x=resultAutosome
    plotCovarQvalueDistribution<-function (covarName, x, dataColumnType="permuted_p", strTitlePrefix="", minScoreThreshold=0.05) {
        l=paste("_",dataColumnType, sep="")
        colName=paste(covarName,l, sep="")
        nBelowThreshold=length(which(x[[colName]]<=minScoreThreshold ))
        strTitle=paste(strTitlePrefix, covarName, " [" ,nBelowThreshold, "/",  dim (x)[1], "] score <=", minScoreThreshold)
        z=hist (x[[colName]], main=strTitle, breaks=100, xlab=paste (covarName, dataColumnType), axes=T, xlim=c(0,1))
        abline (v=minScoreThreshold, col="red")
    }

    covarNames=sub("_qval", "", colnames(result)[grep("qval", colnames(result))])

    sapply(covarNames, plotCovarQvalueDistribution, resultAutosome, dataColumnType="permuted_p", strTitlePrefix="Autosomes")
    sapply(covarNames, plotCovarQvalueDistribution, resultX, dataColumnType="permuted_p", strTitlePrefix="Chr X")


}

#buildScatterPlot(permutedResultsFile, contigGroupsFile, colNames=c("SEX_permuted_p", "FRACTION_X_permuted_p"), dataLabels=c("sex permuted qval [-log10]", "fraction x permuted pval [-log10]"), useLog10=T, useNegative=T, strTitle="permuted pvalues")
#buildScatterPlot(permutedResultsFile, contigGroupsFile, colNames=c("SEX_qval", "FRACTION_X_qval"), dataLabels=c("sex permuted qval [-log10]", "fraction x permuted qval [-log10]"), useLog10=T, useNegative=T, strTitle="permuted qvalues")


#' Builds a scatter plot
#'
#' @param permutedResultsFile permuted eQTL results file
#' @param contigGroupsFile A YAML file containing each contig name and groups it belongs to.  In this case the contigs should belong to the groups autosome and chrX.  If provided, renormalizes autosome expression data to not include the x chromosomne.
#' @param colNames Column names to plot
#' @param dataLabels Plot labels for column names
#' @param useLog10 Set to true to log10 transform the data
#' @param useNegative Set true to multiply the results by -1
#' @param strTitle The title for plots.
#' @import graphics
buildScatterPlot<-function (permutedResultsFile, contigGroupsFile, colNames=c("SEX_permuted_p", "FRACTION_X_permuted_p"), dataLabels=c("sex permuted pval [-log10]", "fraction x permuted pval [-log10]"), useLog10=T, useNegative=T, strTitle="") {
    if (length(colNames)!=2 || length(dataLabels)!=2)
        stop("Must provide exactly 2 columns to plot against each other, and 2 names for those columns.")

    result=read.table(permutedResultsFile, header=T, stringsAsFactors=F)

    z=splitResultsByContigGroup (result, contigGroupsFile)
    resultAutosome=z$autosome
    resultX=z$X

    getData<-function (x, colName, useLog10=F, useNegative=F) {
        if (useLog10) {
            xx=log10(x[[colName]])
        } else {
            xx=x[[colName]]
        }
        if (useNegative) xx=xx*-1
        return (xx)
    }

    xRange=range(getData(result, colNames[1], useLog10, useNegative))
    yRange=range(getData(result, colNames[2], useLog10, useNegative))

    xAut=getData(resultAutosome, colNames[1], useLog10, useNegative)
    yAut=getData(resultAutosome, colNames[2], useLog10, useNegative)
    xX=getData(resultX, colNames[1], useLog10, useNegative)
    yX=getData(resultX, colNames[2], useLog10, useNegative)

    parOld=graphics::par()
    par(mar=c(5,5,5,8))
    plot (xRange, yRange, xlab=dataLabels[1], ylab=dataLabels[2], type='n', main=strTitle)
    points(xAut, yAut, cex=0.25, pch=16, col="blue")
    points(xX, yX, cex=0.5, pch=16, col="green")

    legend("topright", inset=c(-0.25, 0), legend=c("Autosome", "X"), fill=c("blue", "green"), xpd=TRUE)
    suppressWarnings(par(parOld))

}

#very specific for our current problem.
buildSexByFractionXScatterPlot<-function (permutedResultsFile, contigGroupsFile, colNames=c("SEX_permuted_p", "FRACTION_X_permuted_p"), dataLabels=c("sex permuted pval [-log10]", "fraction x permuted pval [-log10]"), useLog10=T, useNegative=T, strTitle="") {
    result=read.table(permutedResultsFile, header=T, stringsAsFactors=F)

    z=splitResultsByContigGroup (result, contigGroupsFile)
    resultAutosome=z$autosome
    resultX=z$X

    runScatterPlot<-function (xRange, yRange, xlab, ylab, strTitle, xAut, yAut, xX, xY) {
        old.par <- par(no.readonly = TRUE)
        par(mar=c(5,5,5,8))
        graphics::plot (xRange, yRange, xlab=xlab, ylab=ylab, type='n', main=strTitle)
        graphics::points(xAut, yAut, cex=0.25, pch=16, col="blue")
        graphics::points(xX, yX, cex=0.5, pch=16, col="green")

        graphics::legend("topright", inset=c(-0.2, 0), legend=c("Autosome", "X"), fill=c("blue", "green"), xpd=TRUE, title="Chromosome")
        suppressWarnings(graphics::par(old.par))

    }

    #pvalues
    xAut=-log10(resultAutosome$SEX_permuted_p)
    yAut=-log10(resultAutosome$FRACTION_X_permuted_p)
    xX=-log10(resultX$SEX_permuted_p)
    yX=-log10(resultX$FRACTION_X_permuted_p)
    xyRange=c(0, max(xAut,yAut,xX, yX))
    runScatterPlot(xyRange, xyRange, xlab="Sex permuted pval [-log10]", ylab="fraction x permuted pval [-log10]",
                   strTitle="Permuted pvalues", xAut, yAut, xX, yX)

    #qvalues
    xAut=-log10(resultAutosome$SEX_qval)
    yAut=-log10(resultAutosome$FRACTION_X_qval)
    xX=-log10(resultX$SEX_qval)
    yX=-log10(resultX$FRACTION_X_qval)
    xyRange=c(0, max(xAut,yAut,xX, yX))

    runScatterPlot(xyRange, xyRange, xlab="Sex permuted qval [-log10]", ylab="fraction x permuted qval [-log10]",
                   strTitle="Permuted qvalues ", xAut, yAut, xX, yX)

    #qvalues recalculated
    resultAutosomeFixed=recalculateQValues(resultAutosome)
    resultXFixed=recalculateQValues(resultX)

    xAut=-log10(resultAutosomeFixed$SEX_qval)
    yAut=-log10(resultAutosomeFixed$FRACTION_X_qval)
    xX=-log10(resultXFixed$SEX_qval)
    yX=-log10(resultXFixed$FRACTION_X_qval)

    xyRange=c(0, max(xAut,yAut,xX, yX))
    runScatterPlot(xyRange, xyRange, xlab="Sex permuted pval [-log10]", ylab="fraction x permuted pval [-log10]",
                   strTitle="Permuted qvalues (recalculated)", xAut, yAut, xX, yX)

}

#PVALUES are the same as running the long slow linear regression via lm().
runSingleCovarFast<-function (ee, s, cvrt=MatrixEQTL::SlicedData$new()) {
    #z<-MatrixEQTL::Matrix_eQTL_engine(snps=s, gene=ee, cvrt=cvrt, output_file_name=NULL, pvOutputThreshold=1, verbose=F)
    invisible(capture.output(z<-MatrixEQTL::Matrix_eQTL_engine(snps=s, gene=ee, cvrt=cvrt, output_file_name=NULL, pvOutputThreshold=1, verbose=F),type="message"))
    zz=z$all$eqtls
    colnames(zz)[match(c("snps", "statistic", "pvalue"), colnames(zz))]=c("SNP", "t-stat","p-value")
    return (zz)
}



#result=runPermutations (resultFast, ee, s, numPermutations=numPermutations)
runSingleCovarPermutations<-function (resultFast, ee, s, cvrt, numPermutations=1000, random.seed=1) {
    set.seed(random.seed)
    resultFast$exceeds_empiric_p=0
    for (i in 1:numPermutations) {
        if (i%%100==0)cat (paste("Permutation [", i, "] of [", numPermutations, "]\n"))
        #permute genders.
        s$ColumnSubsample(sample (1:s$nCols()))
        permResult=runSingleCovarFast(ee, s, cvrt)
        idx=match(resultFast$gene, permResult$gene)
        permResult=permResult[idx,]
        idxExeeds=which(permResult$`p-value`<=resultFast$`p-value`)
        resultFast[idxExeeds,]$exceeds_empiric_p=resultFast[idxExeeds,]$exceeds_empiric_p+1
    }
    resultFast$permuted_p=(resultFast$exceeds_empiric_p+1)/(numPermutations+1)
    q=qvalue::qvalue(resultFast$permuted_p)
    resultFast$fdr=p.adjust(resultFast$permuted_p, method="BH")
    resultFast$qvalue=q$qvalues
    return (resultFast)
}


plotSexBiasedGeneExpression<-function (expression_file_name, covariates_file_name, permutedResultsFile, contigGroupsFile, qvalueThreshold=0.05, outPDFExemplars=NULL) {
    validateSingleFileExists(expression_file_name)
    validateSingleFileExists(covariates_file_name)

    expData=read.table(expression_file_name, header=T, stringsAsFactors = F, sep="\t")
    covars=read.table(covariates_file_name, header=T, stringsAsFactors = F, sep="\t")
    result=read.table(permutedResultsFile, header=T, stringsAsFactors=F)

    z=splitResultsByContigGroup (result, contigGroupsFile)
    resultAutosome=z$autosome
    resultX=z$X

    covarNames=sub("_qval", "", colnames(result)[grep("qval", colnames(result))])
    expectedCovarNames=c("SEX", "FRACTION_X")
    missingCovars=setdiff(expectedCovarNames, covarNames)
    if (length(missingCovars)>0)
        stop(paste("Must have the expected covariants to use this plot", paste(expectedCovarNames, collapse=" ")))

    #XX/XY biased but not x reactivation
    resultAutosomeSex=resultAutosome[which(resultAutosome$SEX_qval<=qvalueThreshold & resultAutosome$FRACTION_X_qval>qvalueThreshold),]
    resultXSex=resultX[which(resultX$SEX_qval<=qvalueThreshold & resultX$FRACTION_X_qval>qvalueThreshold),]

    #x reactivation but not XX/XY biased
    resultAutosomeFrac=resultAutosome[which(resultAutosome$SEX_qval>qvalueThreshold & resultAutosome$FRACTION_X_qval<qvalueThreshold),]
    resultXFrac=resultX[which(resultX$SEX_qval>qvalueThreshold & resultX$FRACTION_X_qval<=qvalueThreshold),]

    #BOTH
    resultAutosomeSex=resultAutosome[which(resultAutosome$SEX_qval<=qvalueThreshold & resultAutosome$FRACTION_X_qval<=qvalueThreshold),]
    resultXSex=resultX[which(resultX$SEX_qval<=qvalueThreshold & resultX$FRACTION_X_qval<=qvalueThreshold),]

    summaryCounts=data.frame(X_Reactivation_Autosome=dim (resultAutosomeFrac)[1], X_Reactivation_X=dim (resultXFrac)[1],
                             X_Dosage_Autosome=dim(resultAutosomeSex)[1], X_Dosage_X=dim(resultXSex)[1])

    #lots of X reactivation genes, let's do the top 20.
    resultAutosomeFrac=resultAutosomeFrac[order(resultAutosomeFrac$FRACTION_X_qval, resultAutosomeFrac$median_expression),]
    resultX=resultX[order(resultX$FRACTION_X_qval, resultX$median_expression),]

    if (!is.null(outPDFExemplars)) pdf(outPDFExemplars)

    z=barplot (log10(as.numeric(summaryCounts+1)), names.arg=c("X Reactivation", "X Reactivation", "X Dosage", "X Dosage"), ylab="counts [log10]", col=c("light blue", "light green"), main=paste("Summary of genes with qvalue <=", qvalueThreshold))
    avg=mean(range(log10(summaryCounts+1)))
    text(z[,1], rep(avg, dim(z)[1]), labels=as.numeric(summaryCounts), cex=1.5)
    legend("topright", legend=c("Autosome", "X"), fill=c("light blue", "light green"), title="Chromosome")
}

convertToXReactivationPhenotype<-function (covars) {
    #fracX/(mean (fracX[s==1]))
    fracX=as.numeric(covars[covars$id=="FRACTION_X",-1])
    s=as.numeric (covars[covars$id=="SEX",-1])
    reactivation=fracX/(mean (fracX[s==1]))
    covars[covars$id=="FRACTION_X",-1]=reactivation
    return (covars)

}



#' Find genes with biased expression based SEX or FRACTION_X
#'
#' Generates expression with SEX or FRACTION_X removed.
#' Tests covariate against adjusted expression, runs permutations.
#'
#' Converts FRACTION_X into X reactivation.  This normalizes fraction_X by the mean FRACTION_X of individuals with SEX=1.
#' This makes all SEX=1 individuals close to 1.
#'
#' @param expression_file_name A matrix of expression data, 1 row per gene, 1 column per sample.
#' @param gene_location_file_name The location of each gene.  Same number of lines as the expression_file_name, encodes the gene name, chromosome, start, end.  Optional.
#' @param covariates_file_name A matrix of 1 or more covariates.  Each row is a covariate, each column a sample.  Contains a row with the donors' gender and fraction X expression
#' @param contigGroupsFile A YAML file containing each contig name and groups it belongs to.  In this case the contigs should belong to the groups autosome and chrX.
#' If provided, renormalizes autosome expression data to not include the x chromosomne.
#' @param outFile The report of each gene tested and the p-values and effect sizes.
#' @param numPermutations How many permutations should be run to generate permuted p-values
#' @param numThreads The number of threads to run processing on.  Defaults to 1.
#' @param renormAutosomeExpression If a contigGroupsFile and this is true, then renormalize expression of the autosome to the sum of the autosome genes (instead of all genes)
#' @param random.seed A random seed to ensure reproducibility across runs of the same data.
#' @import data.table qvalue utils MatrixEQTL
#' @export
# sexBiasedGeneExpressionStepWise<-function (expression_file_name, covariates_file_name, gene_location_file_name=NULL,
#                                     contigGroupsFile=NULL, outFile=NULL, numPermutations=10000, numThreads=1,
#                                     renormAutosomeExpression=F, random.seed=1234) {
#
#     validateSingleFileExists(expression_file_name)
#     validateSingleFileExists(covariates_file_name)
#
#     expData=read.table(expression_file_name, header=T, stringsAsFactors = F, sep="\t")
#     covars=read.table(covariates_file_name, header=T, stringsAsFactors = F, sep="\t")
#
#     geneLoc=read.table(gene_location_file_name, header=T, stringsAsFactors = F, sep="\t")
#
#     if (length(which(colnames (expData)!=colnames(covars)))>0) warning ("expression and covariance matrixes out of frame")
#     expectedCovarNames=c("SEX", "FRACTION_X")
#     missingCovars=setdiff(expectedCovarNames, covars$id)
#     if (length(missingCovars)>0)
#         stop("Missing required covariates [", paste(missingCovars, collapse=","), "]")
#     covars=convertToXReactivationPhenotype(covars)
#
#     if (!is.null(contigGroupsFile) & renormAutosomeExpression)
#         expData=renormalizeAutosomeExpression(expData, geneLoc, contigGroupsFile)
#
#     #subtract effects of covariate from expression for the two covars
#     expDataSexAdjusted=removeCovarEffectOnExpression("SEX", expData, covars)
#     expDataFractXAdjusted=removeCovarEffectOnExpression("FRACTION_X", expData, covars)
#
#     #run matrix eQTL to get linear regression of SEX ~ FRAC_X_adjusted_expression
#     resultsSex=runOneCovarMatrixeQTL("SEX", covars, expDataFractXAdjusted, numPermutations=numPermutations)
#     #run matrix eQTL to get linear regression of FRAC_X ~ SEX_adjusted_expression
#     resultsFracX=runOneCovarMatrixeQTL("FRACTION_X", covars, expDataSexAdjusted, numPermutations=numPermutations)
#
#     #knit the results together.
#     resultsFracX=resultsFracX[match(resultsSex$gene, resultsFracX$gene),]
#
#     df=data.frame(gene=resultsSex$gene, SEX_pvalue=resultsSex$`p-value`, SEX_beta=resultsSex$beta, SEX_permuted_p=resultsSex$permuted_p,
#                   SEX_qval=resultsSex$qvalue, FRACTION_X_pvalue=resultsFracX$`p-value`, FRACTION_X_beta=resultsFracX$beta,
#                   FRACTION_X_permuted_p=resultsFracX$permuted_p, FRACTION_X_qval=resultsFracX$qvalue, stringsAsFactors = F)
#
#     #add the median expression
#     df=addMedianExpression(df, expData)
#
#     idx=match(df$gene, geneLoc$geneid)
#     loc=data.frame(chr=geneLoc[idx,]$chr, start=geneLoc[idx,]$s1, end=geneLoc[idx,]$s2, stringsAsFactors = F)
#     df=data.frame(loc, df, stringsAsFactors = F)
#
#     #sort by SEX permuted p
#     df=df[order(df$SEX_permuted_p),]
#
#     if (!is.null(outFile)) write.table(df, outFile, row.names=F, col.names=T, quote=F, sep="\t")
#
# }

# plotExampleGeneXXDosage<-function (geneName, expData, covars, result) {
#
#     cex.lab=1.25;cex.axis=1.25
#     s=as.numeric (covars[covars$id=="SEX",-1])
#     m=as.numeric (expData[expData$id==geneName,-1])
#     fracX=as.numeric(covars[covars$id=="FRACTION_X",-1])
#     mNoFracX=residuals(lm(m ~ fracX))+median(m)
#     mNoSex=residuals(lm(m ~ s))+median(m)
#     ylim=range (c(m, mNoFracX, mNoSex))
#
#     old.par <- par(no.readonly = TRUE)
#     par(mfrow=c(2,2), mar=c(5,5,2,1), oma=c(0,0,3,0), xpd=TRUE)
#     x=result[result$gene==geneName,]
#     titleStr=paste(geneName, paste(x$chr, ":", x$start, "-", x$end, sep=""))
#
#     plotSEX(s, m, x, showPvalue=T, ylim=ylim, cex.lab=cex.lab, cex.axis=cex.axis)
#     plotFracX(fracX, m, x, showPvalue=T, ylim=ylim, cex.lab=cex.lab, cex.axis=cex.axis)
#     plotSEX(s, mNoFracX, x, removeFracXEffects=T, showPvalue=F, ylim=ylim, cex.lab=cex.lab, cex.axis=cex.axis)
#     plotFracX(fracX, mNoSex, x, removeSexEffects=T, showPvalue=F, ylim=ylim, cex.lab=cex.lab, cex.axis=cex.axis)
#     legend("topright",ncol=2, legend=c("XX", "XY"), fill=c("red", "blue"), inset=c(0,0), xpd=T)
#     mtext(titleStr, line=1, side=3, outer=T, cex=1.5)
#     invisible(par(old.par))
#
#
# }
#
