#split the result data by the contig group into autosomes and X.

#' Split a data set into a list of data frames by defined contig groups
#'
#' @param df A data frame that has a column defining the contig name
#' @param contigColumn The name of the column that has contig names
#' @param contigGroupsFile A YAML file containing the groups contigs belong to
#'
#' @return A list of data frames containing the input information, split into new groups.
#' @export
splitResultsByContigGroup<-function (df, contigColumn="chr", contigGroupsFile) {
    #df=data.frame(chr=c("1", "2", "3", "X", "MT"), pvalue=rnorm(5));contigColumn="chr";groupName="chrX"
    contigGroups=getContigGroups(contigGroupsFile)

    getOne<-function (groupName, contigGroups, df) {
        cg=contigGroups[[groupName]]
        idxResult=sort(which(!is.na(match(df[[contigColumn]], cg))))
        df[idxResult,]
    }

    r=lapply(names(contigGroups), getOne, contigGroups, df)
    names(r)=names(contigGroups)
    return (r)
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
#' @export
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
#' Get the residual expression after removing the effects of one or more covariates
#'
#' Fit a linear regression to explain how much of expression variance is explained by a covariate
#' IE: lm (expression ~ covariate)
#' Extract the residuals from the regression and add the mean expression.
#'
#' @param covarNames The covariates to subtract from the expression data
#' @param expData The expression data
#' @param covars The covariate data
#' @return A data frame in the same format as the input expression data, with the effects of the covariate removed.
#' @export
removeCovarEffectOnExpression<-function (covarNames, expData, covars) {
    if (dim (expData)[2]!=dim (covars)[2])
        stop("Number of donors in expression data and covariates must be the same.")

    expDataT=t(expData[,-1]); colnames(expDataT)=expData[,1]
    covarsList=lapply(covars$id, function (covarName) {as.numeric(covars[covars$id==covarName,-1])})
    names(covarsList)=covars$id
    covarsList=covarsList[covarNames]

    strTerms=paste(sapply(1:length(covarsList), function (x) paste("covarsList[[",x,"]]", sep="")), collapse=" + ")

    formStr=paste("expDataT ~ ", strTerms)
    z=lm(formula(formStr))

    #extract out the residuals, add the mean expression.
    d=residuals(z)
    dd=t(d)
    expMeans=apply(expData[,-1], 1, mean)
    expNew=dd+expMeans

    #make this look like the original expression data
    expNew=data.frame(id=rownames(expNew), expNew, stringsAsFactors = F, check.names=F)
    rownames (expNew)=NULL
    return (expNew)
}

#

#' Renormalize expression data from gene/sum(all genes) to gene/sum(autosome genes).
#'
#' Renormalize expression of the autosomes to sum to 1, instead of including the X/Y/MT chromosomes.
#'
#' @param expData The expression data.  This must be a data frame, NOT a data table.
#' @param geneLoc The location of each gene.  Encodes the gene name, chromosome, start, end. (optional, can use this or reducedGTF)
#' @param reducedGTF The location of each gene.
#' @param contigGroupsFile A yaml file containing the contig groups
#'
#' @return The expression data after normalization
#' @export
renormalizeAutosomeExpression<-function (expData, geneLoc=NULL, reducedGTF=NULL, contigGroupsFile=NULL) {
    if (is.null(contigGroupsFile)) return (expData)
    if (is.null(geneLoc) && is.null(reducedGTF)) return (expData)

    if (!"id" %in% colnames(expData))
        stop ("Expected id column in expression data!")

    contigGroups=getContigGroups(contigGroupsFile)
    if (!is.null(geneLoc))
        genesAutosome=geneLoc[which(!is.na(match(geneLoc$chr, contigGroups$autosomes))),]$geneid
    if (!is.null((reducedGTF)))
        genesAutosome=reducedGTF[which(!is.na(match(reducedGTF$chr, contigGroups$autosomes))),]$gene_name

    #genes can be in the annotations but not in the expression data.
    idxExpData=which(!is.na(match(expData$id, genesAutosome)))

    #renormalize autosomes.
    expDataAutosomes=expData[idxExpData,-1]
    beforeColSums=colSums(expDataAutosomes)
    #plot (beforeColSums, covars[covars$id=="FRACTION_X",-1], xlab="data before normalization", ylab="fraction expressed X", main="Before normalization to autosome genes")
    expDataAutosomes=sweep(expDataAutosomes,2,colSums(expDataAutosomes),'/')
    expDataAutosomes=expDataAutosomes*100000
    afterColSums=colSums(expDataAutosomes)
    #plot (afterColSums, covars[covars$id=="FRACTION_X",-1], xlab="data after normalization", ylab="fraction expressed X", main="After normalization to autosome genes")

    #plug data back in.
    expData[idxExpData,-1]=expDataAutosomes

    return (expData)
}




##################################################
# eQTL / Permutations Math bits
##################################################


#' Run matrix eQTL on an expression data as the dependent variable to measure the effects of an independent variable encoded by a covariate.
#'
#' @param covarName The covariate name to analyze
#' @param covars A data frame of 1 or more covariates.  Each row is a covariate, each column a sample.  The first column "id" contains the name of the covariate.
#' @param expData A data frame of expression data, 1 row per gene, 1 column per sample. The first column "id" contains the name of the gene.
#' @param numPermutations Runs permutations of the linear regression, permuting the covariate donor labels.
#'
#' @return The empiric results if numPermutations=0, else the permuted results.
#' @export
runOneCovarMatrixeQTL<-function (covarName, covars, expData, numPermutations=1000) {
    ee=DropSeq.eqtl::convertToSlicedData(expData)
    phenotype=DropSeq.eqtl::convertToSlicedData(covars[covars$id==covarName,])

    #run the single covar against the normalized expression.
    resultSingle=runSingleCovarFast(ee,phenotype, cvrt=MatrixEQTL::SlicedData$new())
    #run permutations if requested
    if (numPermutations==0) return (resultSingle)
    resultP=runSingleCovarPermutations (resultSingle, ee, phenotype, cvrt=MatrixEQTL::SlicedData$new(), numPermutations=numPermutations)
    return (resultP)
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
runSingleCovarPermutations<-function (resultFast, ee, phenotype, cvrt, numPermutations=1000, random.seed=1) {
    set.seed(random.seed)
    resultFast$exceeds_empiric_p=0
    for (i in 1:numPermutations) {
        if (i%%100==0)cat (paste(Sys.time(), " Permutation [", i, "] of [", numPermutations, "]\n"))
        #permute phenotype
        phenotype$ColumnSubsample(sample (1:phenotype$nCols()))
        permResult=runSingleCovarFast(ee, phenotype, cvrt)
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

#' Quantile normalize expression data
#' The goal of this is to normally distribute expression data that may not be normally distributed.
#' See: https://en.wikipedia.org/wiki/Quantile_normalization
#'
#' @param expressionData The data.frame or data.table of expression data to normalize.  Rownames contain genes, column
#' names contain donor identifiers.  The values of the data.frame are all numeric.
#' @return A data.table of expression data with the same number of rows and columns as the input data, but with quantile normalized expression values
#' @import data.table
#' @export
quantileNormalizeExpression<-function (expressionData) {
    ranks=apply(expressionData, 2, rank, ties.method = "min")
    sortedMat=apply(expressionData,2,function (x) sort(x, decreasing=F))
    rankValues=apply(sortedMat, 1, mean)
    r=apply(ranks, c(1,2), function (x) rankValues[x])
    r=data.table::data.table(r)
    rownames(r)=rownames(expressionData)
    return (r)
}


