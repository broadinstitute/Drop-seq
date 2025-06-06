#general params

# method="linear";qqPlot=T; filterExpressionToDonorList=T; quantileNormalize=F
# autosomesOnly=F; normalizeAutosomeExpression=T;pctGenesRetained=50;logExpression=T

######################
# Expression or metaCellFile
######################
# metaCellFile="/downloads/group_differences/d5QueenBasic.meta_cell.expression.txt"
# rejectedDonorFile="/downloads/group_differences/rejectedDonors.txt"
#
# donorListOneFile=NULL; donorListOne=NULL;donorListTwo=NULL;donorListTwoFile=NULL
#
# contigGroupsFile="/broad/mccarroll/software/metadata/individual_reference/GRCh38_maskedAlt.89/GRCh38_maskedAlt.contig_groups.yaml"
# reducedGTFFile="/broad/mccarroll/software/metadata/individual_reference/GRCh38_maskedAlt.89/GRCh38_maskedAlt.reduced.gtf"

#####################################################################################
# Wilcox test, use donor files to define donors, no covariate correction.
######################################################################################
# covariatesFile=NULL;residualizeCovariates=c();testCovar=NULL
# donorListOneFile="/downloads/group_differences/QB_FB.txt"


##############################
# Wilcox Test setup with covariates
##############################
# covariatesFile="/broad/mccarroll/dropulation_census_metadata/COVARIATES/CIRMw1w2w3.covariates_plusCollectionSite.txt"
# testCovar="CELL_TYPE"
# residualizeCovariates=c("SEX", "PC1", "PC2", "PC3", "PC4")
#
# testCovar="SEX"
# residualizeCovariates=c("CELL_TYPE", "PC1", "PC2", "PC3", "PC4")


##############################
# Linear regression setup
##############################
# testCovariates=c("CELL_TYPE", "SEX", "PC1", "PC2", "PC3", "PC4")


#' Given a set of covariates, remove the effects of covariates on gene expression before applying group level difference test.
#'
#' This loads in a covariates file containing one more more covariates, then runs a linear regression using the expression as the dependent variable
#' and the selected covariates as independent variables.  From these regressions, the residuals are extracted and added to the median
#' expression of the gene.  This new gene expression accounts for the covariates, and is then used in the groupLevelDifferences test.
#' The gene level differences are calculated using the wilcoxon rank sum test to compare groups.
#' @param covariatesFile A file containing covariates.  The first column is named IID and contains the donor IDs.  Columns after the first contain covariates.
#' Each row contains the donor ID, followed by the values of the covariates.  If null, then no covariates will be accounted for in the analysis, and
#' the method will be unable to partition donors into groups by a covariate.
#' @param testCovar A covariate used to partition the data into two groups.  This covariate must split the data into exactly two groups,
#' for example 0/1, A/B, TRUE/FALSE, etc.  This is optional, and if null this method will fall back on the donorList parameters.
#' @param residualizeCovariates A list of covariates from the covariatesFile which will be accounted for in the regression and subtracted out of
#' the expression data.
#' @param donorListOne A vector containing a list of donor IDs, which must be in the column names of the expressionFile.
#' This will select the donors in the first group.
#' @param donorListTwo An optional vector containing a list of donor IDs for the second group.  If NULL, assume all donors
#' in the expression file not in donorListOne are in this list.
#' @param donorListOneFile  If set, reads in the donorListOne from this file.  Overrides settings to donorListOne.
#' @param donorListTwoFile  If set, reads in the donorListOne from this file.  Overrides settings to donorListTwo.
#' @param reducedGTFFile The location of each gene.  Encodes the gene name, chromosome, start, end.
#' @param contigGroupsFile A yaml file containing the contig groups.
#' @param logExpression If true, use the log(expression) in the linear regression to remove the effects of confounding covariates.
#' @param quantileNormalize If true data is quantile normalized.  This is probably a terrible idea.
#' @param method The linear regression model used to remove  effects of covariates on the expression data..  Model "linear" uses the standard regression model.
#' Model "robustLinear" uses MASS:rlm to find the best weights for each donor in the regression and coefficient for the model that is less influenced by outliers.
#' @param outFile writes the results data frame to this file.
#' @param outPDF If not null, emit plots to this file.
#' @return A data frame containing the gene name, p-value for the group level test, and permuted p-value if requested.
#' @inheritParams readGroupLevelDifferencesData
#' @export
groupLevelDifferencesWithCovariates<-function (metaCellFile, rejectedDonorFile=NULL, covariatesFile=NULL, testCovar=NULL, residualizeCovariates=c(), contigGroupsFile,
                                               reducedGTFFile, donorListOne=NULL, donorListTwo=NULL, donorListOneFile=NULL, donorListTwoFile=NULL,
                                               autosomesOnly=F, normalizeAutosomeExpression=T, logExpression=T, pctGenesRetained=50,
											   quantileNormalize=F, method=c("linear", "robustLinear"), outFile=NULL, outPDF=NULL) {
    method=match.arg(method)

    #optional to evalute covarites.
    #This should default the straight up wilcox test.
    if (!is.null(covariatesFile)) {
        #if (length(residualizeCovariates)==0) stop ("Must supply 1 or more residualizeCovariates to remove from expression data!")

        covarsDF=readCovariateFile(covariatesFile)
        groups=getDonorListsFromCovariates(covarsDF, testCovar)
        donorListOne=groups$donorListOne
        donorListTwo=groups$donorListTwo
    } else {
        covarsDF=NULL
    }

    #debug (readGroupLevelDifferencesData)
	r=readGroupLevelDifferencesData(metaCellFile, rejectedDonorFile, donorListOne, donorListTwo, donorListOneFile, donorListTwoFile,
	                                reducedGTFFile, contigGroupsFile, autosomesOnly=autosomesOnly, normalizeAutosomeExpression=normalizeAutosomeExpression,
	                                pctGenesRetained=pctGenesRetained, validateDonorsInExpData=F, filterExpressionToDonorList=T)
	expressionData=r$expressionData
	donorListOne=r$donorListOne
	donorListTwo=r$donorListTwo
	reducedGTF=r$reducedGTF

	if (quantileNormalize) expressionData=DropSeq.eqtl::quantileNormalizeExpression(expressionData)

	#get covariates and expression data ordered the same way with the same donors.
	if (!is.null(covarsDF) ) {
	    covarsDF=covarsDF[match(colnames(expressionData), covarsDF$IID),]
	    #optionally take the log of the expression data for the linear fitting.
	    system.time(expressionFitted<-getFittedExpressionLinear(expressionData, residualizeCovariates, covarsDF, method=method, logExpression=logExpression))
	}  else {
	    expressionFitted=expressionData
	}


    system.time(g<-getGroupLevelDifferencesFaster(expressionFitted, donorListOne, donorListTwo))
    g=addMedianExpressionThis(expressionData, g)
    g=addGeneLocations(reducedGTF, g, geneColName="gene")
    #reorder columns of g
    g=g[,c("gene", "chr", "start", "end", "medianExpression", "pvalue", 'FDR_pvalue', "meanG1", "meanG2", "effectSize")]
    g$group_one_label=covarsDF[covarsDF$IID==donorListOne[1],][[testCovar]]
    g$group_two_label=covarsDF[covarsDF$IID==donorListTwo[1],][[testCovar]]

    if (!is.null(outPDF)) {
        pdf(outPDF)
        contigGroups=DropSeq.eqtl::getContigGroups(contigGroupsFile)
        strTitle=paste("Wilcox Rank Sum Test")
        if (!is.null(testCovar)) strTitle=paste(strTitle, "[", testCovar, "]")
        strTitle=paste(strTitle, "\nmethod [", method, "] covariates [", paste(residualizeCovariates, collapse=","), "]", sep="")
        if (quantileNormalize) strTitle=paste(strTitle, "\nquantile normalized")
        qqplotSimple(g$pvalue, main=paste(strTitle))
        qqplotSimple(g[g$chr!=contigGroups$chrX,]$pvalue, main=paste(strTitle, "\nautosomes only"))
        plot (g$effectSize, -log10(g$pvalue), xlab="effect size [mean(group1)-mean(group2)/sd(all)", ylab="p-value wilcoxon rank sum test [-log10]")
        dev.off()
    }

    if (!is.null(outFile)) {
        write.table(g, outFile, row.names=F, col.names = T, quote=F, sep="\t")
    }
    #if outFile is null return
    return (g)
}

#' Simultaneously test many covariates for their effects on expression
#'
#' This tests each gene for differences in expression across multiple groups simultaniously.  Each gene's expression is the dependent variable
#' in a linear regression against all listed test covariates which are treated as independent variables.
#'
#' Note that donor selection criteria for list one and two sub-select from the entire list, but there is no explicit donor grouping
#' in this operation, as many independent variables aree tested simultaniously.  If all donor options are NULL (the default) then all
#' donors in the intersect of the covariates and expression data are used.
#'
#' @param testCovariates A list of covariates that will be tested for differential expression.
#' @param method The linear regression model used find effects.  Model "linear" uses the standard regression model.
#' Model "robustLinear" uses MASS:rlm to find the best weights for each donor in the regression, then runs a standard
#' regression with lm using those weights.
#' @inheritParams groupLevelDifferencesWithCovariates
#' @return A data frame containing the gene name, beta, and p-value for each tested covariate
#' @export
groupLevelDifferencesWithRobustRegression<-function (metaCellFile, rejectedDonorFile=NULL, covariatesFile, testCovariates=NULL, reducedGTFFile, contigGroupsFile,
                                                     donorListOne=NULL, donorListTwo=NULL, donorListOneFile=NULL, donorListTwoFile=NULL,
                                                     autosomesOnly=T, normalizeAutosomeExpression=F, pctGenesRetained=50,
                                                     quantileNormalize=F, method=c("linear", "robustLinear"), outFile=NULL, outPDF=NULL) {

    method=match.arg(method)

    covarsDF=readCovariateFile(covariatesFile)

    #r=readGroupLevelDifferencesData(metaCellFile, donorListOne, donorListTwo, donorListOneFile, donorListTwoFile, reducedGTFFile, contigGroupsFile,validateDonorsInExpData=F, filterExpressionToDonorList=T)
    r=readGroupLevelDifferencesData(metaCellFile, rejectedDonorFile, donorListOne, donorListTwo, donorListOneFile, donorListTwoFile,
                                    reducedGTFFile, contigGroupsFile, autosomesOnly=autosomesOnly, normalizeAutosomeExpression=normalizeAutosomeExpression,
                                    pctGenesRetained=pctGenesRetained, validateDonorsInExpData=F, filterExpressionToDonorList=T)

    expressionData=r$expressionData
    reducedGTF=r$reducedGTF

    if (quantileNormalize) expressionData=DropSeq.eqtl::quantileNormalizeExpression(expressionData)

    #get covariates and expression data ordered the same way with the same donors, slim down to the covars we want to test.
    covarsDF=covarsDF[match(colnames(expressionData), covarsDF$IID),c("IID", testCovariates)]

    g=groupDifferencesByLinearRegression (expressionData, testCovariates, covarsDF, method=method)
    g=addMedianExpressionThis(expressionData, g)
    g=addGeneLocations(reducedGTF, g, geneColName="gene")

    #clean up and reorder output column names.
    g=g[,c("covariate", "gene", "chr", "start", "end", "medianExpression", "beta","pvalue")]

    if (!is.null(outPDF)) {
        pdf(outPDF)
        #covName=unique (g$covariate)[1]
        for (covName in unique (g$covariate)) {
            cat (covName,"\n")
            strTitle=paste("Regression of group [", covName, "] \n method [", method,"] covariates [", paste(testCovariates, collapse=","), "]", sep="")
            if (quantileNormalize) strTitle=paste(strTitle, "\nquantile normalized")
            qqplotSimple(g[g$covariate==covName,]$pvalue, main=paste(strTitle))
        }
        dev.off()
    }

    if (!is.null(outFile)) {
        write.table(g, outFile, row.names=F, col.names = T, quote=F, sep="\t")
        return (NULL)
    }
    #if outFile is null return
    return (g)
}

#' Calculate the residuals of the fit of the expression data to the testCovariates.
#'
#' For each gene, run a linear regression with expression as the dependent variable and the listed testCovariates as the independent variables.
#' Extract the residuals, then add the sum of the median expression to the residuals.
#' Emit a new matrix of expression with the same dimensions.
#'
#' @param expressionData a data table of the expression data, with the rownames encoding the gene names, and the column names the donor names.
#' @param residualizeCovariates which covariates will be fit as the independent variables in the regression.  If this argument is of length 0, this function is a no-op.
#' @param covarsDF The data.frame or data table containing the covariates.  This contains the identifier in the first column (IID),
#' followed by one covariate per column, with the covariate name encoded in the column names.  Each row is a donor.
#' @param method Which linear regression method to use.  "linear" uses the default stats::lm function to fit the data, while
#' robustLinear uses MASS::rlm to fit the data and less heavily weigh outliers in the expression data.
#' @param logExpression If true, take the log of the expression data.  Any donor with expression 0 for a gene is set to the minimum non-zero value.
#' @return A data.table with the same dimensions and encoding as the input expression data, with modified expression values.
getFittedExpressionLinear<-function (expressionData, residualizeCovariates, covarsDF, method=c("linear", "robustLinear"), logExpression=F) {
    #gn="BTNL8"; x=expressionData[which (rownames(expressionData)==gn),]

    # Maybe what I need to do for the log fit is:
    # 1) Fit to log.
    # 2) Get the predicted values (in log) and exponentiate
    # 3) Calculate residuals from the data in not-log space.
    # 4) How does that compare to the log of the residuals?
    getResidualsLinear<-function (x, form, covarsList, logExpression=F) {
        x=as.numeric (x)
        if (logExpression) {
            idx0=which(x<=0)
            if (length(idx0)>0) {
                minValue=min (x[-idx0])
                x[idx0]<-minValue
            }
            xx=log(x)
            d=c(x=list(xx), covarsList)
            fitlog=lm (form, data=d)
            #residuals hand calculated with intercept so the two models match better
            r=(x-exp(fitlog$fitted.values))+exp(fitlog$coefficients[1])
        } else {
            d=c(x=list(x), covarsList)
            fit=lm (form, data=d)
            r=residuals(fit)+fit$coefficients[1]
        }
        r=r+median(x)
        df=data.frame(rbind (r))
        return (df)
        #plot (covarsList[[1]], x, xlab="SEX", ylab="expression"); abline(fit); title (gn)
    }

    getResidualsRobustLinear<-function (x, form, covarList, logExpression=F) {
        if (logExpression)
            stop ("Log expression not supported in robust model")
        x=as.numeric (x)
        d=c(x=list(x), covarsList)
        fit2=MASS::rlm(form, data=d, maxit=50)
        z2=residuals(fit2)+median(x)
        df2=data.frame(rbind (z2))
        return (df2)
    }


    if (length(residualizeCovariates)==0) return (expressionData)
    method=match.arg(method)
    #get the covariates you want to account for
    covarsList=lapply(residualizeCovariates, function (covarName) {covarsDF[[covarName]]})
    names(covarsList)=residualizeCovariates

    strTerms=paste(sapply(residualizeCovariates, function (x) paste("covarsList[[\"",x,"\"]]", sep="")), collapse=" + ")

    form=formula(paste("x ~ ", strTerms))

    result<-switch(method,
        linear=expressionData[,getResidualsLinear(.SD, form, covarsList, logExpression=logExpression), by=list(rownames(expressionData))],
        robustLinear=expressionData[,getResidualsRobustLinear(.SD, form, covarsList,logExpression=logExpression), by=list(rownames(expressionData))],
        stop ("Invalid method for regression")
    )

    result[,rownames:=NULL]
    rownames(result)=rownames(expressionData)
    colnames(result)=colnames(expressionData)
    return (result)

    # quick validate
    # idx=123
    # fit=lm(as.numeric (expressionData[idx,]) ~ covarsDF[,residualizeCovariates[1]] + covarsDF[,residualizeCovariates[2]] +  covarsDF[,residualizeCovariates[3]])
    # res=residuals(fit)
    # as.numeric (res)==as.numeric(residuals[idx,])
}


#' Run a multivariate linear regression to explain gene expression by one or more covariates.
#'
#' For each gene, run a linear regression with expression as the dependent variable and the listed testCovariates as the independent variables.
#' Extract the slope and p-value of each covariate's fit peer gene
#'
#' @param expressionData a data table of the expression data, with the rownames encoding the gene names, and the column names the donor names.
#' @param testCovariates which covariates will be fit as the independent variables in the regression
#' @param covarsDF The data.frame or data table containing the covariates.  This contains the identifier in the first column (IID),
#' followed by one covariate per column, with the covariate name encoded in the column names.  Each row is a donor.
#' @param method Which linear regression method to use.  "linear" uses the default stats::lm function to fit the data, while
#' robustLinear uses MASS::rlm to fit the data and less heavily weigh outliers in the expression data.  These weights
#' are then used to fit the standard linear regression with those outliers providing less influence on the regression.
#' @return A data.table with the pvalue and beta per gene and testCovariate.
groupDifferencesByLinearRegression<-function (expressionData, testCovariates, covarsDF, method=c("linear", "robustLinear")) {

    getLabelsFromFit<-function (coef) {
        labelsRaw=strsplit (rownames(coef), '\"', fixed=T)
        l1=sapply(labelsRaw, function (x) x[2])
        l2=sub("]]", "", sapply(labelsRaw, function (x) x[3]))
        l2[l2==""]<-NA
        labelFinal=sub("_NA", "", paste (l1, l2, sep="_"))
        return (labelFinal)
    }

    extractResultsFromFit<-function (fit) {
        z=summary(fit)
        coef=z$coefficients
        coef=coef[-1,]
        labels=getLabelsFromFit(coef)

        m=data.frame(covariate=labels, beta=as.numeric (coef[,1]), pvalue=as.numeric (coef[,4]), stringsAsFactors = F)
        return (m)
        # getOne<-function (index, coef, labels) {
        #     z=coef[index,c(1,4)]
        #     return (z)
        # }
        # m=data.frame(matrix(sapply(1:dim (coef)[1], getOne, coef, labels), byrow=T, nrow=1), stringsAsFactors = F)
        # colnames(m)=as.character (matrix(sapply(labels, function (x) paste(x, c("beta", "pvalue"), sep="_")), byrow=T, nrow=1))
        # return (m)
    }

    #gn="AL627309.1"; x=expressionData[which (rownames(expressionData)==gn),]
    linearRegression<-function (x, form, covarsList) {
        d=c(x=list(as.numeric(x)), covarsList)
        fit=lm (form, data=d)
        extractResultsFromFit(fit)
    }

    robustLinearRegression<-function (x, form, covarList) {
    	w=NULL;
        d=c(x=list(as.numeric(x)), covarsList)
        fit2=MASS::rlm(form, data=d, maxit=50)
        d=c(d, w=list(fit2$w))
        fit3=lm(form, data=d, weights=w)
        extractResultsFromFit(fit3)

    }

    method=match.arg(method)
    #get the covariates you want to account for
    covarsList=lapply(testCovariates, function (covarName) {covarsDF[[covarName]]})
    names(covarsList)=testCovariates

    strTerms=paste(sapply(testCovariates, function (x) paste("covarsList[[\"",x,"\"]]", sep="")), collapse=" + ")

    form=formula(paste("x ~ ", strTerms))
    result<-switch(method,
                      linear=expressionData[,linearRegression(.SD, form, covarsList), by=list(rownames(expressionData))],
                      robustLinear=expressionData[,robustLinearRegression(.SD, form, covarsList), by=list(rownames(expressionData))],
                      stop ("Invalid method for regression")
    )

    colnames(result)[1]="gene"
    return (result)

}


#' Reads in the data neccesary to do group level tests
#'
#' @param metaCellFile Extract expression data from a metacell file.  For each cell, this normalizes expression to sum to 100,000 across all genes.
#' @param donorListOne A vector containing a list of donor IDs, which must be in the column names of the expressionFile.
#' This will select the donors in the first group.
#' @param donorListTwo An optional vector containing a list of donor IDs for the second group.  If NULL, assume all donors
#' in the expression file not in donorListOne are in this list.
#' @param donorListOneFile  If set, reads in the donorListOne from this file.  Overrides settings to donorListOne.
#' @param donorListTwoFile  If set, reads in the donorListOne from this file.  Overrides settings to donorListTwo.
#' @param reducedGTFFile a data frame with column names gene_name, chr, start, end
#' @param contigGroupsFile A yaml file containing the contig groups.
#' @param validateDonorsInExpData Test to see if all donors in the donor lists are in the metaCellFile
#' @param filterExpressionToDonorList If true, subsets the donors in the metaCellFile to the donors in the union of the donor lists.
#' @inheritParams readMetaCellData
#' @return A list with 4 elements: 1) A data frame of expression data 2) a data frame of gene locations
#' 3) A vector donor IDs in group 1. 4) A vector donor IDs in group 2.
#' @export
readGroupLevelDifferencesData<-function (metaCellFile, rejectedDonorFile=NULL, donorListOne, donorListTwo=NULL, donorListOneFile=NULL, donorListTwoFile=NULL,
                                             reducedGTFFile, contigGroupsFile,
                                             autosomesOnly=T, normalizeAutosomeExpression=F, pctGenesRetained=50,
    										 validateDonorsInExpData=T, filterExpressionToDonorList=T) {
	#input validation
	validateInputFilesExist(metaCellFile, donorListOneFile, donorListTwoFile)
	reducedGTF=readReducedGTF(reducedGTFFile)
	expressionData=readMetaCellData(metaCellFile, rejectedDonorFile, reducedGTF, contigGroupsFile, autosomesOnly=autosomesOnly,
	                                normalizeAutosomeExpression=normalizeAutosomeExpression, pctGenesRetained=pctGenesRetained)

	#read in and validate donor lists as neccesary
	if (!is.null(donorListOneFile)) donorListOne=read.table(donorListOneFile, header=F, stringsAsFactors = F)$V1
	if (!is.null(donorListTwoFile)) donorListTwo=read.table(donorListTwoFile, header=F, stringsAsFactors = F)$V1
	if (is.null(donorListTwo)) donorListTwo=setdiff(colnames (expressionData), donorListOne)

	#validate donor lists don't overlap, and are in the expression data.
	if (length(intersect (donorListOne, donorListTwo)) >0)
		stop(paste("Donor groups overlap for donors ["), paste(intersect (donorListOne, donorListTwo), collapse=","), "]")
	if (validateDonorsInExpData) {
		validateDonorsInExpressionData(expressionData, donorListOne)
		validateDonorsInExpressionData(expressionData, donorListTwo)
	}

	#restrict the donor lists to the expression data donors
	donorListOne=intersect(donorListOne, colnames(expressionData))
	donorListTwo=intersect(donorListTwo, colnames(expressionData))

	#restrict expression data to the donors in the donor lists
	if (filterExpressionToDonorList) {
	    allDonors=sort(c(donorListOne, donorListTwo))
	    geneNames=rownames(expressionData)
	    expressionData=subset(expressionData,,allDonors)
	    #for some reason data.table gets rid of rownames
	    rownames(expressionData)=geneNames
	}

	r=list(expressionData=expressionData, reducedGTF=reducedGTF, donorListOne=donorListOne, donorListTwo=donorListTwo)
	return (r)

}


plotManyPermutedDistribution<-function (expressionData, donorListOne, donorListTwo, numPermutations=100) {
	system.time (z<-replicate(numPermutations, generatePermutedResults(expressionData, donorListOne, donorListTwo)$pvalue))
	z=sort(as.vector (z))
	expected <- c(1:dim(expressionData)[1])/((dim(expressionData)[1])+1)
	expected=sort(rep(expected, numPermutations))
	maxXY=max(ceiling(max (c(-log10(expected), -log10(z)))))

	smoothScatter(-log10(expected), -log10(z), xlim=c(0, maxXY), ylim=c(0, maxXY),
				  xlab="expected pvalues [-log10] (uniform distribution)", ylab="observed permuted pvalues [-log10]")
	abline (0, 1, col="black")
	title (paste("Pvalue distribution for ", numPermutations, "permutations"))

	plot(-log10(expected), -log10(z), xlim=c(0, maxXY), ylim=c(0, maxXY), cex=0.25,
				  xlab="expected pvalues [-log10] (uniform distribution)", ylab="observed permuted pvalues [-log10]")
	abline (0, 1, col="red")
	title (paste("Pvalue distribution for ", numPermutations, "permutations"))

}

#' Add the median expression of the gene to the data frame.
#' Look up the gene name in the to-be-annotatated data frame using the geneColName
#' @param expressionData A matrix of expression data.  Each row is a gene, each column a donor.
#' Each rowname is a gene names.
#' @param df The data frame to annotate
#' @param geneColName The column name in df that contains the gene name.
addMedianExpressionThis<-function (expressionData, df, geneColName="gene") {
    #for R CMD CHECK
    medianExpression=NULL
	me=apply(expressionData, 1, median)
	idx=match(df[[geneColName]], rownames(expressionData))
	df[,medianExpression:=me[idx]]
	return (df)
}

#' Add the gene location of the gene to the data frame.
#' Look up the gene name in the to-be-annotatated data frame using the geneColName
#' @param reducedGTF A data.frame of gene locations.  Each row is a gene, columns are geneid (gene symbol), chr (contig), s1 (start), s2 (end).
#' If this parameter is null the operation is a no-op and the original dataframe is returned.
#' @param df The data table to annotate
#' @param geneColName The column name in df that contains the gene name.
addGeneLocations<-function (reducedGTF=NULL, df, geneColName="gene") {
	idx=match(df[[geneColName]], reducedGTF$gene_name)
	df[,c("chr", "start", "end"):=reducedGTF[idx,c("chr", "start", "end")]]
	return (df)
}

# An example of how the wilcoxon rank-sum works.
exemplarRankPlot<-function (geneName="ZNF667", expressionData, donorListOne, donorListTwo) {
	d=as.numeric(expressionData[which(rownames(expressionData)==geneName),])
	idxG1=match(donorListOne, colnames(expressionData))
	idxG2=match(donorListTwo, colnames(expressionData))
	ranks=rank(d)

	l=list(group1=d[idxG1], group2=d[idxG2])
	beeswarm::beeswarm(l, col=c("green", "blue"), pch=16, cex=0.75, ylab="gene expression")

}


#' Compare two sets of expression data using the wilcoxon rank sum test
#'
#' see #https://www.researchgate.net/post/How_can_I_calculate_the_effect_size_for_Wilcoxon_signed_rank_test
#' for discussion of effect size used here.
#' @param expressionData the matrix of expression data
#' @param donorListOne identifiers for group one
#' @param donorListTwo identifiers for group two
#' @return a data frame of genes with p-values and BH corrected FDR.
#' @export
#' @import matrixTests
#'
getGroupLevelDifferencesFaster<-function (expressionData, donorListOne, donorListTwo) {
	m1=as.matrix(expressionData[,donorListOne, with=F])
	m2=as.matrix(expressionData[,donorListTwo, with=F])
	r2=matrixTests::row_wilcoxon_twosample(m1, m2)
	df=data.frame(gene=rownames(expressionData), pvalue=r2$pvalue)
	df$FDR_pvalue=p.adjust(df$pvalue, method="BH")
	#note, the row_wilcoxon_twosample doesn't output the needed Z statistic for effect size.
    #but it can be calculated as qnorm(pvalue/2).  Z-score validated by coin::wilcox_test
	#zz=data.frame(group=factor(!is.na(match(colnames(expressionData), donorListOne))), expression=as.numeric (expressionData[5,]))
	#coin::wilcox_test(expression ~ group, data=zz)
	#zStat=sapply(r2$pvalue/2, qnorm)
	#effectSize
	#r<-z/sqrt(N)
	df$effectSizeR=sapply(r2$pvalue/2, qnorm)/sqrt(length(donorListOne)+length(donorListTwo))
    df$meanG1=apply(m1, 1, mean)
    df$meanG2=apply(m2, 1, mean)
    sd=apply(expressionData, 1, sd)
    df$effectSize=(df$meanG1-df$meanG2)/sd
	data.table::setDT(df)
	df=df[order(df$pvalue, decreasing=F),]
	return (df)
}

generatePermutedResults<-function (expressionData, donorListOne, donorListTwo) {
	n1=sample(colnames(expressionData), length(donorListOne))
	n2=sample (setdiff(colnames(expressionData), n1), length(donorListTwo))
	return (getGroupLevelDifferencesFaster(expressionData, n1, n2))
}

#' A simple QQ-plot
#' @param r A vector of p-values
#' @export
qqplotSimple<-function (r, ...) {
	observed <- sort(r)
	lobs <- -(log10(observed))
	expected <- c(1:length(observed))
	lexp <- -(log10(expected / (length(expected)+1)))

	maxXY=max(ceiling(max (c(lobs, lexp))))
	plot(c(0,maxXY), c(0,maxXY), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,maxXY), ylim=c(0,maxXY), las=1, xaxs="i", yaxs="i", bty="l", ...)
	points(lexp, lobs, pch=23, cex=.4, bg="black")

}

#' A plot of effect sizes
#'
#' @param g The results of the groupLevelDifferencesWithCovariates test
#' @param donorListOne the donors in group one.
#' @param covarsDF The first column is named IID and contains the donor IDs.  Columns after the first contain covariates.
#' @param testCovar The covariate tested by groupLevelDifferencesWithCovariates
#' @export
effectSizePlotSimple<-function (g, testCovar, ...) {
    #needs to extract who's in the first group to understand
    labels=names(table (covarsDF[[testCovar]]))
    groupOneLabel=unique (g$group_one_label)
    groupTwoLabel=unique (g$group_two_label)
    xLabel=paste("mean(", groupOneLabel, ") - mean(", groupTwoLabel, ") / sd (all)", sep="")
    xlim=c(-max(abs(g$effectSize)),max(abs(g$effectSize)))
    plot (g$effectSize, -log10(g$pvalue), xlab=paste("effect size", xLabel, sep="\n"), ylab="p-value wilcoxon rank sum test [-log10]", xlim=xlim, ...)
}


validateInputFilesExist<-function (metaCellFile, donorListOneFile=NULL,donorListTwoFile=NULL) {
    if (!file.exists(metaCellFile))
        stop (paste("MetaCell file does not exist [", metaCellFile, "]"))
	if (!is.null(donorListOneFile) && !file.exists(donorListOneFile))
		stop (paste ("Donor List File One does not exist! [", donorListOneFile, "]"))
	if (!is.null(donorListTwoFile) && !file.exists(donorListTwoFile))
		stop (paste ("Donor List File One does not exist! [", donorListTwoFile, "]"))

}

validateDonorsInExpressionData<-function (e, donorList) {
	if (is.null(donorList)) return()
	s=setdiff(donorList, colnames (e))
	if (length(s)>0)
		stop (paste("Requested donors in list that were not in expression data matrix! ["), paste(s, collapse=" "), "]")
}

#' Read in the meta cells file and normalize to look like the expression data.
#'
#' Convert the meta cell counts data into normalized expression data.
#' @param metaCellFile Extract expression data from a metacell file.  For each cell, this normalizes expression to sum to 100,000 across all genes.  The top 50%
#' of genes are selected for downstream analysis.
#' @param rejectedDonorFile A file containing a list of donors that should be excluded from the analysis. Single column with no header.(Optional)
#' @param reducedGTF a data frame with column names gene_name, chr, start, end
#' @param contigGroupsFile A yaml file containing the contig groups.
#' @param numTranscriptsPerDonor Normalize expression data such that each donor (column of the matrix) sums to this many transcripts
#' @param autosomesOnly  If true filter the metaCellMatrix to the autosome contigs before any normalization occurs
#' @param normalizeAutosomeExpression If true, normalize the expression of genes on the autosomes such that the total expression on the autosomes = numTranscriptsPerDonor.
#' Can not specify this option and autosomesOnly to both be true.
#' @param pctGenesRetained Calculate the median number of transcripts per gene, and select this percent of the highly expressed genes to analyze.
#' This selection takes place after normalization to numTranscriptsPerDonor.
#' @return a data frame of expression data, where each column is a donor and each row is a gene.
#' @export
readMetaCellData<-function (metaCellFile, rejectedDonorFile, reducedGTF, contigGroupsFile, numTranscriptsPerDonor=100000, autosomesOnly=T, normalizeAutosomeExpression=F, pctGenesRetained=50) {

    if (autosomesOnly && normalizeAutosomeExpression)
        stop ("Can't filter data to autosomes AND normalize data to autosomes.  Pick one or neither option.")

    if (is.null(reducedGTF) || is.null(contigGroupsFile))
        stop ("To normalize expression requires non-null arguments for reducedGTF and contigGroupsFile")

    metaCellDF=read.table(metaCellFile, header=T, stringsAsFactors = F, check.names=F)
    if (!is.null(rejectedDonorFile)) {
        r=read.table(rejectedDonorFile, header=F, stringsAsFactors = F)$V1
        idx=na.omit(match(r, colnames(metaCellDF)))
        if (length(idx)>0) metaCellDF=metaCellDF[,-idx]
    }

    if (autosomesOnly) {
        metaCellDF=filterMetaCellsToAutosome(metaCellDF, reducedGTF, contigGroupsFile)
    } else {
        metaCellDF=filterMetaCellsToAutosomePlusX(metaCellDF, reducedGTF, contigGroupsFile)
    }

    #normalize such that each donor has 100,000 transcripts.
    metaCellDFNorm=sweep(metaCellDF[,-1],2,colSums(metaCellDF[,-1]),`/`)*numTranscriptsPerDonor
    rownames(metaCellDFNorm)=metaCellDF$GENE

    #filter to the top 505 of expressed genes
    perGeneExpressionMedian=apply(metaCellDFNorm, 1, median)
    frac=1-(pctGenesRetained/100)
    threshold=quantile(perGeneExpressionMedian, probs=frac)[[1]]
    idx=which(perGeneExpressionMedian>=threshold)

    metaCellDFNorm=metaCellDFNorm[idx,]
    if (normalizeAutosomeExpression) {
        #need id as first column for genes.
        metaCellDFNormTemp=cbind(id=row.names(metaCellDFNorm), metaCellDFNorm)
        metaCellDFNorm=DropSeq.eqtl::renormalizeAutosomeExpression(expData=metaCellDFNormTemp, reducedGTF=reducedGTF, contigGroupsFile=contigGroupsFile)
    }
    data.table::setDT(metaCellDFNorm, keep.rownames = T)
    rownames(metaCellDFNorm)=metaCellDFNorm$rn
    rn=NULL; #R CMD CHECK
    metaCellDFNorm[,rn:=NULL]
    return (metaCellDFNorm)
}

#given a metacell data frame, filter it to the genes that are on the autosomes
filterMetaCellsToAutosome<-function (metaCellDF, reducedGTF, contigGroupsFile) {
    contigGroups=DropSeq.eqtl::getContigGroups(contigGroupsFile)
    genesAutosome=reducedGTF[which(!is.na(match(reducedGTF$chr, contigGroups$autosomes))),]$gene_name
    genesAutosome=intersect(genesAutosome, metaCellDF$GENE)
    result=metaCellDF[match(genesAutosome, metaCellDF$GENE),]
    result=result[order(result$GENE),]
    return (result)
}

filterMetaCellsToAutosomePlusX<-function (metaCellDF, reducedGTF, contigGroupsFile) {
    contigGroups=DropSeq.eqtl::getContigGroups(contigGroupsFile)
    genesAutosome=reducedGTF[which(!is.na(match(reducedGTF$chr, c(contigGroups$autosomes, contigGroups$chrX)))),]$gene_name
    genesAutosome=intersect(genesAutosome, metaCellDF$GENE)
    result=metaCellDF[match(genesAutosome, metaCellDF$GENE),]
    result=result[order(result$GENE),]
    return (result)
}


#' Reads in a reduced GTF file
#' @param reducedGTFFile The location of each gene.  Encodes the gene name, chromosome, start, end.
#' @return a data frame with column names gene_name, chr, start, end
readReducedGTF<-function (reducedGTFFile) {
    b=read.table(reducedGTFFile, header=T, stringsAsFactors = F)
    b=b[b$annotationType=="gene",]
    return (b)
}


plotUncorrectedVsCorrectedPValues<-function (g, g2, residualizeCovariates, testCovar) {
    strTitle=paste("Comparison of [" , testCovar, "] \nnormalized by [", paste(residualizeCovariates, collapse="+"), "]", sep="")
    xyLim=c(0, max(c(-log10(g$pvalue), -log10(g2$pvalue))))
    plot (-log10(g$pvalue), -log10(g2$pvalue), xlab="pvalues", ylab="expression regressed pvalues", xlim=xyLim, ylim=xyLim)
    abline (0,1, col='red')
    title (strTitle)
    qqplotSimple(g, main=testCovar)
    strTitle=paste(testCovar, " normalized by [", paste(residualizeCovariates, collapse="+"), "]", sep="")
    qqplotSimple(g2, main=strTitle)

}


#a nice pretty flat contingency table for the covariates of interest in the covariates file, optionally
#restricted by the donors in the expression data.

#' Generate a flat contingency table for the covariates of interest
#'
#' @param covariatesOfInterest A list of covariates to include from the covariates file in the contingency table
#' @inheritParams groupLevelDifferencesWithCovariates
#' @return a flat contingency table
#' @seealso stats::ftable
#' @export
getCovariateTable<-function (covariatesFile, covariatesOfInterest=c("SEX", "CELL_TYPE"), metaCellFile=NULL, rejectedDonorFile=NULL) {
    if (length(covariatesOfInterest)==0)
        stop ("Need at least one covariatesOfInterest!")

    a=readCovariateFile(covariatesFile)
    #if there's an expression file, filter down to that set.
    if (!is.null(metaCellFile)) {
        b=read.table(metaCellFile, header=T, stringsAsFactors = F, sep="\t", check.names = F, nrows = 1)
        idx=match(colnames (b)[-1], a$IID)
        a=a[idx,]
    }
    #if there's a rejectedDonorFile, remove those donors
    if (!is.null(rejectedDonorFile)) {
        r=read.table(rejectedDonorFile, header=F, stringsAsFactors = F)
        a=a[match(setdiff(a$IID, r$V1), a$IID),]
    }
    rownames(a)=a$IID
    a=a[,-1]
    idx=sort(match(covariatesOfInterest, colnames(a)))
    a=a[,idx,drop=F]

    if (length(covariatesOfInterest)==1) {
        z=table(a,dnn=covariatesOfInterest)
    } else {
        z=ftable(a)
    }

    return (z)
}

plotManyExemplars<-function (g2, expressionData,expressionFitted, donorListOne,donorListTwo) {
    geneNames=g2[g2$FDR_pvalue<0.05,]$gene
    pdf("/downloads/weird_genes.pdf")
    sapply(as.vector(geneNames), plotExpressionExemplar, expressionData,expressionFitted, donorListOne,donorListTwo)
    dev.off()

}

#expressionFitted=expressionFittedLinearRobust;geneName="ELP5"
#geneName="MMP20"
plotExpressionExemplar<-function (geneName="CFAP52", expressionData,expressionFitted, donorListOne,donorListTwo) {
    defaultExpG1=as.numeric (expressionData[rownames(expressionData)==geneName,match(donorListOne, colnames(expressionData)), with=F])
    defaultExpG2=as.numeric (expressionData[rownames(expressionData)==geneName,match(donorListTwo, colnames(expressionData)), with=F])
    normExpG1=as.numeric (expressionFitted[rownames(expressionFitted)==geneName,match(donorListOne, colnames(expressionFitted)), with=F])
    normExpG2=as.numeric (expressionFitted[rownames(expressionFitted)==geneName,match(donorListTwo, colnames(expressionFitted)), with=F])
    d=list(group1=defaultExpG1, group2=defaultExpG2, group1_norm=normExpG1, group2_norm=normExpG2)
    beeswarm::beeswarm (d)
    title (geneName)
    return (d)
}

#' Quantile normalize expression data
#' The goal of this is to normally distribute expression data that may not be normally distributed.
#' See: https://en.wikipedia.org/wiki/Quantile_normalization
#'
#' @param expressionData The data.frame or data.table of expression data to normalize.  Rownames contain genes, column
#' names contain donor identifiers.  The values of the data.frame are all numeric.
#' @return A data.table of expression data with the same number of rows and columns as the input data, but with quantile normalized expression values
quantileNormalizeExpression<-function (expressionData) {
    ranks=apply(expressionData, 2, rank, ties.method = "min")
    sortedMat=apply(expressionData,2,function (x) sort(x, decreasing=F))
    rankValues=apply(sortedMat, 1, mean)
    r=apply(ranks, c(1,2), function (x) rankValues[x])
    r=data.table(r)
    rownames(r)=rownames(expressionData)
    return (r)
}

#' For a given covariate in a covariate data frame that splits the donors into two groups, extract the donor IDs for each group.
#'
#' @param covarsDF The first column is named IID and contains the donor IDs.  Columns after the first contain covariates.
#' To successfully partition the donor IDs into two groups, this column must contain 2 unique values or an error is returned.
#' @param testCovar The covariate name to group donors with.
#' @return A list with 2 vectors of donor IDs, 1 per group.  If testCovar is null the returned list is NULL.
#' @export
getDonorListsFromCovariates<-function (covarsDF, testCovar=NULL) {
    if (is.null(testCovar)) return (NULL)

    z=table (covarsDF[[testCovar]])
    if (length(z)!=2) {
        print (z)
        stop (paste("Provided test covariate column [", testCovar, "] doesn't partition data into 2 groups", sep=""))
    }
    idxGroup1=which(covarsDF[[testCovar]]==names(which.max(z)))
    donorListOne=covarsDF[idxGroup1,]$IID
    donorListTwo=covarsDF[-idxGroup1,]$IID
    return (list(donorListOne=donorListOne, donorListTwo=donorListTwo))
}


readCovariateFile<-function (covariatesFile) {
    covarsDF=read.table(covariatesFile, header=T, stringsAsFactors = F, check.names = F, sep="\t")
    if (! "IID" %in% colnames(covarsDF)) stop ("covariate file must have a donor id column named IID")
    return (covarsDF)
}



