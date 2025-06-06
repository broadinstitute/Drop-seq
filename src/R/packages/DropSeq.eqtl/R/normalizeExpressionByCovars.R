
expressionFile="/downloads/eQTL_Covars/d5QueenBasic.maf_0.20_cisDist_10kb.gene_expression.txt"
covariateFile="/downloads/eQTL_Covars/PEER/d5QueenBasic.full_covars.peer_enhanced.txt"
outExpressionFile="/downloads/eQTL_Covars/PEER/d5QueenBasic.expression_residualized.txt"

peerNormalizedExpFile="/downloads/eQTL_Covars/PEER/d5QueenBasic.maf_0.20_cisDist_10kb.gene_expression.peer_residuals.txt"

#' Remove the effect of covariates from the expression data
#'
#' Given the expression and covariate inputs, regress out the covariates from the expression data
#' and generate new expression data. This is done by running a linear regression of
#' expression ~ covars, extracting the residuals and then adding the mean expression to the residuals.
#'
#' @param expressionFile The input expression data set - columns: id column for gene followed by donor names.
#' rows: gene name followed by expression level per donor
#' @param covariateFile The input covariate data set - columns: id column for covariate followed by donor names.
#' rows: covariate name followed by covariate values per donor.
#' @param outExpressionFile The ouput file, in the same format as the input.
#' @return If outExpressionFile is NULL, return the modified expression data frame.
#' @export
normalizeExpressionByCovariates<-function (expressionFile, covariateFile, outExpressionFile=NULL) {
    exp=read.table(expressionFile, header=T, stringsAsFactors = F, sep="\t", check.names=F)
    covs=read.table(covariateFile, header=T, stringsAsFactors = F, sep="\t", check.names = F)
    covs=covs[1:17,] #get rid of peer factors for now.


    z=lm(as.numeric(exp[1,-1]) ~ t(as.matrix(covs[,-1])))
    rr=residuals(z)

    o=read.table("/downloads/eQTL_Covars/PEER/residuals.txt",sep="\t")
    plot (rr, o[,1], xlab="Jim Residuals", ylab="PEER residuals", main="Known covariates residuals")

    #set up covariates
    rownames (covs)=covs$id
    covs2=t(covs[,-1])

    #set up expression
    rownames(exp)=exp$id
    exp2=t(as.matrix (exp[,-1]))

    z=lm(exp2 ~ covs2)
    r=residuals(z)

    #how correlated are the results?
    corR=sapply(1:dim(r)[2], function (idx) cor (r[,idx], o[,idx]))
    hist (corR, breaks=100, main="Correlation of residuals between Jim and Peer")


    #add mean expression to residuals
    meanExp=apply (exp2, 2, mean)
    expResidualized=t(r)+meanExp

    result=data.frame(id=rownames(expResidualized), expResidualized, stringsAsFactors = F, check.names = F)
    rownames(result)=NULL
    if (!is.null(outExpressionFile)) write.table(result, outExpressionFile, row.names = F, col.names = T, quote=F, sep="\t")
    if (is.null(outExpressionFile)) return (result)


}

testVsPeer<-function () {
    exp=read.table(expressionFile, header=T, stringsAsFactors = F, sep="\t", check.names=F)
    d=normalizeExpressionByCovariates(expressionFile, covariateFile, outExpressionFile=NULL)
    norm=read.table(peerNormalizedExpFile, header=T, stringsAsFactors = F, sep="\t", check.names=F)

    z1=as.numeric(d[1,-1])
    z2=as.numeric(norm[1,-1])

    z3=as.numeric (exp[1,-1])

    plot (z3, z1, xlab="raw expression", ylab="jim normalized expression")
    plot (z3, z2, xlab="raw expression", ylab="PEER normalized expression")

}
