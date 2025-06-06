
#metaCellFile="/downloads/emi/hg19_Cumulative_NucleiDropulation_BA46_2019-10-10.meta_cell.expression.glutamatergic_neurons.txt"
#donorFile="/downloads/emi/22q11_Donor_IDs.txt"


#' Differential expression via wilcox test
#'
#' For two sets of donors, test if the two sets of donors are part of the same continuous distribution.
#' Expression is converted from a continuous variable to a ranked representation.
#' If the given set of donor IDs is ranked higher or lower than the other donors, this will result in a
#' pvalue closer to 0.  In other words, the wilcoxon test determines if the expression values of the
#' selected set of donors shifted compared to all other donors.
#'
#' @param metaCellFile A meta cell file, containing expression counts of donors, where each row contains data for a gene,
#' and each column the expression of a donor for that gene.
#' @param donorFile A file with no header containing the set of donors to test against all other donors in the data set.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less"
#' @return A data frame containing the gene names and thir wilcox.test p-value
#' @export
nonParametricDifferentialExpression<-function (metaCellFile, donorFile, alternative="two.sided") {
    m=readMetaCellFile(metaCellFile)
    d=read.table(donorFile, header=F, stringsAsFactors = F, check.names = F, sep="\t")$V1
    validateDonors(m, d)
    #the index of the donors to test in the metacells.
    idxDonors=match(d, colnames(m))
    #note: this transposed from the original input.
    #calculating the ranks probably isn't neccesary, but does break ties in non-deterministic way.
    #this may affect stability of the output for genes that have many ties (such as low expresssion genes)
    ranks=apply(m, 1, order, decreasing=T)
    pvalues=apply (ranks, 2, function (x, idxDonors) {stats::wilcox.test(x=x[idxDonors], y=x[-idxDonors], alternative=alternative)$p.value}, idxDonors)
    result=data.frame(gene=names(pvalues), pvalue=as.numeric(pvalues))
    result=result[order(result$pvalue, decreasing=F),]
    return (result)
}

validateDonors<-function (m, d) {
    missingDonors=setdiff(d, colnames(m))
    if (length(missingDonors)>0)
        stop(paste("Some requested donors not found in meta cell file [", paste(missingDonors, collapse=","), "]", sep=""))
}

readMetaCellFile<-function (metaCellFile) {
    m=read.table(metaCellFile, header=T, stringsAsFactors = F, check.names = F, sep="\t")
    rownames(m)=m[,1]
    m=m[,-1]
    return (m)
}
