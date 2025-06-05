# baselineFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.baseline.summary.txt"
# #testFile="/broad/mccarroll/nemesh/single_cell_eQTL/old_code_results/BA46.astrocyte.noOutliers.single_cell_eQTL_interaction.summary.txt"
# testFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k100.factor2.summary.txt"
# #testFile="/broad/mccarroll/nemesh/single_cell_eQTL/compressed_data/BA46_Astrocyte.single_cell_eQTLs.k20.factor_2.dynamic_eQTL.summary.txt"
# #testFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k20.all_factor.summary.txt"
# #testFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k20.factor2.qe_05.summary.txt"
# nCellsTested=100

#evaluateVectorQuant(baselineFile, testFile, 20)
evaluateVectorQuant<-function (baselineFile, testFile, nCellsTested) {
    plotOne<-function (z, nCellsTested, lim=NULL) {
        #compute the sum squared error
        sse=sum ((((-log10(z$baseline))- (-log10(z$test))))^2)
        absDif=sum (abs(((-log10(z$baseline))- (-log10(z$test)))))
        plot (-log10(z$baseline), -log10(z$test), xlim=lim, ylim=lim, xlab="baseline", ylab=paste(nCellsTested, "cells per donor"), col="light green", pch=16)
        points (-log10(z[z$status=="FP",]$baseline), -log10(z[z$status=="FP",]$test), col="red", pch=16)
        points (-log10(z[z$status=="FN",]$baseline), -log10(z[z$status=="FN",]$test), col="blue", pch=16)
        points (-log10(z[z$status=="TP",]$baseline), -log10(z[z$status=="TP",]$test), col="dark green", pch=16)
        title (paste("Clusters per donor [", nCellsTested, "] SSE [", round (sse,1), "] abs diff [", round (absDif, 1), "]\n", "pvalue of genotype: factor interaction"))
        abline (0,1, col='black')
        legendStr=c(
            paste("True Positive [", length(which(z$status=="TP")), "]", sep=""),
            paste("True Negative [", length(which(z$status=="TN")), "]", sep=""),
            paste("False Positive [", length(which(z$status=="FP")), "]", sep=""),
            paste("False Negative [", length(which(z$status=="FN")), "]", sep=""))
        legend("topleft", legend = legendStr, fill=c("dark green", "light green", "red", "blue"))
    }

    a=read.table(baselineFile, header=T, stringsAsFactors = F, sep="\t")
    b=read.table(testFile, header=T, stringsAsFactors = F, sep="\t")
    both=intersect(a$gene, b$gene)
    a=a[match(both, a$gene),]
    b=b[match(both, b$gene),]

    if (any(a$gene!=b$gene))
        stop ("Gene names don't match between data sets")

    pvalueThreshold=0.05
    z=data.frame(baseline=a$lrt_pval, test=b$lrt_pval, status=NA)
    z[which(z$baseline<pvalueThreshold & z$test<pvalueThreshold),]$status<-"TP"
    z[which(z$baseline>pvalueThreshold & z$test>pvalueThreshold),]$status<-"TN"
    z[which(z$baseline>pvalueThreshold & z$test<pvalueThreshold),]$status<-"FP"
    z[which(z$baseline<pvalueThreshold & z$test>pvalueThreshold),]$status<-"FN"

    plotOne(z, nCellsTested, lim=c(0, 20))

}

plotROCVectorQuant<-function () {
    s=list.files("/downloads/single_cell_eQTL", pattern="summary", full.names = T)
    testIdx=grep (".k", s)
    baseline=s[-testIdx]
    tests=s[testIdx]

    if (length(baseline)>1)
        stop ("Multiple files could be the baseline?")

    getLabel<-Vectorize(function (x) { if (x<0.05) return (0) else return (1)})

    #test=tests[3]
    plotOne<-function (test, baseline) {
        z=strsplit (test, ".", fixed=T)[[1]]
        k=z[length(z)-2]
        k=as.numeric (sub("k", "", k))

        a=read.table(baseline, header=T, stringsAsFactors = F, sep="")
        b=read.table(test, header=T, stringsAsFactors = F, sep="")
        both=intersect(a$gene, b$gene)
        a=a[match(both, a$gene),]
        b=b[match(both, b$gene),]

        df=data.frame(predictions=b$lrt_pval, labels=getLabel(a$lrt_pval))
        pred <- prediction(df$predictions, df$labels)
        perf <- performance(pred,"tpr","fpr")
        auc=performance(pred, measure = "auc")
        auc=round (auc@y.values[[1]],3)

        plot(perf,colorize=TRUE, main=paste("K [", k, "] AUC [", auc, "]"))
        df2=data.frame(baseline=a$lrt_pval<0.05, compressed=b$lrt_pval<0.05)
        t=table (df2)
    }



}


compareMetricsDF<-function () {
    o=read.table("/downloads/single_cell_eQTL/old_code/metricsDF.txt", header=T, stringsAsFactors = F, sep="\t")
    n=read.table("/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k20.cell_features.txt", header=T, stringsAsFactors = F, sep="\t")

    setdiff(o$CELL_BARCODE_FINAL, n$CELL_BARCODE_FINAL)

    colsBoth=intersect (colnames(o), colnames(n))
    o=o[,colsBoth]
    n=n[,colsBoth]

    all.equal(o,n)
}

compareFactorScore<-function () {
    o=read.table("/downloads/single_cell_eQTL/old_code/metricsDF.txt", header=T, stringsAsFactors = F, sep="\t")
    n=read.table("/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k20.usages.txt", header=T, stringsAsFactors = F, sep="\t", row.names=1)
    all (rownames (n)==o$CELL_BARCODE_FINAL)


}

compareExpressionDF<-function () {
    oldExp="/downloads/single_cell_eQTL/old_code/expression.txt"
    h5SeuratFile="/downloads/single_cell_eQTL/BA46_Astrocyte.single_cell_eQTLs.k20.h5seurat"

    data <- SeuratDisk::LoadH5Seurat(h5SeuratFile, assays = "RNA", verbose =TRUE, mode="r")
    n=Seurat::GetAssayData(data)
    o=read.table(oldExp, header=T, stringsAsFactors=F, sep="\t", check.names=F)

    bothGenes=intersect(rownames (o), rownames(n))
    o=o[bothGenes,]
    n=n[bothGenes,]

    bothCols=intersect(colnames(o), colnames(n))
    o=o[,bothCols]
    n=n[,bothCols]

    n=as.matrix(n)



}

