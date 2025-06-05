#https://towardsdatascience.com/introducing-distance-correlation-a-superior-correlation-metric-d569dc8900c7
# https://stats.stackexchange.com/questions/1580/regression-coefficients-that-flip-sign-after-including-other-predictors

#library(reshape2)
library (heatmaply)

# GOALS:
# Are genes significant on multiple individual factors?  Is that reflected in the multi-factor fit?

# DATA GENERATION:
# Run discovery on the uncompressed data for 8 gene:snp interactions.
# Run on each factor individually, as well as all factors together, then compare results.


# VISUALIZATION:
# For each gene, plot the individual factor pvalues, as well as the joint model pvalue.
# For each gene, plot the individual factor effect sizes, and the joint model effect sizes.


singleFits=list.files("/downloads/single_cell_eQTL/8genes", full=T, pattern="factor\\d*.summary.txt")
allFits=list.files("/downloads/single_cell_eQTL/8genes", full=T, pattern="factors.summary.txt")
pdfOut="/downloads/single_cell_eQTL/8genes/fit_plots.pdf"

plotRegressionResults<-function (allFits, singleFits) {

    readSingleFit<-function (x) {
        b=read.table(x, header=T, stringsAsFactors = F, sep="\t")
        coef=b[grep ("genotype:", b$Feature),]
        fm=b[b$Feature=="FULL_MODEL_LRT",]
        fm=fm[match(fm$gene, coef$gene),]
        coef$full_model_lrt=fm$pval
        return (coef)
    }

    a=read.table(allFits, header=T, stringsAsFactors = F, sep="\t")
    b=do.call(rbind, lapply(singleFits, readSingleFit))
    b=b[order(b$gene, b$snp, b$Feature),]

    genes=unique (a$gene)

    pdf(pdfOut)
    sapply(genes, plotGene, a, b)
    dev.off()
}

#gene="RERE"
plotGene<-function (gene, a, b) {
    aa=a[a$gene==gene,]
    bb=b[b$gene==gene,]
    featuresBoth=intersect(aa$Feature, bb$Feature)
    aa=aa[match(featuresBoth, aa$Feature),]
    bb=bb[match(featuresBoth, bb$Feature),]

    plotPvalues<-function (aa, bb) {
        featureName="pval"
        m=rbind(aa[[featureName]], bb[[featureName]])
        rownames (m)=c("all_factors", "single_factor")
        colnames(m)=aa$Feature
        colnames(m)=sub("genotype:FACTOR_", "", colnames(m))
        gene=unique (aa$gene)
        barplot(-log10(m), main=gene, ylab="p-value [-log10]", xlab='latent factor',beside = TRUE, col=c("red", "blue"), las=1)
        legend("topright", legend=c("all factors", "single factor fit"), fill=c("red", "blue"))
    }

    plotEstimate<-function (aa, bb) {
        featureName="Estimate"
        m=rbind(aa[[featureName]], bb[[featureName]])
        rownames (m)=c("all_factors", "single_factor")
        colnames(m)=aa$Feature
        colnames(m)=sub("genotype:FACTOR_", "", colnames(m))

        barplot(m, main='', ylab="coefficient", xlab='latent factor',beside = TRUE, col=c("red", "blue"), las=1)
        #legend("topright", legend=c("all factors", "single factor fit"), fill=c("red", "blue"))
    }

    par(mfrow=c(2,1), mar=c(3,4.5, 1, 2))
    plotPvalues(aa, bb)
    plotEstimate(aa, bb)
}

plotFactorCorrelation<-function () {
    m=read.table("/downloads/single_cell_eQTL/8genes/RERE_data_FACTOR_1:FACTOR_2:FACTOR_3:FACTOR_4:FACTOR_5:FACTOR_6:FACTOR_7:FACTOR_8:FACTOR_9:FACTOR_10.txt", header=T, stringsAsFactors = F, sep="\t")
    m=m[,grep ("FACTOR", colnames(m))]
    heatmaply_cor(abs(cor (m)), distfun="pearson", limits=c(0,1), colors = viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis"))

    par(mar=c(5,5,4,4))
    z=cor (m[,2], m[,3])
    strTitle=paste("Cell Latent Factor Score\ncor [", round (z, 3), "]", sep="")
    smoothScatter(m[,2], m[,3], xlab="Factor 2", ylab="Factor 3", main=strTitle, cex.axis=1.5, cex.lab=1.5)
}

plotLatentFactorScore<-function () {
    cellMetaDataFile="/broad/mccarroll/dropulation/analysis/latent_factor/metadata/BA46/BA46.n191/BA46_cell_level_metadata.txt"
    usageFile="/broad/mccarroll/dropulation/analysis/latent_factor/BA46_n180_Astrocytes/BA46_n180_Astrocytes/BA46_n180_Astrocytes.usages.k_11.dt_0_05.consensus.txt"

    cm=getSingleCellMetrics (cellMetaDataFile, donorList=NULL, quantizationErrorThreshold=NULL)
    u=read.table(usageFile, header=T, stringsAsFactors=F, sep="\t", row.names=1)

    b=intersect(cm$CELL_BARCODE_FINAL, rownames(u))
    cm=cm[match(b, cm$CELL_BARCODE_FINAL),]
    u=u[match(b, rownames(u)),]

    plotFactorVsLibrarySize<-function (factorIdx, cm, u) {
        z=cor (u[,factorIdx], cm$NUM_TRANSCRIPTS)
        titleStr=paste("Factor", factorIdx, "\ncor [", round(z,3), "]")
        smoothScatter (cm$NUM_TRANSCRIPTS, u[,factorIdx], xlab="library size", ylab="cNMF factor score", main=titleStr, cex.lab=1.25, cex.axis=1.25)
    }

    sapply(1:dim(u)[2], plotFactorVsLibrarySize, cm, u)

    factorList=sort(c(2,8,3,9,11))
    par(mfrow=c(2,3))
    sapply(factorList, plotFactorVsLibrarySize, cm, u)
    sapply(1, plotFactorVsLibrarySize, cm, u)
}

#validate that the inputs for a single gene are correct.
#this means that all features including the latent factor are the same in the single latent factor and reference multi-latent factor data.
# for gene RERE, this validates all data. [DONE]
checkInputs<-function () {
    #get the list of individual latent factor runs.
    otherFiles=list.files("/downloads/single_cell_eQTL/8genes/", pattern="RERE_data", full.names = T)
    otherFiles=otherFiles[-grep(":", otherFiles)]

    #get the joint latent factor run
    m=read.table("/downloads/single_cell_eQTL/8genes/RERE_data_FACTOR_1:FACTOR_2:FACTOR_3:FACTOR_4:FACTOR_5:FACTOR_6:FACTOR_7:FACTOR_8:FACTOR_9:FACTOR_10.txt", header=T, stringsAsFactors = F, sep="\t")

    #evaluate one latent factor
    testOne<-function (index, otherFiles, m) {
        o=read.table(otherFiles[index], header=T, stringsAsFactors = F, sep="\t")
        colsBoth=intersect(colnames(o), colnames(m))
        all.equal(m[,colsBoth], o[,colsBoth])
    }


    sapply(1:length(otherFiles), testOne, otherFiles, m)

}


