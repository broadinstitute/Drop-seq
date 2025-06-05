#' Plot the cis eQTL region for a single gene.
#'
#' This is a plot of the cis eQTL pvalues and the genes.
#' @param geneName the gene to plot, must be in the cisEQTLs data frame gene column.
#' @param start The start position to plot from
#' @param end The end position to plot to.
#' @param cisEQTLs is a data frame that includes the following columns:  gene, chr, pos, p-value.
#' @param geneLoc is a data frame of gene locations.  Has the following columns: geneid, chr, s1, s2.
#' s1 is the start of the gene, and s2 is the end.
#' @param ratioSNPsToGenes what fraction of the figure should be used by the SNP part of the plot?  [0-1].
#' @param indexSNP A SNP selected from the eQTL data (SNP column) as the best SNP for the region,
#' which will be highlighted in red on the SNP plot.  (Optional)
#' @param mafColumn if provided, will color eQTL snps by minor allele frequency found in this column
#' @export
#' @import RColorBrewer
ploteQTLRegion<-function (geneName, start, end, cisEQTLs, geneLoc, ratioSNPsToGenes=0.6, indexSNP=NULL, mafColumn=NULL) {
    cisEQTL=cisEQTLs[cisEQTLs$gene==geneName,]
    if (dim (cisEQTLs)[1]==0) warning ("Could not find gene", geneName, "in cisEQTLs")
    strTitle=paste("Gene", geneName)

    if (!is.null(indexSNP)) {
        indexSNPName=paste(indexSNP$chr, indexSNP$pos, indexSNP$a1, indexSNP$a2, sep=":")
        strTitle=paste(strTitle, "Index SNP", indexSNPName)
    }
    def.par <- par(no.readonly = TRUE)
    par(mar=c(0,5,2,2), mgp=c(2,0.5,0))
    mat=matrix(c(1,2),nrow=2)
    layout(mat, heights = c(ratioSNPsToGenes, (1-ratioSNPsToGenes)), respect = FALSE)
    #layout.show(n=2)
    #plot (z$pos, -log10(z$`p-value`), xlab="position", ylab="Observed (-logP)", main=strTitle, axes=F)
    cols=rep("black", 4)
    mafBreaks=c(0.05, 0.1, 0.15, 0.2, 1)
    mafLabels=paste(mafBreaks[-length(mafBreaks)], mafBreaks[-1], sep="-")
    if (!is.null(mafColumn)) cols=brewer.pal("Blues", n=8)[c(2,4,6,8)]
    mafCuts=cut(cisEQTL[[mafColumn]], breaks=mafBreaks)

    plot (cisEQTL$snp_pos, -log10(cisEQTL$`p-value`), xlab="", ylab="Observed (-logP) eQTL", main=strTitle, axes=F, cex=0.35, xlim=c(start, end), col=cols[as.numeric(mafCuts)])
    axis(2)
    if (!is.null(mafColumn))
        legend("topright", legend=mafLabels, fill=cols, title="MAF")
    if (!is.null(indexSNP)) {
        idxSNPPos=indexSNP$pos
        ii=cisEQTL[cisEQTL$snp_pos==idxSNPPos]
        if (dim (ii)[1]>0) points(ii$snp_pos, -log10(ii$`p-value`), col="red", pch=16, cex=1.5)
    }
    chrName=unique(geneLoc[geneLoc$geneid==geneName,]$chr)
    genes=geneLoc[geneLoc$chr==chrName & geneLoc$s2>=start & geneLoc$s1<=end,]

    targetGene=genes[genes$geneid==geneName,]
    colnames(targetGene)[3:4]=c("start", "end")
    par(mar=c(4,5,0,2), mgp=c(2,0.5,0))
    plotGenes(genes, start, end, highlightFeatureRegion=targetGene, highlightColor="red", plotXAxis=T)
    par(mfrow=c(1,1))
    par(def.par)
}

#' Plot a GWAS region
#'
#' Plots the SNP pvalues [-log10] as well as the genes in the region.
#' @param geneName the gene to plot.  Only affects the plot title.
#' @param gwas The data frame of gwas data to draw SNPs from.  expected columns chr, pos, and a pvalue column that is by default 'p'
#' @param chrName The name of the chromosome to draw data from.
#' @param start The start position to plot from
#' @param end The end position to plot to.
#' @param geneLoc A dataframe of gene locations.  Has 4 columns: geneid, chr, s1, s2.  S1 is the start, S2 is the end of the gene.
#' @param ratioSNPsToGenes what fraction of the figure should be used by the SNP part of the plot?  [0-1].
#' @param indexSNP A SNP selected from the GWAS data as the best SNP for the region.
#' @param pvalColumn The column name to look for the SNP pvalue in.
#' @export
plotGWASRegion<-function (geneName, gwas, chrName, start, end, geneLoc, ratioSNPsToGenes=0.6, indexSNP=NULL, pvalColumn="p") {

    #if the indexSNP has a better pvalue than the GWAS, you need to set the Y lim to make it visible.

    a=gwas[gwas$CHR==chrName & gwas$POS>=start & gwas$POS<=end,]
    strTitle=paste("Gene", geneName)
    yLim=c(0, max(-log10(a[[pvalColumn]])))
    if (!is.null(indexSNP)) {
        indexSNPName=paste(indexSNP$chr, indexSNP$pos, indexSNP$a1, indexSNP$a2, sep=":")
        strTitle=paste(strTitle, "Index SNP", indexSNPName)
        maxYIndex=max(-log10(a[[pvalColumn]]))
        if (maxYIndex>yLim[2]) yLim[2]=maxYIndex
    }
    def.par <- par(no.readonly = TRUE)
    par(mar=c(0,5,2,2), mgp=c(2,0.5,0))
    mat=matrix(c(1,2),nrow=2)
    layout(mat, heights = c(ratioSNPsToGenes, (1-ratioSNPsToGenes)), respect = FALSE)
    #layout.show(n=2)

    plot (a$POS, -log10(a[[pvalColumn]]), xlab="", ylab="Observed (-logP) GWAS", main=strTitle, axes=F, cex=0.35, ylim=yLim, xlim=c(start, end))
    axis(2)
    if (!is.null(indexSNP)) {
        points(indexSNP$POS, -log10(indexSNP[[pvalColumn]]), col="red", pch=16, cex=1.5)
    }
    par(mar=c(4,5,0,2), mgp=c(2,0.5,0))
    #why do we re-define the start and end?
    #start=min(as.numeric(a$pos))
    #end=max(as.numeric(a$pos))
    genes=geneLoc[geneLoc$chr==chrName & geneLoc$s2>=start & geneLoc$s1<=end,]
    targetGene=genes[genes$geneid==geneName,]
    colnames(targetGene)[3:4]=c("start", "end")
    plotGenes(genes, start, end, highlightFeatureRegion=targetGene, highlightColor="red", plotXAxis=T)
    par(mfrow=c(1,1))
    par(def.par)
}


#' Plot genes.
#'
#' Plot a set of gene locations as a series of lines.
#' @param  genes a data frame of gene locations. Has the following columns: geneid, chr, s1, s2.
#' s1 is the start of the gene, and s2 is the end.
#' @param s the start of the region
#' @param e the end of the region
#' @param highlightFeatureRegion a data frame with two columns start and end, that defines regions to highlight a different color
#' @param highlightColor What color should the highlights be?
#' @param plotXAxis should the genomic locations with tick marks be plotted?
#' @export
plotGenes<-function (genes, s, e, highlightFeatureRegion=NULL, highlightColor="red", plotXAxis=T) {
    start=as.numeric(s)
    end=as.numeric(e)

    #short circuit for no genes.
    if (is.null(genes) || dim(genes)[1]==0) {
        plot(c(start, end), c(0,1), type='n', axes=F, ylab="", xlab="")
        text((start+end)/2, 0.5, "No Genes")
        return (NULL)
    }

    #if there are genes, make sure they are truncated to fit in the interval
    idx=which(genes$s1<start)
    if (length(idx)>0) genes[idx,]$s1=start
    idx=which(genes$s2>end)
    if (length(idx)>0) genes[idx,]$s2=end

    yRanges=seq(from= 0, to= 1, length=(dim(genes)[1]+1)*2)
    xRanges=c(start,end)
    plot(xRanges, c(0,1), type='n', axes=F, ylab="", xlab="")
    if (dim(genes)[1]==0) {
        text(mean(xRanges), 0.5, labels="No Genes Available")
        return(NULL)
    }
    par(xpd=NA)
    if (plotXAxis)  {
        pos=axTicks(1)
        #convert to 1kb.
        axis(1, at=pos, labels=pos/1e3)
        title(sub=paste("Chromosome", unique (genes$chr), "(kb)"), line=2, cex.sub=1.5)
    }
    arrowLocsX=seq(xRanges[1], xRanges[2], length=10)
    #color="blue"

    #plot genes
    labelXLoc=xRanges[1]-((xRanges[2]-xRanges[1])*0.1)
    labelXLocS=xRanges[1]-((xRanges[2]-xRanges[1])*0.02)
    labelXLocE=xRanges[2]+((xRanges[2]-xRanges[1])*0.02)
    for (i in 1:dim(genes)[1]) {
        gene= genes[i,]
        color=getColorForFeature(gene$s1, gene$s2,
                                 highlightFeatureRegion=highlightFeatureRegion, highlightColor=highlightColor)


        yStart= yRanges [(2*i)]
        yEnd= yRanges[((2*i)+1)]
        yMid=(yStart+yEnd)/2
        #only need a few columns from exons - start, end, strand.
        text(labelXLoc, ((yStart+yEnd)/2), gene$geneid, cex=0.75, col=color)
        lines(x=c(gene$s1,gene$s2), y=c(yMid,yMid), col=color)

    }

}


#check the start/end coordinates for what color the feature should be.
getColorForFeature<-function (start, end, defaultColor="black", highlightFeatureRegion=NULL, highlightColor="red") {
    if (is.null(highlightFeatureRegion)) return (defaultColor)
    #if (highlightFeatureRegion$start<end & highlightFeatureRegion$end>start) return (highlightColor)
    if (highlightFeatureRegion$start==start & highlightFeatureRegion$end==end) return (highlightColor)
    return (defaultColor)
}
