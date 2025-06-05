#inFiles=c("/downloads/eQTL/rNPC_pilot/rNPC_pilot.eQTL_permuted_results.txt", "/downloads/eQTL/rNPC_Alive/rNPC_Alive_AllPools.auto.eQTL_permuted_results.txt")
#experimentNames=c("rNPC Pilot", "rNPC Alive")

#inFiles=c("/downloads/eQTL/heat1-6/heat1-6.eQTL_permuted_results.txt", "/downloads/eQTL/villageC/out_nextseq_novaseq/VillageC.maf_0.2.eQTL_permuted_results.txt")
#experimentNames=c("heat 1-6", "Village C")

#inFiles=c("/downloads/eQTL/rNPC_pilot/rNPC_pilot.eQTL_permuted_results.txt", "/downloads/eQTL/rNPC_Alive/rNPC_Alive_AllPools.auto.eQTL_permuted_results.txt", "/downloads/eQTL/heat1-6/heat1-6.eQTL_permuted_results.txt", "/downloads/eQTL/villageC/out_nextseq_novaseq/VillageC.maf_0.2.eQTL_permuted_results.txt")
#experimentNames=c("rNPC Pilot", "rNPC Alive", "heat 1-6", "Village C")
#library (data.table);library (VennDiagram);library (RColorBrewer); library (vegan)

#' Plot Venn Diagrams of genes found in multiple experiments
#'
#' @param inFiles A list of eQTL result files that have a qvalue column.
#' @param experimentNames The names of the experiments for the plot
#' @param qvalueThreshold the maximum q-value genes that will be retained.
#' @import RColorBrewer VennDiagram data.table futile.logger grid
#' @export
plotGeneIntersection<-function (inFiles, experimentNames, qvalueThreshold=0.05) {
    cols=brewer.pal(3, "Set2")[1:length(experimentNames)]

    d=lapply(inFiles, fread)
    d=lapply(d, function (x) x[x$qvalue<=qvalueThreshold,])
    genes=lapply (d, function (x) unique (x$gene))
    par(mar=c(5,6,5,6))
    #fix so venn diagram doesn't output random log message files.
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
    v=venn.diagram(genes, filename=NULL,
                   fill=cols, alpha=rep(0.5, length(d)), cat.cex=2, cex = 2, main.cex=2.5, cat.fontface=4, category.names=experimentNames, main="eQTL overlap", cat.pos = c(180,180))
    grid.newpage();grid.draw(v)
}

#metaCellFile="/downloads/meta_cell/d21NGN2Live.meta_cell.expression.txt"

#' Plot diversity of UMIs/Donors in an eQTL experiment.
#' @param metaCellFile A single file of metacells.
#' @param outPDF The PDF file to write to, or NULL.
#' @import vegan
#' @export
plotDiversity <-function (metaCellFile, outPDF=NULL) {
    a=read.table(metaCellFile, header=T, stringsAsFactors = F, check.names=F)
    d=colSums(a[,-1])
    d=d[order(d, decreasing = T)]
    div=vegan::diversity(d)
    expName=sub (".meta_cell.expression.txt", "", basename (metaCellFile))
    strTitle=paste(expName, "\nDistribution of UMIs across donors\n",length(d), " donors; ", round (sum (d)/1e6,1), "M UMIs; SW Div: ", sprintf("%.2f",round(div,2)), "; SW Eq: ", sprintf("%.2f", round (div/log(length(d)),2)), sep="")
    if (!is.null(outPDF)) pdf(outPDF, width=8, height=6)
    par(mar=c(8,4,4,2))
    barplot ((d/1e6), las=3, cex.axis=0.5, cex.names=0.75, cex.lab=1, ylab="Total UMIs [millions]", names.arg = names(d), main=strTitle, col="light blue")
    par(mar=c(4,4,4,4))
    if (!is.null(outPDF)) dev.off()
}

#metaCellFiles=c("/downloads/meta_cell/NGN2d06Fixed.meta_cell.expression.txt", "/downloads/meta_cell/NGN2d14Live.meta_cell.expression.txt", "/downloads/meta_cell/d21NGN2Live.meta_cell.expression.txt", "/downloads/meta_cell/NGN2d28Live.meta_cell.expression.txt","/downloads/meta_cell/NGN2d35Live.meta_cell.expression.txt")
#experimentNames=c("d6 fixed", "d14", "d21", "d28", "d35");plotTitle="NGN2"

#metaCellFiles=c("/downloads/meta_cell/MNd06Fixed.meta_cell.expression.txt", "/downloads/meta_cell/d14MNLive.meta_cell.expression.txt", "/downloads/meta_cell/d28MNLive.meta_cell.expression.txt")
#experimentNames=c("d6 fixed", "d14", "d28");plotTitle="Motor Neurons"

#' A stacked area plot by donor assignable UMIs.
#' @param metaCellFiles A vector of metaCell files
#' @param experimentNames a vector of matching label names
#' @param plotTitle A plot title.
#' @import areaplot
#' @export
plotDonorBalanceByAssignableUMIs<-function (metaCellFiles, experimentNames, plotTitle="") {
    getFractionOfDonor<-function (x) {
        file.exists(x)
        a=read.table(x, header=T, stringsAsFactors = F, check.names=F)
        d=colSums(a[,-1])
        d=d/sum(d)
        return (d)
    }

    r=lapply(metaCellFiles, getFractionOfDonor)
    donorNames=unique (unlist((lapply(r, names))))

    #order all of the outputs by donor name.
    f2<-function (x, donorNames) {
        z=x[donorNames]
        z[is.na(z)]<-0
        names (z)=donorNames
        return (z)
    }

    rr=lapply(r, f2, donorNames)
    rrr=do.call(cbind, rr)
    colnames(rrr)=experimentNames

    #cols=rainbow(dim (rrr))
    #cols=RColorBrewer::display.brewer.all()
    cols=brewer.pal(9, "Set3")

    plotOne<-function (x, plotTitle="") {
        suppressWarnings(areaplot(t(x), col=cols, xlab="Passage", axes=F))
        axis(1, at=1: dim (x)[2], labels=colnames(x))
        title (paste(plotTitle, "\n", "Donor Representation by assignable UMIs"))
    }
    plotOne(rrr, plotTitle=plotTitle)

    #pairwise correlation?
    getCor<-function (index, x) {
        z=cor (x[,index], x[,(index+1)])
        label=paste(colnames (x)[index], "-", colnames (x)[index+1])
        names (z)=label
        return (z)
    }

    corelations=sapply(1:(dim(rrr)[2]-1), getCor, rrr)
    barplot (corelations, ylim=c(0,1), ylab="correlation across passages", main=plotTitle, col="light blue")

    #rrr=rrr[,-1]
    #plotOne(rrr, plotTitle=plotTitle)
    #write.table(rrr, "/downloads/text.txt", row.names=T, col.names=T, quote=F, sep="\t")
}

#' Scatter plot of p-values for (SNP, gene) that overlap between 2 sets of eQTL_permuted_results
#' @param method1Path tab-separated file with columns (SNP, gene, p.value) (and maybe other columns).  Plotted on x-axis.
#' @param method2Path tab-separated file with columns (SNP, gene, p.value) (and maybe other columns).  Plotted on y-axis.
#' @param method1Label x-axis label.
#' @param method2Label y-axis label.
#' @param title plot title
#' @param outFile If specified, create PDF with plot
#' @export
scatterPlotPvaluesForTwoMethods<-function(method1Path, method2Path,
    method1Label=basename(method1Path), method2Label=basename(method2Path),
    title=NA, outFile=NULL) {
    # read.table converts column p-value into p.value
    desired_columns = c("SNP", "gene", "gene_permuted_pvalue")
    method1 = read.table(method1Path, sep="\t", header=TRUE, stringsAsFactors = FALSE)[, desired_columns]
    method2 = read.table(method2Path, sep="\t", header=TRUE, stringsAsFactors = FALSE)[, desired_columns]
    colnames(method1)[which(colnames(method1) == 'gene_permuted_pvalue')] = 'p.value1'
    colnames(method2)[which(colnames(method2) == 'gene_permuted_pvalue')] = 'p.value2'
    merged = merge(method1, method2)

    if (!is.null(outFile)) {
        pdf(outFile)
    }

    plot(-log10(merged$p.value1), -log10(merged$p.value2), xlab=method1Label, ylab=method2Label, main=title)

    if (!is.null(outFile)) {
        dev.off()
    }
}