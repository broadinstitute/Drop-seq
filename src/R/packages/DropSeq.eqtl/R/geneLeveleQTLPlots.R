# MIT License
#
# Copyright 2025 Broad Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#
#' Plot many eQTL expression vs genotype plots on one page.
#'
#' @param snps_location_file_name The location of each SNP.  Same number of lines as the SNP_file_name, encodes the snp name, chromosome, position.
#' @param expression_file_name A matrix of expression data, 1 row per gene, 1 column per sample.
#' @param SNP_file_name A matrix of genotypes, 1 row per SNP, 1 column per sample.  Encoded as 0,1,2 copies of the alternate allele.
#' @param eQTLPermutationResultFile The permuted eQTL result file
#' @param outPDF The file location to output the PDF of SNP/Expression plots.
#' @param outIndexSNPs The file location of the gene/SNPs with a q-value <= qValueThreshold
#' @param qValueThreshold Genes with a score <= this threshold will be emitted in the outIndexSNPs and outPDF.
#' @param geneList Override the qValueThreshold and plot these genes instead.
#' @param gene_location_file_name Obsolete, will be removed in the future.
#' @suggest DropSeq.utilities
#' @export
plotGeneQTLs<-function (snps_location_file_name, expression_file_name, SNP_file_name, eQTLPermutationResultFile,
                        outPDF=NULL, outIndexSNPs=NULL, qValueThreshold=0.05, geneList=NULL,
                        gene_location_file_name=NULL) {
    cisEQTLs=DropSeq.utilities::fastRead(eQTLPermutationResultFile)
    if ("qvalue" %in% colnames(cisEQTLs)) {
      cisEQTLs=cisEQTLs[cisEQTLs$qvalue<=qValueThreshold,]
      if (nrow(cisEQTLs) == 0) {
        stop("No significant cisEQTLs found; consider changing qValueThreshold")
      }
      cisEQTLs=cisEQTLs[order(cisEQTLs$qvalue, decreasing=F),]
    }
    snpLoc=DropSeq.utilities::fastRead(snps_location_file_name)
    expressionMatrix=DropSeq.utilities::fastRead(expression_file_name)
    geneotypes=DropSeq.utilities::fastRead(SNP_file_name)
    genesToPlot = unique(cisEQTLs$gene)
    if (!is.null(geneList)) genesToPlot=geneList

    if (length(genesToPlot) == 0) {
      warning("No genes to plot! No outPDF will be generated")
    } else {
      pdf(outPDF, width=16, height=9)
      par(mfrow=c(2,4))
      for (i in 1:length(genesToPlot)) {
          geneName=genesToPlot[i]
          cat (geneName, "\n")
          #indexSNP=plotRegion(geneName, cisEQTLs, snpLoc, ratioSNPsToGenes=0.6, plot=F)
          plotGenotypesVsExpression (geneName, cisEQTLs, snpLoc, expressionMatrix, geneotypes)
      }
      dev.off()
    }

    #write out index SNPs report.
    if (!is.null(outIndexSNPs)) {
        write.table(cisEQTLs, outIndexSNPs, row.names = F, col.names = T, quote=F, sep="\t")
    }

}






#' Plot the genotype by expression regression analysis plot for a single SNP and gene.
#'
#' @param geneName The name of the gene to plot
#' @param cisEQTLs The list of eQTL results (can be permuted or all hits.)
#' @param snpLoc A data frame of snp locations
#' @param expressionMatrix A data frame of expression data per gene/donor
#' @param geneotypes A matrix of genotypes (copies of alternate allele) per snp/donor
#' @param strTitlePrefix Put this at the front of the plot
#' @param indexSNP A SNP name to plot.  If not supplied, and you passed in the permuted hits, uses the best hit SNP for the gene.
#' @param ylim Limit the expression of the Y axis (optional)
#' @param eQTLResult All SNP by gene results for this experiment (optional).  If provided, get the empiric pvalue from this dataframe
#' if the SNP/Gene pair does not exist in the cisEQTLs result.  Useful when cisEQTLs is the permuted result and an index SNP is defined.
plotGenotypesVsExpression<-function (geneName, cisEQTLs, snpLoc, expressionMatrix, geneotypes, strTitlePrefix="", indexSNP=NULL, ylim=NULL, eQTLResult=NULL) {
    z=cisEQTLs[cisEQTLs$gene==geneName,]
    #select the index SNP if not provided as the best hit for this gene.
    idx=match(z$SNP, snpLoc$snp)
    z=cbind(z, snpLoc[idx,c("chr", "pos")])
    if (is.null(indexSNP)) indexSNP=z$SNP
    if (!is.null(indexSNP)) z=z[match(indexSNP, z$SNP),]
    #add in the empiric p-value if it's passed in and it's not in the cis data.
    empiric_pvalue=z$`p-value`
    if (!is.null(eQTLResult) & is.na(empiric_pvalue)) {
        x=eQTLResult[eQTLResult$gene==geneName & eQTLResult$SNP==indexSNP,]
        empiric_pvalue=x$'p-value'
    }
    e=as.numeric(expressionMatrix[expressionMatrix$id==geneName,-1,with=F])
    g=as.numeric(geneotypes[geneotypes$id==indexSNP,-1])
    refAllele=strsplit (indexSNP, ":")[[1]][3]
    altAllele=strsplit (indexSNP, ":")[[1]][4]
    #if the SNP is missing from the data set, don't run the regression, and generate an empty plot.
    r=NULL
    if (any(!is.na(g)) & any(!is.na(e)))
      r=lm(e ~ g)
    #r2=round (summary(r)$adj.r.squared,3)

    #set up the complicated title string.
    strTitle=paste(strTitlePrefix, geneName, sep=" ")
    if (!is.null(indexSNP)) strTitle=paste(strTitle, indexSNP, sep="\n")

    #I think z always has 1 row because R.
    if (dim (z)[1]>0) {
        if (is.null(z$qvalue) || is.na(z$qvalue)) {
            strTitle=paste(strTitle, "\np-value", format(empiric_pvalue ,digits=3, scientific=T))
        } else {
            strTitle=paste(strTitle, "\nq-value", format(z$qvalue,digits=3, scientific=T))
        }
    }

    par(mar=c(5,5,6,1))
    if (is.null(ylim)) ylim=c(min(0, e), max(0, e))
    plot (e ~ g, col="black", main=strTitle, ylab="EXPRESSION", xlab="GENOTYPE", xlim=c(-0.2,2.2), axes=F, cex.lab=1.5, cex.main=2, ylim=ylim)
    #axis(1, at=c(0,1,2), labels=c("RR", "RA", "AA"), cex.lab=2, cex.axis=1.5)
    axis(1, at=c(0,1,2), labels=c(paste(refAllele,refAllele,sep=""), paste(refAllele,altAllele,sep=""), paste(altAllele,altAllele,sep="")), cex.lab=2, cex.axis=1.5)
    axis(2, cex.lab=2, cex.axis=1.5)
    if (!is.null(r)) {
      yStart=r$coefficients[[1]]
      yEnd=2*r$coefficients[[2]]+yStart
      lines (x=c(0,2), y=c(yStart, yEnd), lty=2, col="black")
    }
    par(mar=c(5, 4, 4, 2) + 0.1)
}
