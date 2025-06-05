# MIT License
#
# Copyright 2017 Broad Institute
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

# source ("/Users/nemesh/dropseqrna/transcriptome/R/packages/DropSeq.eqtl/R/Matrix_eQTL_CommonFunctions.R")
# library (data.table)
# library (IRanges)

## Settings

# data.dir="/downloads/d21STORMI"
# eQTL.dir="/downloads/d21STORMI"
# base.name="d21Ngn2.maf_0.05_cisDist_250kb"
# qValueThreshold=0.05
# #feed in the non permuted results to plot arbitrary SNPs.
# eQTLPermutationResultFile="/downloads/d21STORMI/d21Ngn2.maf_0.05_cisDist_250kb.eQTL_results.txt"

#plotGeneQTLsFromBaseFileNames(data.dir="/downloads/d06NGN2/data", base.name="NGN2d06Fixed_MAF_0.2", eQTL.dir="/downloads/d06NGN2/eQTL")
#' Plot expression by genotype plots for eQTLs below some qvalue threshold.
#'
#' @param data.dir The directory data resides in.
#' @param eQTL.dir The directory the eQTL permuted results reside in.  Files will be written out here.
#' @param base.name The prefix of all files for this data set.
#' @param eQTLPermutationResultFile a set of results for the eQTLs that have been permuted and have qvalues.  If not supplied,
#' will attempt to infer the proper file name given the base name and eQTL.dir where it usually resides.
#' @param qValueThreshold the threshold at which eQTLs are selected.  Selects all eQTLs with <=threhold value.
#' @return The results of adaptive permutation (ALL SNPs per gene.)
#'
#' @export
plotGeneQTLsFromBaseFileNames<-function (data.dir, base.name, eQTL.dir, eQTLPermutationResultFile=NULL,qValueThreshold=0.05) {

    # Genotype file name
    snps_location_file_name = getSNPLocationFileName(data.dir, base.name)
    SNP_file_name= getSNPFileName(data.dir, base.name)

    # Gene expression file name
    expression_file_name = getExpressionFileName(data.dir, base.name)
    gene_location_file_name = getGeneLocationFileName(data.dir, base.name)

    if (is.null(eQTLPermutationResultFile))
      eQTLPermutationResultFile=paste(paste(eQTL.dir, "/", base.name, ".eQTL_eigenMT_results.txt", sep=""))

    # Output file names
    outPDF=paste(eQTL.dir, "/", sub (".txt", "", basename (eQTLPermutationResultFile)), ".pdf", sep="")
    outIndexSNPs=paste(eQTL.dir, "/", sub (".txt", "", basename (eQTLPermutationResultFile)), ".index_snps.txt", sep="")

    plotGeneQTLs(gene_location_file_name, snps_location_file_name, expression_file_name, SNP_file_name, eQTLPermutationResultFile,
                 outPDF, outIndexSNPs, qValueThreshold)

}


#
#' Plot many eQTL expression vs genotype plots on one page.
#'
#' @param gene_location_file_name The location of each gene.  Same number of lines as the expression_file_name, encodes the gene name, chromosome, start, end.
#' @param snps_location_file_name The location of each SNP.  Same number of lines as the SNP_file_name, encodes the snp name, chromosome, position.
#' @param expression_file_name A matrix of expression data, 1 row per gene, 1 column per sample.
#' @param SNP_file_name A matrix of genotypes, 1 row per SNP, 1 column per sample.  Encoded as 0,1,2 copies of the alternate allele.
#' @param eQTLPermutationResultFile The permuted eQTL result file
#' @param outPDF The file location to output the PDF of SNP/Expression plots.
#' @param outIndexSNPs The file location of the gene/SNPs with a q-value <= qValueThreshold
#' @param qValueThreshold Genes with a score <= this threshold will be emitted in the outIndexSNPs and outPDF.
#' @param geneList Override the qValueThreshold and plot these genes instead.
#' @export
plotGeneQTLs<-function (gene_location_file_name, snps_location_file_name, expression_file_name, SNP_file_name, eQTLPermutationResultFile, outPDF=NULL, outIndexSNPs=NULL, qValueThreshold=0.05, geneList=NULL) {
    cisEQTLs=fread(eQTLPermutationResultFile)
    if ("qvalue" %in% colnames(cisEQTLs)) {
      cisEQTLs=cisEQTLs[cisEQTLs$qvalue<=qValueThreshold,]
      if (nrow(cisEQTLs) == 0) {
        stop("No significant cisEQTLs found; consider changing qValueThreshold")
      }
      cisEQTLs=cisEQTLs[order(cisEQTLs$qvalue, decreasing=F),]
    }
    geneLoc=fread(gene_location_file_name)
    snpLoc=fread(snps_location_file_name)
    expressionMatrix=fread(expression_file_name)
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


# dataDirList=c("/broad/mccarroll/dropulation/analysis/eQTL/Levy_Perturbation/Levy_village_iPSC.permissive_44_donor_WGS/maf_0.20_cisDist_10kb",
#               "/broad/mccarroll/dropulation/analysis/eQTL/Levy_Perturbation/Levy_village_NPC.permissive_44_donor_WGS/maf_0.20_cisDist_10kb",
#               "/broad/mccarroll/dropulation/analysis/eQTL/Levy_Perturbation/Levy_village_Neuron.permissive_44_donor_WGS/maf_0.20_cisDist_10kb")
#
# basenames=c("Levy_village_iPSC.permissive_44_donor_WGS.maf_0.20_cisDist_10kb", "Levy_village_NPC.permissive_44_donor_WGS.maf_0.20_cisDist_10kb", "Levy_village_Neuron.permissive_44_donor_WGS.maf_0.20_cisDist_10kb")
# expNames=sub("Levy_village_", "", sub(".permissive_44_donor_WGS.maf_0.20_cisDist_10kb", "", basenames))
# outPDF="/downloads/Levy_gwas_eQTLs.pdf"
# mode="union"; qValueThreshold=0.05;
# geneList=c("EMB", "HYI","KANSL1", "LINC01068", "LSM1", "MOB4", "MSI2", "NEGR1", "NRXN1", "PDIA3", "PSMA4", "SETD6", "THOC7")
# geneList=c("EMB", "HYI","KANSL1", "LINC01068", "LSM1", "MOB4", "MSI2", "NEGR1", "NRXN1", "PDIA3", "PSMA4", "THOC7")

#' For multiple experiments, plot the eQTL results for the same gene same SNP across multiple experiments.
#'
#' Each input parameter is a list of files. Results are plotted for genes that are significant (q-value < threshold) in any experiment.
#' @param dataDirList The directory data resides in - this is a vector of entries.
#' @param basenames The prefix of all files for this data set - this is a vector of entries of the same length as dataDirList
#' @param expNames The experiment names for the data, which will be in the header of the plots.
#' @param outPDF Where to write a PDF output
#' @param qValueThreshold Genes with a score <= this threshold will be emitted in the outIndexSNPs and outPDF.
#' @param geneList Override the qValueThreshold and plot these genes instead.
#' @param mode One of union/intersect/disagree.  Union picks all genes in any experiment.  Intersect picks genes in all experiments.  Disagree selects genes that are positive in one experiment and negative in another.
#' @export
plotGeneQTLsForMultipleExperimentsBaseNames<-function (dataDirList, basenames, expNames, outPDF=NULL, qValueThreshold=0.05, geneList=NULL, mode=c("union", "intersect", "disagree")) {
    applyFunc<-function (index, dataDirList, basenames, FUNC) {
        r=FUNC(dataDirList[index], basenames[index])
        if (any(!file.exists(r)))
            stop("One more more files missing ", paste(r, collapse=","))
        return (r)
    }
    # Genotype file name
    snps_location_file_name = sapply(1:length(dataDirList), applyFunc, dataDirList, basenames, FUNC=getSNPLocationFileName)
    SNP_file_name= sapply(1:length(dataDirList), applyFunc, dataDirList, basenames, FUNC=getSNPFileName)

    # Gene expression file name
    expression_file_name = sapply(1:length(dataDirList), applyFunc, dataDirList, basenames, FUNC=getExpressionFileName)
    gene_location_file_name = sapply(1:length(dataDirList), applyFunc, dataDirList, basenames, FUNC=getGeneLocationFileName)

    #snp x gene results files
    eqtl_result_file_name=sapply(1:length(dataDirList), applyFunc, dataDirList, basenames, FUNC=getEQTLResults)

    #permutation files
    eQTLPermutationResultFile = sapply(1:length(dataDirList), applyFunc, dataDirList, basenames, FUNC=getIndexSNPsFileName)


    plotGeneQTLsForMultipleExperiments(gene_location_file_name, snps_location_file_name, expression_file_name, SNP_file_name,
                                       eqtl_result_file_name, eQTLPermutationResultFile, expNames=expNames, outPDF=outPDF,
                                       qValueThreshold=qValueThreshold, geneList=geneList, mode=mode)

}

#' For multiple experiments, plot the eQTL results for the same gene same SNP across multiple experiments.
#'
#' Each input parameter is a list of files. Results are plotted for genes that are significant (q-value < threshold) in any experiment.
#' @param expNames the experiment names for the input files
#' @param mode One of union/intersect/disagree.  Union picks all genes in any experiment.  Intersect picks genes in all experiments.  Disagree selects genes that are positive in one experiment and negative in another.
#' @inheritParams plotGeneQTLs
#' @export
plotGeneQTLsForMultipleExperiments<-function (gene_location_file_name, snps_location_file_name, expression_file_name, SNP_file_name, eqtl_result_file_name, eQTLPermutationResultFile, expNames, outPDF=NULL, qValueThreshold=0.05, geneList=NULL, mode=c("union", "intersect", "disagree")) {
    n=length(gene_location_file_name)
    cisEQTLs=lapply(eQTLPermutationResultFile, fread)
    eqtl_result=lapply(eqtl_result_file_name, fread)
    geneLoc=lapply(gene_location_file_name, fread)
    snpLoc=lapply(snps_location_file_name, fread)
    expressionMatrix=lapply(expression_file_name, fread)
    geneotypes=lapply(SNP_file_name, fread)
    genesToPlot = getGenesToPlot(eQTLPermutationResultFile, qValueThreshold, mode)
    if (!is.null(geneList)) genesToPlot=geneList

    indexSNPs=getIndexSNPs(eQTLPermutationResultFile, genesToPlot)

    pdf(outPDF, width=16, height=9)
    if (n<4) par(mfrow=c(1,n))
    if (n==3) par(mfrow=c(2,n))
    #if (n==4) par(mfrow=c(2,2))
    if (n>=4) par(mfrow=c(2,(n/2)))

    #quick function to get the maximum expression of a gene across experiments.
    getExpression <- function(geneName, expressionMatrix, min_or_max) {
        gm<-function (e, geneName) {
          ee=e[e$id==geneName,-1]
          if (dim(ee)[1]==0) return (0)
            return(min_or_max(ee))
        }
        min_or_max(unlist(lapply(expressionMatrix, function(e, geneName) gm(e, geneName), geneName)))
    }

    getMaxExpression <- function(geneName, expressionMatrix) {
        getExpression(geneName, expressionMatrix, max)
    }

    getMinExpression <- function(geneName, expressionMatrix) {
        getExpression(geneName, expressionMatrix, min)
    }

    for (i in 1:length(genesToPlot)) {
        geneName=genesToPlot[i]
        maxExpression <- max(0, round(ceiling(getMaxExpression(geneName, expressionMatrix)) * 1.2))
        minExpression <- min(0, round(floor(getMinExpression(geneName, expressionMatrix)) * 1.2))
        cat (geneName, "\n")
        for (j in 1:n)
            plotGenotypesVsExpression (geneName, cisEQTLs[[j]], snpLoc[[j]], expressionMatrix[[j]], geneotypes[[j]], expNames[j], indexSNP=indexSNPs[geneName], ylim=c(minExpression, maxExpression), eqtl_result[[j]])
    }
    dev.off()
}

#get the index SNP (the best result out of many eQTL results) to plot for each listed gene.
getIndexSNPs<-function (eQTLPermutationResultFileList, genesToPlot) {
    r=lapply(eQTLPermutationResultFileList, fread)
    getBestIndexSNP<-function (geneName, r) {
        bestPerExp=lapply (r, function (x) x[x$gene==geneName,])
        idxBest=which.min(lapply (bestPerExp, function (x) x$qvalue))
        bestPerExp[[idxBest]]$SNP
    }
    indexSNPs=sapply(genesToPlot,getBestIndexSNP, r)
    return (indexSNPs)
}
#eQTLPermutationResultFileList=eQTLPermutationResultFile
getGenesToPlot<-function (eQTLPermutationResultFileList, qValueThreshold=0.05, mode=c("union", "intersect", "disagree")) {
    mode=match.arg(mode, choices=c("union", "intersect", "disagree"))
    r=lapply(eQTLPermutationResultFileList, fread)
    rPass=lapply (r, function (x) x[x$qvalue<=qValueThreshold,])
    rFail=lapply (r, function (x) x[x$qvalue>qValueThreshold,])
    genesPass=lapply (rPass, function (x) x$gene)
    genesFail=lapply (rFail, function (x) x$gene)
    if (mode=="union") result=Reduce(union, genesPass)
    if (mode=="intersect") result=Reduce (intersect, genesPass)
    if (mode=="disagree") {
        #for disagree, we want genes that fail in any experiment and pass in any experiment.
        #if a gene passes in all experiments, or fails in all, it's not interesting.
        genesPassAny=Reduce(union, genesPass)
        genesFailAny=Reduce(union, genesFail)
        result=intersect(genesPassAny, genesFailAny)
    }
    return (result)
}

#plotGenotypesVsExpression(geneName="IFITM3", cisEQTLs, snpLoc, expressionMatrix, geneotypes, strTitlePrefix="11:320394:C:T", indexSNP="11:320394:C:T", ylim=NULL)
#plotGenotypesVsExpression(geneName="IFITM3", cisEQTLs, snpLoc, expressionMatrix, geneotypes, strTitlePrefix="11:320836:C:T", indexSNP="11:320836:C:T", ylim=NULL)


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
#' @export
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


plotGenotypesVsExpressionSimple<-function (geneName, indexSNP, snpLoc, expressionMatrix, geneotypes, strTitlePrefix="", qvalue=NULL, ylim=NULL) {

    e=as.numeric(expressionMatrix[expressionMatrix$id==geneName,-1,with=F])
    g=as.numeric(geneotypes[geneotypes$id==indexSNP,-1])
    refAllele=strsplit (indexSNP, ":")[[1]][3]
    altAllele=strsplit (indexSNP, ":")[[1]][4]
    r=lm(e ~ g)
    s=summary(r)

    r2=round (s$adj.r.squared,3)
    pval=s$coefficients[2,4]

    if (is.null(qvalue)) {
        strTitle=paste(strTitlePrefix, geneName, "\np-value", format(pval, digits=3, scientific=T))
    } else {
        strTitle=paste(strTitlePrefix, geneName, "\nq-value", format(qvalue,digits=3, scientific=T))
    }

    par(mar=c(5,5,4,1))
    if (is.null(ylim)) ylim=c(min(0, e), max(0, e))
    plot (e ~ g, col="black", main=strTitle, ylab="EXPRESSION", xlab="GENOTYPE", xlim=c(-0.2,2.2), axes=F, cex.lab=1.5, cex.main=2, ylim=ylim)
    #axis(1, at=c(0,1,2), labels=c("RR", "RA", "AA"), cex.lab=2, cex.axis=1.5)
    axis(1, at=c(0,1,2), labels=c(paste(refAllele,refAllele,sep=""), paste(refAllele,altAllele,sep=""), paste(altAllele,altAllele,sep="")), cex.lab=2, cex.axis=1.5)
    axis(2, cex.lab=2, cex.axis=1.5)
    yStart=r$coefficients[[1]]
    yEnd=2*r$coefficients[[2]]+yStart
    #abline (r, lty=2)
    lines (x=c(0,2), y=c(yStart, yEnd), lty=2, col="black")
}



# #geneName="CRYZ"
# plotRegion<-function (geneName, cisEQTLs, snpLoc, ratioSNPsToGenes=0.6, plot=T) {
#
#   if (plot) {
#     #"FDR", format(z[idx,]$FDR, digits=3)
#     pdf ("/downloads/CHRNA5_eQTL.pdf", width=8.7, height=6)
#     strTitle=paste("Gene", geneName, "Index SNP", indexSNP)
#     def.par <- par(no.readonly = TRUE)
#     par(mar=c(0,3,2,2), mgp=c(2,0.5,0))
#     mat=matrix(c(1,2),nrow=2)
#     layout(mat, heights = c(ratioSNPsToGenes, (1-ratioSNPsToGenes)), respect = FALSE)
#     #layout.show(n=2)
#     #plot (z$pos, -log10(z$`p-value`), xlab="position", ylab="Observed (-logP)", main=strTitle, axes=F)
#     plot (z$pos, -log10(z$`p-value`), xlab="", ylab="Observed (-logP)", main=strTitle, axes=F, cex=0.35)
#     axis(2)
#     points(z[z$SNP==indexSNP]$pos, -log10(z[z$SNP==indexSNP,]$`p-value`), col="red", pch=16, cex=1.5)
#
#     chrName=as.numeric(unique(geneLoc[geneLoc$geneid==geneName,]$chr))
#     genes=geneLoc[geneLoc$chr==chrName & geneLoc$s2>=min(z$pos) & geneLoc$s1<=max(z$pos),]
#     par(mar=c(4,7,0,2), mgp=c(2,0.5,0))
#     plotGenes(genes, min(z$pos), max(z$pos), color="blue", highlightFeatureRegion=NULL, highlightColor="red", plotXAxis=T)
#     par(mfrow=c(1,1))
#     par(def.par)
#     dev.off()
#   }
#
#   return (indexSNP)
# }
