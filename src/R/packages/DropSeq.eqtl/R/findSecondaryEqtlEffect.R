#library (MatrixEQTL)
# options(warn=2, error=recover)

#  10000

# dropSeqToolsLocation="/Users/nemesh/jn_branch"
# pythonPath="/Users/nemesh/opt/anaconda3/bin/python"
#
# rootDir="/broad/mccarroll/dropulation/analysis/eQTL/Levy_Perturbation/Levy_village_Neuron_INFa_7hr.permissive_44_donor_WGS/maf_0.20_cisDist_10kb"
# rootName="Levy_village_Neuron_INFa_7hr.permissive_44_donor_WGS.maf_0.20_cisDist_10kb"
#
# genotypeFile=paste(rootDir, "/", rootName, ".genotype_matrix.txt.gz", sep="")
# geneLocationFile=paste(rootDir, "/", rootName, ".gene_locations.txt", sep="")
# variantLocationFile=paste(rootDir, "/", rootName, ".variant_locations.txt", sep="")
# expressionFile=paste(rootDir, "/", rootName, ".gene_expression_peer.txt", sep="")
# indexSNPFile=paste(rootDir, "/", rootName, ".eQTL_eigenMT_results.txt", sep="")
# snpByGeneResultFile=paste(rootDir, "/", rootName, ".eQTL_results.txt.gz", sep="")
# outFile=paste(rootDir, "/", rootName, ".secondary_scan.txt", sep="")
# outPDF=paste(rootDir, "/", rootName, ".secondary_scan.1.pdf", sep="")
# outFileEigenMT=paste(rootDir, "/", rootName, ".secondary_scan.eigenMT.txt", sep="")
#
# cisDist=10000; fdrThreshold=0.05; qvalueThreshold=0.05

# findSecondaryEqtlEffects(indexSNPFile, expressionFile, variantLocationFile, genotypeFile, outFile, outFileEigenMT, outPDF, fdrThreshold=0.1, qvalueThreshold=0.05, cisDist=250000, dropSeqToolsLocation="/Users/nemesh/jn_branch", pythonPath="/Users/nemesh/opt/anaconda3/bin/python")


#' Search for secondary eQTL effects at eGenes
#'
#' Given a set of eGenes, test to see if there are additional eQTL effects from SNPs that
#' are not highly correlated with the index SNP.
#'
#' This runs a linear regression via MatrixEQTL at each eGENE against every SNP in the window
#' except for the original index SNP.  The formula for each SNP in the window is expression ~ SNP + INDEX_SNP
#'
#' Given those results, an estimate of the FDR is calculated for secondary eQTLs based on the number of tests generated
#' by the original eigenMT analysis of the window.  New SNPs with an FDR are included if they are below fdrThreshold.  EigenMT
#' is then re-run to include both primary and secondary effects.
#'
#' If the outPDF argument is not null, 2 pages of plots per eGene will be generated.  The first page has the
#' region plot for the index SNP on the top, and the region plot after the index SNP has been accounted for on the bottom.
#'
#'  The second page contains plots of the effect of the primary and secondary effects on expression of the gene, their
#'  individual and combined betas and how correlated the genotypes of the two SNPs are (R2).
#' @param indexSNPFile A file containing the index SNP for each gene.  This contains the columns: SNP, gene, gene_start_pos, gene_end_pos, TESTS, qvalue.
#' @param snpByGeneResultFile A file containing all SNP by gene interactions for the original eQTL scan.
#' @param geneLocationFile The locations of genes in the scan.  This must have 4 columns in the header: geneid, chr, s1, s2.
#' @param expressionFile The expression values of the genes.  A matrix with an id column, followed by 1 column per donor, 1 row per gene.
#' @param variantLocationFile The location of the SNPs in the scan.  This must have 3 columns in the header: snp, chr, pos
#' @param genotypeFile The matrix of genotypes for the SNPs in the scan.  A matrix with an id column, followed by 1 column per donor, 1 row per SNP.
#' Genotypes are encoded as the number of copies of the alternate allele.
#' @param outFile Location to write the outfile to.  This contains all the new SNP x Gene interactions.
#' @param outFileEigenMT Location to write the final eigenMT result to.  This contains both the original index SNP results and new secondary effects
#' after correction by eigenMT.  There may be more than one SNP listed per gene, on seperate lines.
#' @param outPDF An optional file to write plots to.
#' @param fdrThreshold A FDR threshold to include secondary effects into the eigenMT testing.  Setting this to a high value will
#' result in all genes having a secondary effect included, which will increase multiple hypothesis testing burden.  Suggested value of 0.05.
#' Not all SNPs with an FDR threshold estimated <= fdrThreshold will have a final permuted p-value or q-value less than the threshold.
#' @param qvalueThreshold This threshold selects the set of eGenes to scan for secondary effects.  Suggested value 0.05.
#' @param cisDist What window size should be searched for secondary effects.  Suggested to use the same window size as the original search.
#' @param dropSeqToolsLocation Location where dropseq tools are installed.
#' @param pythonPath Path to the python executable.  This requires python 3 to run.  If python is on your path, this can be "python".
#' @import MatrixEQTL
#' @export
findSecondaryEqtlEffects<-function (indexSNPFile, snpByGeneResultFile, geneLocationFile, expressionFile, variantLocationFile, genotypeFile, outFile, outFileEigenMT, outPDF=NULL, fdrThreshold=0.1, qvalueThreshold=0.05, cisDist=1000000,
                                    dropSeqToolsLocation="/broad/mccarroll/software/dropseq/priv", pythonPath="/broad/software/free/Linux/redhat_7_x86_64/pkgs/python_3.9.2/bin/python") {

    #after reading in the data, validate the column names.
    validateInputsExist(genotypeFile, geneLocationFile, variantLocationFile, expressionFile, indexSNPFile, snpByGeneResultFile)
    eGenes=data.table::fread(indexSNPFile)
    validateDFColumns(indexSNPFile, eGenes, expectedColumns=c("SNP", "gene", "gene_start_pos", "gene_end_pos", "TESTS", "qvalue"))
    eGenes[,chr:=getContigFromSNP(SNP), by=rownames(eGenes)]

    exp=data.table::fread(expressionFile)
    validateDFColumns(expressionFile, exp, expectedColumns=c("id"))

    varLoc=data.table::fread(variantLocationFile, data.table = F)
    validateDFColumns(variantLocationFile, varLoc, expectedColumns=c('snp', 'chr', 'pos'))
    varLoc=varLoc[,c('snp', 'chr', 'pos')]

    geneLoc=data.table::fread(geneLocationFile, data.table = F)
    validateDFColumns(geneLocationFile, geneLoc, expectedColumns=c('geneid', 'chr', 's1', 's2'))

    geno=data.table::fread(genotypeFile)
    validateDFColumns(genotypeFile, geno, expectedColumns=c("id"))

    #set up a temporary directory to write outputs.
    rootOutDir=paste(dirname (outFile), "/secondary_scan", sep="")
    if (!dir.exists(rootOutDir)) dir.create(rootOutDir)

    #eGenes=eGenes[eGenes$chr=="chr1",]
    genes=eGenes[eGenes$qvalue<=qvalueThreshold,]$gene
    if (length(genes)==0)
      stop (paste("0 genes significant at qvalue <="), qValueThreshold)

    #this is much faster if it done by contig.
    contigs=sort(unique (varLoc$chr))
    #contigs=contigs[1:10]
    #contigs=c("chrX")
    r=lapply(contigs, findSeondaryEffectByContig, genes, eGenes, exp, varLoc, geno, geneLoc, rootOutDir=NULL, cisDist=cisDist, fdrThreshold=fdrThreshold, verbose=T)
    r=do.call(rbind, r)

    #this is the slower way.
    #r=do.call(rbind, lapply(genes, findSecondaryEffect, eGenes, exp, varLoc, geno, cisDist=cisDist, fdrThreshold=fdrThreshold, verbose=T))

    annotateOutput(r, outFile, genotypeFile, variantLocationFile, expressionFile, geneLocationFile)

    cat(paste("[", length(unique(r$gene)), "] of [", length(genes), "] have a secondary QTL with an FDR <= [", fdrThreshold, "]\n"))
    #debug(runEigenMT)
    runEigenMT(indexSNPFile, eQTLFile=outFile, genotypeFile, variantLocationFile, geneLocationFile, outFileEigenMT, cisDist=cisDist, dropSeqToolsLocation, pythonPath)

    #plots!
    if (!is.null(outPDF)) plotQC (outFileEigenMT, outFile, outPDF, snpByGeneResultFile, expressionFile, variantLocationFile, genotypeFile, qvalueThreshold=qvalueThreshold, cisDist=cisDist)

}

validateEgeneDF<-function (indexSNPFile, eGenes) {
    expectedColumns=c("SNP", "gene", "gene_start_pos", "gene_end_pos", "TESTS", "qvalue", "foo")
    validateDFColumns("eGenes", )
}

validateDFColumns<-function (inFile, df, expectedColumns) {
    idxProblem=which(!expectedColumns %in% colnames(df))
    if (length(idxProblem)>0)
        stop (paste(inFile, "data frame missing column", expectedColumns[idxProblem]))
}

validateInputsExist<-function (genotypeFile, geneLocationFile, variantLocationFile, expressionFile, indexSNPFile, snpByGeneResultFile) {
  DropSeq.eqtl::validateSingleFileExists(genotypeFile)
  DropSeq.eqtl::validateSingleFileExists(geneLocationFile)
  DropSeq.eqtl::validateSingleFileExists(variantLocationFile)
  DropSeq.eqtl::validateSingleFileExists(expressionFile)
  DropSeq.eqtl::validateSingleFileExists(indexSNPFile)
  DropSeq.eqtl::validateSingleFileExists(snpByGeneResultFile)

}

getContigFromSNP<-function (snpID) {
    strsplit (snpID, ":", fixed=T)[[1]][1]
}


findSeondaryEffectByContig<-function (contig, genes, eGenes, exp, varLoc, geno, geneLoc, rootOutDir=NULL, cisDist=1000000, fdrThreshold=0.05, corThreshold=0.5, verbose=F) {
    cat (paste("Contig [", contig, "]\n", sep=""))
    idx=which(varLoc$chr==contig)

    varLocChr=varLoc[idx,]
    genoChr=geno[idx,]
    eGenesChr=eGenes[eGenes$chr==contig,]
    #additionally filter on the genes that were decided to be significant.
    eGenesChr=eGenesChr[eGenesChr$gene %in% genes,]

    genesChr=sort(unique (eGenesChr$gene))
    r=do.call(rbind, lapply(genesChr, findSecondaryEffect, eGenes=eGenesChr, exp, varLoc=varLocChr,
                            geno=genoChr, geneLoc=geneLoc, cisDist=cisDist, fdrThreshold=fdrThreshold, verbose=verbose))

    if (!is.null(rootOutDir)) {
      oFile<-paste(rootOutDir, "/secondary_scan.", contig, ".txt", sep="")
      write.table(r, oFile, row.names=F, col.names = T, quote=F, sep="\t")
    }
    return (r)

}

#geneName="AC016394.1"
#geneName="AC013652.1"
#geneName="WDR45"
#geneName="APIP"
#corThreshold is how correlated the index SNP is to any other SNP
findSecondaryEffect<-function (geneName, eGenes, exp, varLoc, geno, geneLoc, cisDist=1000000, fdrThreshold=0.05, corThreshold=0.5, verbose=F) {
    if (verbose) cat (geneName, "\n")
    eGene=eGenes[eGenes$gene==geneName,]

    #read in the genotypes for the index SNP
    indexSNPGenotype=as.numeric (geno[geno$id==eGene$SNP,-1])

    cvrt = MatrixEQTL::SlicedData$new();
    cvrt$CreateFromMatrix(matrix (indexSNPGenotype, nrow=1))

    startPos=eGene$gene_start_pos-cisDist
    endPos=eGene$gene_end_pos+cisDist
    chr=getContigFromSNP(eGene$SNP)

    #get variants in window, but leave out the eGene SNP
    vars=varLoc[varLoc$chr==chr & varLoc$pos>=startPos & varLoc$pos<=endPos,]
    vars=vars[-match(eGene$SNP, vars$snp),c('snp', 'chr', 'pos')]
    gene=geneLoc[geneLoc$geneid==geneName,]

    #NO SNPS TO TEST FOR GENE.  Happens if the eQTL index SNP was the only SNP in the search space.
    if (dim(vars)[1]==0) return (NULL)

    #RUN ON ONLY THE DONORS WITH COMPLETE GENOTYPES FOR THE LEAD INDEX SNP
    #zold=runMatrixEQTL (eGene, geno, exp, vars, gene, cvrt = cvrt,  dropMissingIndexSNPDonors=T, corThreshold=1.1)
    z=runMatrixEQTL (eGene, geno, exp, vars, gene, cvrt = cvrt,  cisDist, dropMissingIndexSNPDonors=T, corThreshold=0.5)

    if (is.null(z)) return (NULL)

    z$FDR=z$pvalue*eGene$TESTS
    idx=which(z$FDR>1)
    if (length(idx)>0) z[idx,]$FDR<-1

    if (min (z$FDR)<=fdrThreshold) {
        return (z)
    }
    return (NULL)

    #some hacky plots
    # e=as.numeric (exp[exp$id==geneName,-1])
    # g1=as.numeric (geno[geno$id==eGene$SNP,-1])
    # g2=as.numeric (geno[geno$id==z[which.min(z$pvalue),]$snps,-1])
    # g3=as.numeric (geno[geno$id==z[2,]$snps,-1])
    # Yes, the covariate does encode the G1 data:
    # g1==as.numeric (cvrt$FindRow("row1")$row)
    # summary (lm (e ~ g1 + g2))
    #

    # idxNotNA=which(!is.na(indexSNPGenotype))
    # e=as.numeric (exp[exp$id==geneName,(idxNotNA+1),with=F])
    # g1=as.numeric (geno[geno$id==eGene$SNP,(idxNotNA+1),with=F])
    # g2=as.numeric (geno[geno$id==z[which.min(z$pvalue),]$snps,(idxNotNA+1),with=F])
    # summary (lm (e ~ g1 + g2))

    #what about if there's missing data in the 2nd SNP?
    # idxNotNA=which(!is.na(g1) & !is.na(g2))
    # e=as.numeric (exp[exp$id==geneName,(idxNotNA+1),with=F])
    # g1=as.numeric (geno[geno$id==eGene$SNP,(idxNotNA+1),with=F])
    # g2=as.numeric (geno[geno$id==z[which.min(z$pvalue),]$snps,(idxNotNA+1),with=F])
    # summary (lm (e ~ g1 + g2))
}

#snpByGeneResultFile the full snp x gene output from the primary original scan
#outFile the full snp x gene output from the secondary scan
plotQC<-function (outFileEigenMT, outFile, outPDF, snpByGeneResultFile, expressionFile, variantLocationFile, genotypeFile, qvalueThreshold=0.05, cisDist=1000000) {

    getLD<-function (primaryIndexSNP, secondaryIndexSNP, geno, varLoc) {
        g1=as.numeric (geno[geno$id==primaryIndexSNP,-1])
        g2=as.numeric (geno[geno$id==secondaryIndexSNP,-1])
        g1=getGenotypeLabels(getSNPAlleles(primaryIndexSNP))[g1+1]
        g2=getGenotypeLabels(getSNPAlleles(secondaryIndexSNP))[g2+1]
        g1=genetics::as.genotype(g1)
        g2=genetics::as.genotype(g2)
        r=genetics::LD(g1, g2)
        return (r)

    }

    plotLinearRegressionSummary<-function (geneName, primaryIndexSNP, secondaryIndexSNP, exp, geno, varLoc, cex=1) {
        e=as.numeric (exp[exp$id==geneName,-1])
        g1=as.numeric (geno[geno$id==primaryIndexSNP,-1])
        g2=as.numeric (geno[geno$id==secondaryIndexSNP,-1])
        fit=summary(lm (e ~ g1 + g2))
        #if there aren't at least 2 coefficent outputs besides the intercept, get out of here.
        if (dim (fit$coefficients)[1]==2) {
            plot (0:1, 0:1, type='n', axes=F, xlab="", ylab="")
            text(x=0.5, y=0.8, "Problems running multivariate linear regression")
        } else {
            coof=fit$coefficients[-1,]
            #numeric formatting for plot
            coof[,4]=format(coof[,4], digits=3, scientific=T)
            for (i in 1:3)
              coof[,i]=format(as.numeric (coof[,i]), digits=3)

            #rownames (coof)=c(primaryIndexSNP, secondaryIndexSNP)
            rownames (coof)=c("SNP1", "SNP2")
            plot.new()
            plotrix::addtable2plot(-0.2,0.7,coof,bty="o",display.rownames=T,hlines=TRUE, cex=cex, xpad=0.25, ypad=0.75)
        }
        ld=getLD(primaryIndexSNP, secondaryIndexSNP, geno, varLoc)
        df=data.frame(r2=round(ld$r^2,3), pval=format(ld$`P-value`, digits=3, scientific = T))
        plotrix::addtable2plot(0.2,0.3,df,bty="o",display.colnames =T,hlines=TRUE, cex=cex, xpad=0.5, ypad=0.75, title="R2 of SNP Pair")
    }

    plotR2Summary<-function (geneName, primaryIndexSNP, secondaryIndexSNP, exp, geno) {
        e=as.numeric (exp[exp$id==geneName,-1])
        g1=as.numeric (geno[geno$id==primaryIndexSNP,-1])
        g2=as.numeric (geno[geno$id==secondaryIndexSNP,-1])
        r2G1=summary(lm (e ~ g1))$adj.r.squared
        r2G2=summary(lm (e ~ g2))$adj.r.squared
        r2G12=summary(lm (e ~ g1 + g2))$adj.r.squared
        ylim=c(0, max(c(r2G1, r2G2, r2G12)))
        barplot (c(r2G1, r2G2, r2G12), names.arg=c("SNP1", "SNP2", "SNP1+2"), ylab="r2 expression", ylim=ylim)
    }


    #2 page plotting, lays out the regions and SNP x Gene interactions on two pages.
    # use 4 different snp x gene interactions with expression data normalized in one version.
    #geneName="ICMT"
    plotTwo<-function (geneName, a, originalFull, newFull, exp, varLoc, geno, cisDist) {
        cat (paste(geneName, "\n"))
        eGene=a[a$gene==geneName,]
        #eGene=eGene[which.min(eGene$`p-value`),]
        startPos=unique (eGene$gene_start_pos-cisDist)
        endPos=unique (eGene$gene_end_pos+cisDist)
        chr=strsplit (eGene$SNP, ":", fixed=T)[[1]][1]

        zOriginal=originalFull[originalFull$gene==geneName,]
        zNew=newFull[newFull$gene==geneName,]

        primaryIndexSNP=zOriginal[which.min (zOriginal$`p-value`),]$SNP
        secondaryIndexSNP=zNew[which.min (zNew$`p-value`),]$SNP

        #CHECK THE LINEAR REGRESSION with both SNPs in plotPrimarySecondarySNPByGene
        #WHY does the slope look flat?
        old.par <- par(no.readonly = TRUE)
        layout(mat=matrix(c(1,2), ncol=1, byrow=T))
        parRegion=c(3,4,2,0)
        parSNP=c(3,4,2,1)
        plotRegionSecondaryEffects(zOriginal, eGene, danCistDist=cisDist,  ylab="p-value [-log10]", strTitleSuffix=paste("\nq-value [", format (eGene[eGene$SNP==primaryIndexSNP,]$qvalue, digits = 3, scientific = T), "]"), cex=0.5)
        plotRegionSecondaryEffects(zNew, eGene, danCistDist=cisDist, ylab="p-value after removing lead index SNP [-log10]", strTitleSuffix=paste("\nq-value [", format (eGene[eGene$SNP==secondaryIndexSNP,]$qvalue, digits = 3, scientific = T), "]"), cex=0.5)
        layout(mat=matrix(c(1,2,3,4,5,6), ncol=2, byrow=T))
        par(mar=parSNP)
        plotPrimarySecondarySNPByGene (primaryIndexSNP, secondaryIndexSNP, geneName, eGene, exp, varLoc, geno, regressSecondarySNPExpression=F)
        plotPrimarySecondarySNPByGene (primaryIndexSNP, secondaryIndexSNP, geneName, eGene, exp, varLoc, geno, regressSecondarySNPExpression=T)
        plotPrimarySecondarySNPByGene (secondaryIndexSNP, primaryIndexSNP, geneName, eGene, exp, varLoc, geno, regressSecondarySNPExpression=F)
        plotPrimarySecondarySNPByGene (secondaryIndexSNP, primaryIndexSNP, geneName, eGene, exp, varLoc, geno, regressSecondarySNPExpression=T)
        plotLinearRegressionSummary (geneName, primaryIndexSNP, secondaryIndexSNP, exp, geno, varLoc, cex=1.1)
        plotR2Summary (geneName, primaryIndexSNP, secondaryIndexSNP, exp, geno)

        #plot the primary index SNP genotypes colored by the secondary
        par(old.par)
        return (NULL)

    }

    #debug (plotOne)
    #sapply(names(genes), plotOne, a, originalFull, newFull, exp, varLoc, geno)

    #get the genes that have more than one q-value < threshold

    a=data.table::fread(outFileEigenMT)
    a=a[a$qvalue<=qvalueThreshold,]
    genes=which (table (a$gene)>1)
    genes=names (genes)
    if (length(genes)<1) {
        cat ("No new eQTLs passed qvalue threshold, not generating plots \n")
        return (NULL)
    }
    originalFull=data.table::fread(snpByGeneResultFile)
    newFull=data.table::fread(outFile)

    exp=data.table::fread(expressionFile)
    varLoc=data.table::fread(variantLocationFile, data.table = F)
    geno=data.table::fread(genotypeFile)

    #genes=genes[1:10]
    #genes="AHSG";

    #debug (plotPrimarySecondarySNPByGene)

    if (!is.null(outPDF)) pdf(outPDF)
    z=sapply(genes, plotTwo, a, originalFull, newFull, exp, varLoc, geno, cisDist)
    if (!is.null(outPDF)) dev.off()

}

#this is driven by the variants in the vars dataframe.
#corThreshold drop genotypes that are correlated to the index SNP > corThreshold
runMatrixEQTL<-function (eGene, geno, exp, vars, gene, cvrt = NULL, cisDist, dropMissingIndexSNPDonors=F, corThreshold=0.5) {
    if (dropMissingIndexSNPDonors) {
        indexSNPGenotype=as.numeric (cvrt$FindRow("row1")$row)
        idxNotNA=which(!is.na(indexSNPGenotype))
        cvrt = MatrixEQTL::SlicedData$new();
        cvrt$CreateFromMatrix(matrix (indexSNPGenotype[idxNotNA], nrow=1))
        genoThis=geno[match(vars$snp, geno$id),c(1, (idxNotNA+1)),with=F]
        expThis=exp[exp$id==eGene$gene,c(1, (idxNotNA+1)),with=F]
    } else {
        if (is.null(cvrt)) cvrt = MatrixEQTL::SlicedData$new();
        genoThis=geno[match(vars$snp, geno$id),]
        expThis=exp[exp$id==eGene$gene,]
    }

    #if (eGene$gene=="AC004656.1") browser()

    #find correlated SNPs to lead SNP
    indexGenotype=as.numeric (cvrt$FindRow("row1")$row)
    varCor= findVariantsCorrelatedToLead(indexGenotype, genoThis)
    idx=match(varCor[abs(varCor$V1)>=corThreshold,]$id, vars$snp)
    if (length(idx)>0) {
        vars=vars[-idx,]
        genoThis=genoThis[-idx,]
    }

    #if there are no genotypes left, return NULL
    if (dim (genoThis)[1]==0) return (NULL)

    #convert to matrix eQTL format
    genoThis=DropSeq.eqtl::convertToSlicedData(genoThis)
    expThis=DropSeq.eqtl::convertToSlicedData(expThis)

    invisible(capture.output(me<- MatrixEQTL::Matrix_eQTL_main(
        snps = genoThis,
        gene = expThis,
        cvrt = cvrt,
        pvOutputThreshold     = 0,
        useModel = MatrixEQTL::modelLINEAR,
        errorCovariance = numeric(),
        verbose = F,
        pvOutputThreshold.cis = 1,
        snpspos = vars,
        genepos = gene,
        cisDist = cisDist,
        pvalue.hist = FALSE,
        min.pv.by.genesnp = T,
        noFDRsaveMemory = F), type="message"))

    z=data.frame(me$cis$eqtls, stringsAsFactors = F)
    z$pos=vars[match(z$snps, vars$snp),]$pos
    z$dfFull=me$param$dfFull
    #z$corToLead=varCor[match(z$snps, varCor$id),]$V1
    return (z)
}

findVariantsCorrelatedToLead<-function (indexGenotype, genoThis) {
    f<-function (x, indexGenotype) cor(as.numeric(x), indexGenotype, use="complete")
    result=genoThis[,f(.SD, indexGenotype),by="id"]
    return (result)
}

annotateOutput<-function (r, outFile, genotypeFile, variantLocationFile, expressionFile, geneLocationFile) {
    tempCisFile=tempfile(fileext=".cis.txt")
    tempParamFile=tempfile(fileext=".params.txt")
    #it's possible that different tests have different numbers of samples, so penalize the maximum number.
    params=data.frame(dfFull=max(unique (r$dfFull)))

    r=cleanupOutputColumnNames(r)

    write.table(params, tempParamFile, row.names=F, col.names = T, quote=F, sep="\t")
    write.table(r, tempCisFile, row.names=F, col.names = T, quote=F, sep="\t")

    #need to construct the degrees of freedom as a file.
    #params$dfFull


    DropSeq.eqtl::annotateEQTL(output_file_name_cis=outFile, cis_eqtl_file_name=tempCisFile, params_file_name=tempParamFile, SNP_file_name=genotypeFile,
                               snps_location_file_name=variantLocationFile, expression_file_name=expressionFile,
                               gene_location_file_name=geneLocationFile)
    file.remove(tempCisFile)
    file.remove(tempParamFile)
}

cleanupOutputColumnNames<-function (r) {
    colnames(r)[1]="SNP"
    colnames(r)[3]="t-stat"
    colnames(r)[4]="p-value"
    r[c(1,2,6,3,4)]
}



#z is the output of matrixEQTL $cis$eqtls
plotRegionSecondaryEffects<-function (z, eGene, danCistDist=10000, strTitleSuffix="", ...) {
    ylim=c(0, max(-log10(c(z$`p-value`, eGene$`p-value`))))
    plot (z$snp_pos/1000, -log10(z$`p-value`), xlab="position [kb]", ylim=ylim, ...)
    abline (v=(eGene$gene_start_pos-danCistDist)/1000, col='red', lty=2)
    abline (v=(eGene$gene_end_pos+danCistDist)/1000, col='red', lty=2)
    title(paste(unique (z$gene), strTitleSuffix))
}

getGenotypeLabels<-function (snpAlleles) {
    c(paste(snpAlleles$ref, snpAlleles$ref, sep="/"),
      paste(snpAlleles$ref, snpAlleles$alt, sep="/"),
      paste(snpAlleles$alt, snpAlleles$alt, sep="/"))
}

#indexSNP=
plotPrimarySecondarySNPByGene<-function (primaryIndexSNP, secondaryIndexSNP, geneName, eGene, exp, varLoc, geno, regressSecondarySNPExpression=F) {

    plotPrimaryVsSecondary<-function (primaryIndexSNP, secondaryIndexSNP, e, geno, regressSecondarySNPExpression=F) {
        gP=as.numeric (geno[geno$id==primaryIndexSNP,-1])
        gS=as.numeric (geno[geno$id==secondaryIndexSNP,-1])
        cols=c("red", "blue", "green")
        snpAllelesP=getSNPAlleles(primaryIndexSNP)
        labelsP=getGenotypeLabels(snpAllelesP)

        eFinal=e
        ylab="expression"
        if (regressSecondarySNPExpression) {
            fit=lm (e ~ gS, na.action=na.exclude)
            eFinal=residuals(fit) + mean (e)
            ylab="fitted expression"
        }

        ylim=range(eFinal, na.rm=T)
        ylim[2]=ylim[2]*1.2

        #shift the outer edge genotypes slightly towards the center.
        centers=c(0.1, 1, 1.9)
        gpFixed=gP
        gpFixed[gpFixed==0]<-centers[1]
        gpFixed[gpFixed==2]<-centers[3]
        #restrict jitter to bounds.
        gPJitter=jitter(gpFixed, factor=0.4)
        gPJitter[gPJitter>2]<-1.99
        gPJitter[gPJitter<0]<-0.01
        plot (eFinal ~ gPJitter, xlab="", ylab=ylab, axes=F, col=cols[gS+1], pch=16, ylim=ylim, xlim=c(0,2))
        axis(1, at=centers, labels=labelsP)
        axis(2)
        title (paste(primaryIndexSNP))
        abline (lm(eFinal ~ gP))
        labelsS=getGenotypeLabels(getSNPAlleles(secondaryIndexSNP))
        legend('topleft', legend=labelsS, ncol=3, fill=cols)
    }

    #eGene=eGenes[eGenes$gene==geneName,]
    e=as.numeric (exp[exp$id==unique (eGene$gene),-1])

    plotPrimaryVsSecondary (primaryIndexSNP, secondaryIndexSNP, e, geno, regressSecondarySNPExpression)
    #plotPrimaryVsSecondary (secondaryIndexSNP ,primaryIndexSNP, e, geno)

}

#eQTLFile=outFile
runEigenMT<-function (indexSNPFile, eQTLFile, genotypeFile, variantLocationFile, geneLocationFile, outFileEigenMT, cisDist=10000,
                      dropSeqToolsLocation, pythonPath="/Users/nemesh/opt/anaconda3/bin/python") {

    eigenMTCommand=paste(pythonPath, " ", dropSeqToolsLocation, "/3rdParty/eigenMT/eigenMT.py", sep="")
    a=data.table::fread(eQTLFile)
    contigs=sort(unique (sapply(strsplit (a$SNP, ":", fixed=T), function (x) x[1])))

    #TODO: get rid of this.
    #contigs=contigs[1]

    #make a new directory for eigenMT secondary scan outputs.
    rootDir=paste(dirname (outFileEigenMT), "/secondary_scan", sep="")
    if (!dir.exists(rootDir)) dir.create(rootDir)

    getOutFile<-function (chrom) paste(rootDir, "/eigenMT.", chrom, sep="")
    outFileList=sapply(contigs, getOutFile)

    for (index in 1:length(contigs)) {
        contig=contigs[index]
        cat ("Running eigenMT on contig", contig, "\n")
        o=outFileList[index]
        cmd=paste(eigenMTCommand, "--QTL", eQTLFile, "--GEN", genotypeFile, "--GENPOS", variantLocationFile,
                  "--PHEPOS", geneLocationFile, "--OUT ", o, "--var_thresh 0.99 --window 200 --cis_dist",
                  format(cisDist, scientific=F), "--CHROM", contig)
        cat (paste(cmd, "\n"))
        system(cmd)

        #mimic fix_tsv, which I'm having trouble calling.
        b=data.table::fread(as.character (o))
        data.table::setnames(b, "BF", "gene_permuted_pvalue")
        write.table(b, o, row.names=F, col.names=T, quote=F, sep="\t")
    }

    #need to add back in the first pass eQTL scan
    tempScan=tempfile(".original.txt")
    z=data.table::fread(indexSNPFile)
    #remove the qvalue column, move the gene_permuted_pvalue column back to "BF"
    z[,c("qvalue"):=NULL]

    write.table(z, tempScan, row.names=F, col.names = T, quote=F, sep="\t")

    outFileList=c(outFileList, tempScan)
    DropSeq.eqtl::combinePermutationResults(outFileList, outFileEigenMT)
    sapply(outFileList, file.remove)

}

getSNPAlleles<-function (var) {
    z=strsplit (var, split=":")[[1]][3:4]
    data.frame(ref=z[1], alt=z[2])
}

plotNumTestsViaEigenMTPerGene<-function () {
    a=data.table::fread(outFileEigenMT)
    a=a[a$qvalue<=qvalueThreshold,]
    genes=which (table (a$gene)>1)
    genes=names (genes)

    #get results with 2 hits.
    a=a[which(!is.na(match(a$gene, genes))),]
    a=a[order(a$gene),]

    #x=a[a$gene=="AAMDC",]
    getTestsPerGene<-function (x) {
        best=which.min(x$`p-value`)
        second=setdiff(c(1,2), best)
        data.frame(snp1_tests=x[best,]$TESTS, snp2_tests=x[second,]$TESTS, stringsAsFactors = F)
    }


    r=a[,getTestsPerGene(.SD),by="gene"]
    plot (r$snp1_tests, r$snp2_tests, xlab="Number of tests primary index SNP", ylab="Number of tests secondary SNP", main="EigenMT number of tests")
    plot (table(r$snp1_tests-r$snp2_tests), xlab="difference in number of tests [primary - secondary]", ylab="number of eGenes")

    getEffectDirection<-function (x) {
        best=which.min(x$`p-value`)
        second=setdiff(c(1,2), best)
        data.frame(snp1_tests=x[best,]$beta, snp2_tests=x[second,]$beta, stringsAsFactors = F)
    }

    r=a[,getEffectDirection(.SD),by="gene"]
    rr=r


}

# foo<-function () {
#     z=c(0.9, 0.8, 0.7, 0.8, 0.5, 0.99, 0.4, 0.2)
#     mean (z)
#     s=sum (log10(z))
#     10^(s/length(z))
#
#     z=runif(10000)*0.75
#     mean (z)
#     s=sum (log10(z))
#     10^(s/length(z))
#
#     z=runif(10000)*0.5
#     mean (z)
#     s=sum (log10(z))
#     10^(s/length(z))
#
#
#
# }
