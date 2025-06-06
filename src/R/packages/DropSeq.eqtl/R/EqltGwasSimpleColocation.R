# windowSize=100000
# gwasHitThreshold=5e-8
# eQTLEmpiricThreshold=1e-4
#
# gwasFile="/downloads/GWAS_correlation/daner_PGC_SCZ_1016_plus_INDEX_SNPS.txt"
# permutedQTLFile="/downloads/cell_selection/eQTL/d21Ngn2plusGliaLive.exonic+intronic.whole_cell_plus_nuclei.maf_0.20_cisDist_10kb.eQTL_permuted_results.txt"
# unfilteredQTLFile="/downloads/cell_selection/eQTL/d21Ngn2plusGliaLive.exonic+intronic.whole_cell_plus_nuclei.maf_0.05_cisDist_10kb.eQTL_results.txt"
# outFile=paste("/downloads/cell_selection/eQTL/d21Ngn2plusGliaLive.exonic+intronic.whole_cell_plus_nuclei.gwas_hits.", windowSize/1000, "kb.txt", sep="")
#
# library (data.table)

#'  For all GWAS hits above some pvalue threshold, find GWAS SNPs near or within a gene of interest.
#'
#' 	Take the GWAS SNP hit and look up all nearby genes within some window size.
#'  Output the maximum qvalue gene within this region, as well as the empiric p-value for the SNP if tested.
#'  Order results by q-values.
#'
#' @param gwasFile A file containing GWAS SNP Hits
#' @param permutedQTLFile A file containing permuted eQTL hits
#' @param unfilteredQTLFile A file containing all unfiltered eQTL SNP tests.
#' @param outFile Output file.
#' @param gwasHitThreshold The threshold pvalue to select a GWAS SNP as a hit.
#' @param eQTLEmpiricThreshold The threshold pvalue to select an eQTL SNP as a hit.
#' @param windowSize How far around genes to search for an eQTL hit
#' @return A list of genes with both a GWAS and EQTL hit.
#' @export
findGwasHits<-function (gwasFile, permutedQTLFile, unfilteredQTLFile, outFile, gwasHitThreshold=5e-8, eQTLEmpiricThreshold=1e-5, windowSize=10000) {
    g=fread(gwasFile, header=T, stringsAsFactors = F)
    g=g[g$P<=gwasHitThreshold,]

    pQTL=parseEQTLData(permutedQTLFile)
    uQTL=parseEQTLData(unfilteredQTLFile)
    #gg=g[1:100,]
    #r=gg[,getOneSNPResult(.SD, pQTL, uQTL, windowSize),by=rownames(gg)]
    chromosomes=unique (g$CHR)
    system.time (r<-do.call(rbind, lapply(chromosomes, getSNPResultsPerChromosome, g, pQTL, uQTL, windowSize)))
    #system.time (r<-g[,getOneSNPResult(.SD, pQTL, uQTL, windowSize),by=rownames(g)])
    r[,rownames:=NULL]
    r=r[order(r$qvalue, r$gwas_pvalue, r$gene),]
    write.table(r, outFile, row.names=F, col.names=T, quote=F, sep="\t")

    #some interesting SNPs?

    cols=c("blue", "green")
    plot (-log10(r$gwas_pvalue), -log10(r$eQTL_empiric_P), xlab="GWAS pvalue [-log10]", ylab="eQTL empiric pvalue [-log10]", col=cols[as.numeric (factor (is.na(r$qvalue)))], cex=0.5)
    legend("topright", legend=c("permutation tested", "not permutation tested"), fill=cols)
    title(paste("eQTLs within [", windowSize/1000, "] kb of GWAS hits", sep=""))
    abline (h=-log10(eQTLEmpiricThreshold), lty=2, col="grey")

    #how many GWAS snps were tested directly?
    format(dim (r[!is.na(r$eQTL_empiric_P),])[1]/dim (r)[1]*100, digits=3, scientific=F)

    #how many genes were tested directly?
    length(unique (r$gene))
    length(unique(pQTL$gene))

    genes1=unique (r[r$eQTL_empiric_P<eQTLEmpiricThreshold & !is.na(r$gene)]$gene)
    r[r$gene==genes1[1],]

    #some eQTL hits near GWAS hits
    unique (r[r$qvalue<0.05,]$gene)

    #

}

#run one chromosome at a time.
getSNPResultsPerChromosome<-function (chr, g, pQTL, uQTL, windowSize=10000) {
    p=pQTL[pQTL$chromosome==chr,]
    u=uQTL[uQTL$chromosome==chr,]
    gg=g[g$CHR==chr,]
    system.time (r<-gg[,getOneSNPResult(.SD, p, u, windowSize),by=rownames(gg)])
    return (r)
}


#snp=g[1,]
#snp=g[g$SNP=="rs999867",]
getOneSNPResult<-function (snp, pQTL, uQTL, windowSize=10000) {
    #genes within windowSize
    genes=pQTL[pQTL$chromosome==snp$CHR & pQTL$gene_start_pos>=(snp$POS-windowSize) & pQTL$gene_end_pos<=(snp$POS+windowSize),]
    genes=genes[order(genes$qvalue),]
    bestGene=genes[1,]
    empiricEQTLP=uQTL[uQTL$chromosome==snp$CHR & uQTL$pos==snp$POS,]$`p-value`
    if (length(empiricEQTLP)==0) empiricEQTLP=NA_real_
    if (length(empiricEQTLP)>1) empiricEQTLP=min(empiricEQTLP)
    result=data.frame(chr=snp$CHR, pos=snp$POS, snpID=snp$SNP, a1=snp$A1, a2=snp$A2, gwas_pvalue=snp$P, gwas_OR=snp$OR, eQTL_empiric_P=empiricEQTLP,
                      gene=bestGene$gene, gene_start=bestGene$gene_start_pos, gene_end=bestGene$gene_end_pos, qvalue=bestGene$qvalue, stringsAsFactors = F)
    return (result)
}

parseEQTLData<-function (eQTLFile) {
    #for R CMD CHECK
    chromosome=NULL;SNP=NULL;pos=NULL;

    getChromosome<-function (x) strsplit (x, ":")[[1]][1]
    getPosition<-function (x) as.numeric(strsplit (x, ":")[[1]][2])
    x=fread(eQTLFile)
    x[,chromosome:=getChromosome(SNP),by=rownames(x)]
    x[,pos:=getPosition(SNP),by=rownames(x)]
    setkey(x, chromosome, pos)
    return (x)
}

