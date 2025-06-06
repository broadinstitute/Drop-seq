# 1) Read in the input score file (1.Zscores.txt.gz) and the manifest.
# 2) Get the genomic region of interest (chr:start-end)
# 3) Query the VCF for the LD block of the region.  Output as a .gzip file.
# 4) Read in the block and the gwas/eQTL scores file, and look up the IDs vs the input score file
# a. Intersect the SNPs
# 5) Output for the intersect of SNPs:
#     a. LD block file
# b. GWAS scores file
# c. eQTL scores file
# 6) Run ecaviar on 4a,4b,4c.  OR caviar on 4a,4b.
# 7) How do we gather up the summary info from eCaviar?  Maybe plotting where we color in the 95% confidence SNPs, color in the co-located SNP(s)?
#     a. Do we need output file(s) to support LDScore regression later?
#     i. 95% confidence set
# ii. 95% confidence CPP from caviar

# library (data.table);library (DropSeq.utilities)
#
#
# inSNPFile="/downloads/GWAS_correlation/CAVIAR_FILES/snpBlocks_daner_NGN2/115.Zscores.txt.gz"
# manifestFile="/downloads/GWAS_correlation/CAVIAR_FILES/snpBlocks_daner_NGN2/gene_manifest.txt"
# outDir="/downloads/GWAS_correlation/CAVIAR_FILES/snpBlocks_daner_NGN2/out"
# vcfFile="/downloads/GWAS_correlation/EUR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
# gene_location_file_name="/downloads/d21NGNG2/d21NGN2/eQTL/d21NGN2_MAF_0.2.gene_locations.txt"
# ldScoreBcftoolsPluginLocation="/home/unix/nemesh/bin/LDmatrix.so"
# eCAVIARPath="/broad/mccarroll/nemesh/caviar/CAVIAR-C++/CAVIAR"
# bcfToolsLocation="/Users/nemesh/bin/bcftools"
# ldScoreBcftoolsPluginLocation="/Users/nemesh/bin/LDmatrix.so"
# CAVIARPath="/Users/nemesh/caviar/CAVIAR-C++/CAVIAR"
# eCAVIARPath="/Users/nemesh/caviar/CAVIAR-C++/eCAVIAR"
# mode="eCAVIAR";dryRun=F; verbose=T; zScoreAbsValues=F;


#eQTLPermuatedFile="/downloads/d21NGNG2/d21NGNG2_MAF_0.2.eQTL_permuted_results.txt"; eQTLLargeWindowFile="/downloads/d21NGNG2/d21NGNG2/data/d21NGNG2_MAF_0.05.eQTL_results_250kb.txt"; gwasFile="/downloads/GWAS_correlation/ckqny.scz2snpres.gz";outFileDir="/downloads/GWAS_correlation/CAVIAR_FILES"; qValueThreshold=0.05; outManifestFile=NULL; verbose=T; geneListFile=NULL
#eQTLPermuatedFile="/downloads/d21NGNG2/d21NGN2/eQTL/d21NGN2_MAF_0.2.eQTL_permuted_results.txt"; eQTLLargeWindowFile="/downloads/d21NGNG2/d21NGN2/data/d21NGNG2_MAF_0.05.eQTL_results_250kb.txt"; gwasFile="/downloads/GWAS_correlation/ckqny.scz2snpres.prepped.txt";outFileDir="/downloads/GWAS_correlation/CAVIAR_FILES"; qValueThreshold=0.05; outManifestFile=NULL; verbose=T; geneListFile=NULL

#eQTLPermuatedFile="/downloads/d21NGNG2/d21NGNG2_MAF_0.2.eQTL_permuted_results.txt"; eQTLLargeWindowFile="/downloads/d21NGNG2/d21NGNG2/data/d21NGNG2_MAF_0.05.eQTL_results_250kb.txt"; gwasFile="/downloads/GWAS_correlation/daner_PGC_SCZ_1016_plus_INDEX_SNPS.txt"; outFileDir="/downloads/GWAS_correlation/CAVIAR_FILES"; qValueThreshold=0.05; outManifestFile=NULL; verbose=T; geneListFile=NULL


#UNION TEST
#eQTLPermuatedFile="/downloads/d21NGNG2/d21NGN2/eQTL/d21NGN2_MAF_0.2.eQTL_permuted_results.txt"; eQTLLargeWindowFile="/downloads/d21NGNG2/d21NGN2/data/d21NGNG2_MAF_0.05.eQTL_results_250kb.txt"; gwasFile="/downloads/GWAS_correlation/ckqny.scz2snpres.prepped.txt";outFileDir="/downloads/GWAS_correlation/CAVIAR_FILES"; qValueThreshold=0.05; outManifestFile=NULL; verbose=T; geneListFile=NULL; mode="union"
#eQTLPermuatedFile="/downloads/d21NGNG2/d21NGN2/eQTL/d21NGN2_MAF_0.2.eQTL_permuted_results.txt"; eQTLLargeWindowFile="/downloads/d21NGNG2/d21NGN2/data/d21NGNG2_MAF_0.05.eQTL_results_250kb.txt"; gwasFile="/downloads/GWAS_correlation/daner_PGC_SCZ_1016_plus_INDEX_SNPS.txt"; outFileDir="/downloads/GWAS_correlation/CAVIAR_FILES"; qValueThreshold=0.05; gwasPValueThreshold=1e-6; outManifestFile=NULL; verbose=T; geneListFile="/downloads/GWAS_correlation/NGN2_plus_MN_genelist.txt"; mode="union"

#THE INPUT / OUTPUT FILES
getGWASFile<-function (outputBaseFileName) paste(outputBaseFileName, ".gwas.Z.txt", sep="")
geteQTLFile<-function (outputBaseFileName) paste(outputBaseFileName, ".eQTL.Z.txt", sep="")
getBothGWASFile<-function (outputBaseFileName) paste(outputBaseFileName, ".both.gwas.Z.txt", sep="")
getBotheQTLFile<-function (outputBaseFileName) paste(outputBaseFileName, ".both.eQTL.Z.txt", sep="")

getLDScoreMatrix<-function (outputBaseFileName) paste(outputBaseFileName, ".both.LD.txt", sep="")
getLDScoreMatrixGWAS<-function (outputBaseFileName) paste(outputBaseFileName, ".gwas.LD.txt", sep="")
getLDScoreMatrixeQTL<-function (outputBaseFileName) paste(outputBaseFileName, ".eQTL.LD.txt", sep="")
getLDFullRegion<-function (outputBaseFileName) paste(outputBaseFileName, ".LD.region.txt.gz", sep="")
getECaviarJointResult<-function (outputBaseFileName) paste(outputBaseFileName, ".both_col", sep="")
getConfSetEQTL<-function (outputBaseFileName) paste(outputBaseFileName, ".eQTL_set", sep="")
getConfSetGWAS<-function (outputBaseFileName) paste(outputBaseFileName, ".GWAS_set", sep="")
confSetBothEQTL<-function (outputBaseFileName) paste(outputBaseFileName, ".both_1_set", sep="")
confSetBothGWAS<-function (outputBaseFileName) paste(outputBaseFileName, ".both_2_set", sep="")
getPDFOutput<-function (outputBaseFileName) paste(outputBaseFileName, ".pdf", sep="")
getSummaryFileOneGene<-function(outputBaseFileName) paste(outputBaseFileName,".summary.txt", sep="")

#' For each eQTL that has a q-value <= qValueThreshold, get the gwas results and eQTL results and build a dataframe of
#' SNP results for the intersecting SNP set.
#' Calculates Z-scores for each SNP in both data sets.
#' @param eQTLPermuatedFile the permuted eQTL results containing the q-value.
#' @param eQTLLargeWindowFile the empiric eQTL results at a large window size (250kb around gene)
#' @param gwasFile The GWAS results file, containing the following columns:
#' chr, snpid, a1, a2, pos, p.  These are chromosome, snp name, allele 1, allele 2, and p-value.
#' @param qValueThreshold only look at genes with qvalues <= this threshold.
#' @param gwasPValueThreshold only look at blocks that have a GWAS SNP with a value <= this threshold.
#' @param geneListFile override the qValueThreshold and output snp blocks for all genes in the list.
#' @param mode If union, builds a data frame containing all eQTL and GWAS SNPs, with the common snps "in frame".  If intersect, only use the SNPs present in both data sets.
#' @param outFileDir The directory to deposit files.
#' @param outManifestFile The location to write the manifest file, which contains the gene name and file name (relative to the outFileDir) containing the SNP dataframe.
#' If set to null, outputs the file name "Gene_Manifest.txt" to the outFileDir
#' @param verbose Write progress report to standard output.
#' @export
buildGWASeQTLSNPSets<-function (eQTLPermuatedFile, eQTLLargeWindowFile, gwasFile, qValueThreshold=0.05, gwasPValueThreshold=1e-6, geneListFile=NULL, mode=c("intersect","union"), outFileDir, outManifestFile=NULL, verbose=F) {
    mode=match.arg(mode, c("intersect","union"))
    eQTLPermuted=getEQTLResults(eQTLPermuatedFile ,getBestGenePerSNP=F)
    eQTLPermuted=eQTLPermuted[order(eQTLPermuted$qvalue),]
    if (!is.null(geneListFile)) {
        gl=read.table(geneListFile, header=F, stringsAsFactors = F)$V1
        both=intersect(gl, eQTLPermuted$gene)
        cat (paste("Found", length(both), "genes in permuted eQTLs + gene list"))
        eQTLPermuted=eQTLPermuted[match(both, eQTLPermuted$gene),]
        #filter eQTL Permuted list by the gene list instead of by qValue.
    } else {
        eQTLPermuted=eQTLPermuted[eQTLPermuted$qvalue<=qValueThreshold,]
    }

    eQTL=getEQTLResults(eQTLLargeWindowFile, getBestGenePerSNP=F)

    #limit eQTL large window results to q-value thresholded permuted genes for efficiency
    idx=which(!is.na(match(eQTL$gene, eQTLPermuted$gene)))
    eQTL=eQTL[idx,]
    snpsEQTLBeforeFiltering=length(unique(eQTL$SNP))

    gwas=getGWasResults(gwasFile)
    if (mode=="intersect") {
        #intersect of eQTL and GWAS results
        r=intersectGwasEQTL(gwas, eQTL)
        gwas=r$gwas
        eQTL=r$eQTL
        snpsEQTLAfterFiltering=length(unique(eQTL$SNP))
        paste ("Number of eQTL SNPs after intersect with GWAS results [", snpsEQTLAfterFiltering, "]")
        paste("Fraction of eQTL SNPs also in GWAS results [", round(snpsEQTLAfterFiltering/snpsEQTLBeforeFiltering,3), "]")
    } else if (mode=="union") {
        #Union of eQTL and GWAS results.
        r=unionGwasEQTL(gwas, eQTL)
        gwas=r$gwas
        eQTL=r$eQTL
    }

    #for speed, convert the GWAS results to a data table
    setDT(gwas, key=c("CHR", "POS"))
    setDT(eQTL, key="gene")

    #run over all genes.
    geneList=eQTLPermuted$gene
    geneListFinal=c()
    #make room for the GWAS pvalue for the best SNP [minimum p-val] in the region.
    eQTLPermuted$GWAS_P=vector(mode="numeric", length=dim(eQTLPermuted)[1])
    #pb=txtProgressBar(min=1, max=length(geneList), title="Genes Processed", style=3)
    for (i in 1:length(geneList)) {
        #setTxtProgressBar(pb, i)
        if (verbose) cat (paste("Gene Processed [", i, "] of [", length(geneList), "]\n", sep=""))
        geneName=geneList[i]
        if (mode=="intersect") df=getGwasEQTLOverlappingSNPs(geneName, gwas, eQTL)
        if (mode=="union") df=getGwasEQTLUnionSNPs(geneName, gwas, eQTL, verbose=F)

        #check conditions to proceed: DF has SNPs, has at least 1 gwas SNP with pvalue < threshold.
        if (is.null(df) ) next;
        #if there are no GWAS SNPs, go to next block.
        numGWASSNPs=length(which(!is.na(df$gwas_P)))
        if (numGWASSNPs==0) next
        #test GWAS SNPs for minimum p-value.
        minGwasPval=min (df$gwas_P, na.rm=T)
        #if no suitably small gwas pvalue found, go to the next block.
        if (minGwasPval> gwasPValueThreshold) next

        #you made it, process and write out snp block.
        geneListFinal=c(geneListFinal, geneName)
        #update the eQTL permuted results with the GWAS minimum P.
        eQTLPermuted[eQTLPermuted$gene==geneName,]$GWAS_P=minGwasPval
        ii=length(geneListFinal) #the new i.
        outFile=gzfile(paste(outFileDir, "/", ii, ".Zscores.txt.gz",sep=""))
        write.table(df, outFile, row.names=F, col.names=T, quote=F, sep="\t")
    }

    #write manifest
    #include the Z and P values in the manifest.
    if (is.null(outManifestFile)) outManifestFile=paste(outFileDir, "/", "gene_manifest.txt", sep="")
    idx=match(geneListFinal, eQTLPermuted$gene)
    manifestDF=data.frame(iter=1:length(geneListFinal), gene=geneListFinal, gwas_P=format(eQTLPermuted[idx,]$GWAS_P,digits=4, scientific=T), eQTL_P=format(eQTLPermuted[idx,]$`p-value`, digits=4, scientific=T),  qvalue=format(eQTLPermuted[idx,]$qvalue, digits=4, scientific=T), stringsAsFactors = F)
    write.table(manifestDF, outManifestFile, row.names=F, col.names=T, quote=F, sep="\t")
}

# runMany(inSNPFilePath="/downloads/GWAS_correlation/CAVIAR_FILES", manifestFile="/downloads/GWAS_correlation/CAVIAR_FILES/gene_manifest.txt",
#        vcfFile="/downloads/GWAS_correlation/EUR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
#		 outDir="/downloads/GWAS_correlation/CAVIAR_FILES/out",
#		 bcfToolsLocation="/Users/nemesh/bin/bcftools", ldScoreBcftoolsPluginLocation="/Users/nemesh/bin/LDmatrix.so",
#		 CAVIARPath="/Users/nemesh/caviar/CAVIAR-C++/CAVIAR", eCAVIARPath="/Users/nemesh/caviar/CAVIAR-C++/eCAVIAR")

#' Title
#'
#' @param inSNPFilePath A directory containing one or more combined SNP data frames of eQTL + GWAS data, generated by buildGWASeQTLSNPSets
#' @inheritParams runPipeline
#' @export
runPipelineMany<-function (inSNPFilePath, manifestFile, vcfFile, outDir, bcfToolsLocation, ldScoreBcftoolsPluginLocation, CAVIARPath, eCAVIARPath, dryRun=F, verbose=F) {
	snpFiles=list.files(inSNPFilePath, pattern="Zscores.txt.gz", full.names = T)
	snpFiles=snpFiles[1:10]

	for (i in 1:length(snpFiles)) {
		cat(i, "\n")
		runPipeline(snpFiles[i], manifestFile, vcfFile, outDir,  bcfToolsLocation, ldScoreBcftoolsPluginLocation, CAVIARPath, eCAVIARPath, dryRun, verbose)
	}
}


#TODO filter out A/T, C/G SNPs that could have snp-flip effects that can't easily be detected?

#' For a data frame of SNPs in both eQTL and GWAS data, run the caviar fine mapping pipeline.
#'
#' @param inSNPFile A file containing a combined SNP data frame of eQTL + GWAS data, generated by buildGWASeQTLSNPSets
#' @param manifestFile A file containing each inSNPFile number, and the gene it corresponds to.
#' @param vcfFile A VCF file to extract LD data from for this region.  Suggest using a population from thousand genomes or similar.
#' @param gene_location_file_name A file name containing the coordinates of genes.  Has 4 columns: geneid, chr, s1, s2.  S1 is the start, S2 the end of the gene.
#' @param outDir Where should data be output to?
#' @param bcfToolsLocation The full path to the bcftools executable.
#' @param ldScoreBcftoolsPluginLocation The ful path to the ldScore plugin for bcftools.
#' @param CAVIARPath The full path to the CAVIAR executable
#' @param eCAVIARPath The full path to the eCAVIAR executable
#' @param dryRun Print out the commands that would be run for 3rd party tools, but don't execute anything.
#' @param verbose More logging during run.
#' @param zScoreAbsValues Use the unsigned Z scores for inputs to CAVIAR/eCAVIAR
#' @export
#' @seealso buildGWASeQTLSNPSets
runPipeline<-function (inSNPFile, manifestFile, vcfFile, gene_location_file_name, outDir, bcfToolsLocation, ldScoreBcftoolsPluginLocation, CAVIARPath, eCAVIARPath, dryRun=F, verbose=T, zScoreAbsValues=F) {
    #mode <- match.arg(mode)

    outputBaseFileName=getOutputFileNameBase(inSNPFile, outDir)
    df=read.table(inSNPFile, header=T, stringsAsFactors = F, sep="\t")
    #df=df[!is.na(df$eQTL_Z) & !is.na(df$gwas_Z),]
    #filter out A/T, C/G SNPs.

    ldBlock=buildLDBlock(df, vcfFile, outputBaseFileName, bcfToolsLocation, ldScoreBcftoolsPluginLocation, cleanUpFiles=F, dryRun, verbose=verbose)
    #ldBlock=buildLDBlockPlink(df, vcfFile, outputBaseFileName, ldScoreBcftoolsPluginLocation, cleanUpFiles=F)

    #get the GWAS/eQTL intersect SNPs, GWAS ONLY, EQTL ONLY

    #with the full region LD block in hand and a SNP data frame, intersect the two data sets.
    dfBoth=df[!is.na(df$eQTL_P) & !is.na(df$gwas_P),]
    resultBoth=intersectGwasAndLDBlock(dfBoth, ldBlock, checkReversedAlleles=F)
    dfBoth=resultBoth$df
    ldBlockBoth=resultBoth$ldBlock

    #with the full region LD block in hand and a GWAS SNP data frame
    #in this case we check for reversed alleles because GWAS only hits were not previously flipped
    #to the reference allele.
    dfGWAS=df[!is.na(df$gwas_P),]
    resultGWAS=intersectGwasAndLDBlock(dfGWAS, ldBlock, checkReversedAlleles=T)
    dfGWAS=resultGWAS$df
    ldBlockGWAS=resultGWAS$ldBlock

    #with the full region LD block in hand and a eQTL SNP data frame
    dfEQTL=df[!is.na(df$eQTL_P),]
    resulteQTL=intersectGwasAndLDBlock(dfEQTL, ldBlock, checkReversedAlleles=F)
    dfEQTL=resulteQTL$df
    ldBlockEQTL=resulteQTL$ldBlock

    writeECAVIARInputs(dfBoth, ldBlockBoth, dfGWAS, ldBlockGWAS, dfEQTL, ldBlockEQTL, outputBaseFileName, zScoreAbsValues=zScoreAbsValues)
    runECAVIAR(CAVIARPath, eCAVIARPath, outputBaseFileName, dryRun, verbose)
    df=summarizeResults(manifestFile, outputBaseFileName, df)
    outPDFFile=getPDFOutput(outputBaseFileName)
    plotResults(outputBaseFileName, gene_location_file_name, outPDFFile)

}

summarizeResults<-function (manifestFile, outputBaseFileName, df) {
    m=read.table(manifestFile, header=T, stringsAsFactors = F, sep="\t")

    #for each SNP block:
    #add markers for which SNPs are in the 95% confident list in GWAS, eQTL, BOTH.
    #add the co-location score.
    confSetEQTL=read.table(getConfSetEQTL(outputBaseFileName), header=F, stringsAsFactors = F)$V1
    confSetGWAS=read.table(getConfSetGWAS(outputBaseFileName), header=F, stringsAsFactors = F)$V1

    df$GWAS_95=!is.na(match(df$snp, confSetGWAS))
    df$EQTL_95=!is.na(match(df$snp, confSetEQTL))

    #95% confidence set for the co-location data can be found by:
    #1) sort data by Prob_in_pCausalSet
    #2) cumulative sum of Prob_in_pCausalSet
    #3) cutoff all cumsum > 0.05.
    caviarOut=read.table(getECaviarJointResult(outputBaseFileName), header=T, stringsAsFactors = F)
    caviarOut=caviarOut[order(caviarOut$Prob_in_pCausalSet, decreasing=F),]
    caviarOut$cumulativeP=cumsum(caviarOut$Prob_in_pCausalSet)
    idx95Conf=which(caviarOut$cumulativeP>=0.05)

    df$CLPP_95=!is.na(match(df$snp, caviarOut[idx95Conf,]$SNP_ID))
    #add the co-location score
    df$CLPP=caviarOut[match(df$snp, caviarOut$SNP_ID),]$CLPP

    outFile=getSummaryFileOneGene(outputBaseFileName)
    write.table(df, outFile, row.names=F, col.names=T, quote=F, sep="\t")
    return (df)
}

# inSNPFile="/downloads/GWAS_correlation/CAVIAR_FILES/1.Zscores.txt.gz"; outputBaseFileName="/downloads/GWAS_correlation/CAVIAR_FILES/out/1";outPDFFile=NULL

#' Plot the CAVIAR/eCAVIAR fine mapping results on a region.
#' @param outputBaseFileName The prefix for all files for this region, including the full directory path.
#' @param gene_location_file_name File containing the coordinates of all genes in the eQTL study.
#' @param outPDFFile If not null, create a new PDF file and plot to it.
plotResults<-function (outputBaseFileName, gene_location_file_name, outPDFFile=NULL) {

	a=read.table(getSummaryFileOneGene(outputBaseFileName), header=T, stringsAsFactors = F)
    geneName=unique(a$gene)

	if (!is.null(outPDFFile)) pdf(outPDFFile)

	# caviarOut=read.table(getECaviarJointResult(outputBaseFileName), header=T, stringsAsFactors = F)
	# idx=match(caviarOut$SNP_ID, a$snp)
	# caviarOut$pos=a[idx,]$pos
	# caviarOut=caviarOut[order(caviarOut$Prob_in_pCausalSet, decreasing=T),]

	#confSetEQTL=read.table(getConfSetEQTL(outputBaseFileName), header=F, stringsAsFactors = F)$V1
	#confSetGWAS=read.table(getConfSetGWAS(outputBaseFileName), header=F, stringsAsFactors = F)$V1

	confSetEQTL=a[a$EQTL_95==T,]$snp
	confSetGWAS=a[a$GWAS_95==T,]$snp

	#GWAS
	par(mfrow=c(3,1))
	par (mar=c(1,6,2,2))
	xlim=range(a$pos)
	gw=a[!is.na(a$gwas_P),]
	plot (gw$pos, -log10(gw$gwas_P), main=paste(geneName, "GWAS Fine Mapping [Caviar]"), ylab="GWAS pvalue [-log10]", xlim=xlim, axes=F, xlab="", cex=0.5)
	idxGwasConf=which(gw$GWAS_95==T)
	points(gw[idxGwasConf,]$pos, -log10(gw[idxGwasConf,]$gwas_P), col="red", pch=16, cex=1.25)
    axis(2)
	#eQTL
	eq=a[!is.na(a$eQTL_P),]
	plot (eq$pos, -log10(eq$eQTL_P), main=paste(geneName,"eQTL Fine Mapping [Caviar]"), ylab="eQTL pvalue [-log10]", xlim=xlim, axes=F, xlab="", cex=0.5)
	idxeQTLConf=which(eq$EQTL_95==T)
	points(eq[idxeQTLConf,]$pos, -log10(eq[idxeQTLConf,]$eQTL_P), col="red", pch=16, cex=1.25)
	axis(2)

	#GENES
	geneLoc=read.table(gene_location_file_name, header=T, stringsAsFactors = F)
	s=min (a$pos)
	e=max (a$pos)
	genes=geneLoc[geneLoc$chr==unique (a$chr) & geneLoc$s2>=s & geneLoc$s1<=e,]
	targetGene=genes[genes$geneid==geneName,]
	colnames(targetGene)[3:4]=c("start", "end")
	par (mar=c(4,6,0,2))
	plotGenes(genes, s, e, highlightFeatureRegion=targetGene, highlightColor="blue", plotXAxis=T)


	#LD plot?
    # ldBlock=read.table(getLDScoreMatrix(outputBaseFileName), header=F, stringsAsFactors=F)
    # ldBlock=abs(as.matrix(ldBlock))
    # ldBlock[upper.tri(ldBlock)] <- NA
    # image(ldBlock, frame=F,xaxt="n",yaxt="n")

    #library(LDheatmap)
    #LDheatmap(gdat=abs(as.matrix(ldBlock[1:100,1:100])), flip=T, add.map=F, newpage=F, title="")
    #LDheatmap(gdat=abs(as.matrix(ldBlock)), flip=T, add.map=F, newpage=F, title="")
    # library (gaston)
    # par(mfrow=c(2,1), mar=c(2,4,2,4))
    # plot (a$pos, -log10(a$gwas_P), main=paste(geneName, "GWAS Fine Mapping [Caviar]"), ylab="GWAS pvalue [-log10]")
    # idxGwasConf=match(confSetGWAS, a$snp)
    # points(a[idxGwasConf,]$pos, -log10(a[idxGwasConf,]$gwas_P), col="red", pch=16, cex=1.25)
    #
    # ldBlock=read.table(getLDScoreMatrix(outputBaseFileName), header=F, stringsAsFactors=F)
    # ldBlock=abs(as.matrix(ldBlock))
    # LD.plot(ldBlock, write.ld=NULL)


	#BOTH
	# confSetEQTLBoth=read.table(confSetBothEQTL(outputBaseFileName), header=F, stringsAsFactors = F)$V1
	# confSetGWASBoth=read.table(confSetBothGWAS(outputBaseFileName), header=F, stringsAsFactors = F)$V1
	#
	# par(mfrow=c(3,1), mar=c(2,4,2,4))
	# plot (a$pos, -log10(a$gwas_P), main="GWAS Fine Mapping with eQTL [GWAS 95% confidence SNPs]", ylab="GWAS pvalue [-log10]")
	# idxGwasConf=match(confSetGWASBoth, a$snp)
	# points(a[idxGwasConf,]$pos, -log10(a[idxGwasConf,]$gwas_P), col="red", pch=16, cex=1.25)
	#
	# plot (a$pos, -log10(a$eQTL_P), main="eQTL Fine Mapping with GWAS [eQTL 95% confidence SNPs]", ylab="eQTL pvalue [-log10]")
	# idxeQTLConf=match(confSetEQTLBoth, a$snp)
	# points(a[idxeQTLConf,]$pos, -log10(a[idxeQTLConf,]$eQTL_P), col="red", pch=16, cex=1.25)
	#
	# plot (caviarOut$pos, log10(caviarOut$CLPP), main="eCAVIAR Colocation [CLPP] [95% confidenece SNPs]", xlab="pos", ylab="Colocation pvalue [log10]")
	# idxBothCausalSet=which(caviarOut[,2]>=0.05)
	# points(caviarOut[idxBothCausalSet,]$pos, log10(caviarOut[idxBothCausalSet,]$CLPP), col="red", pch=16, cex=1.25)
	# pos=caviarOut[idxBothCausalSet,]$pos

	#plot just the CAVIAR 95% causal SNPs as red.
	#These are the INTERSECT.
	a=a[!is.na(a$gwas_P) & !is.na(a$eQTL_P) & !is.na(a$CLPP),]

	# idxBothCausalSet=which(caviarOut[,2]>=0.05)
	# pos=caviarOut[idxBothCausalSet,]$pos
	# idx95=match(pos, a$pos)

	par(mfrow=c(3,1), mar=c(1,6,2,2))

	#[95% confidenece SNPs]
	idx95=which(a$CLPP_95==T)
	plot (a$pos, -log10(a$gwas_P), main="GWAS Fine Mapping with eQTL [eCaviar]", ylab="GWAS pvalue [-log10]", xlab="", axes=F, cex=0.5)
	points(a[idx95,]$pos, -log10(a[idx95,]$gwas_P), col="red", pch=16, cex=1.25)
	axis(2)

	plot (a$pos, -log10(a$eQTL_P), main="eQTL Fine Mapping with GWAS [eCaviar]", axes=F, ylab="eQTL pvalue [-log10]", xlab="", cex=0.5)
	points(a[idx95,]$pos, -log10(a[idx95,]$eQTL_P), col="red", pch=16, cex=1.25)
	axis(2)
	maxCLPP=format (max (a$CLPP, na.rm=T), digits=4, scientific=T)
	par(mar=c(4,6,2,2))
	plot (a$pos, log10(a$CLPP), main=paste("CAVIAR Colocation: Max CLPP [", maxCLPP, "]", sep="") , ylab="CLPP [log10]", xlab="", axes=F, cex=0.5)
	idxBothCausalSet=which(a$CLPP_95)
	points(a[idxBothCausalSet,]$pos, log10(a[idxBothCausalSet,]$CLPP), col="red", pch=16, cex=1.25)
	pos=axTicks(1)
	#convert to 1kb.
	axis(1, at=pos, labels=pos/1e3)
	title(sub=paste("Chromosome", unique (a$chr), "(kb)"), line=2, cex.sub=1.5)
    axis(2)


	par(mfrow=c(2,1))
	par (mar=c(1,5,3,1))
	plot (a$pos, log10(a$CLPP), main=paste("CAVIAR Colocation: Max CLPP [", maxCLPP, "]", sep="") , ylab="CLPP [log10]", xlab="", axes=F, cex=0.5)
	axis(2)
	idxBothCausalSet=which(a$CLPP_95)
	points(a[idxBothCausalSet,]$pos, log10(a[idxBothCausalSet,]$CLPP), col="red", pch=16, cex=1.25)
	par (mar=c(4,5,0,1))
	plotGenes(genes, s, e, highlightFeatureRegion=targetGene, highlightColor="blue", plotXAxis=T)


	if (!is.null(outPDFFile)) dev.off()
}

#' Run CAVIAR/eCAVIAR on a single region as defined by the outputBaseFileName file set.
#' Generates outputs to the outputBaseFileName with appropriate suffixes.
#' @param outputBaseFileName The prefix for all files for this region, including the full directory path.
#' @inheritParams runPipeline
runECAVIAR<-function (CAVIARPath, eCAVIARPath, outputBaseFileName, dryRun=F, verbose=F) {

	runOne<-function (command, verbose, dryRun) {
		if (verbose) cat (command, "\n")
		if (!dryRun) system(command, ignore.stdout=T, ignore.stderr = T, wait=T)
	}

	#GWAS ONLY
	#/broad/mccarroll/nemesh/caviar/CAVIAR-C++/CAVIAR -o 470.GWAS -l 470.LD.txt -z 470.gwas.Z.txt
	gwasCall=paste(CAVIARPath, " -o ", paste(outputBaseFileName, ".GWAS -l ", sep=""), getLDScoreMatrixGWAS(outputBaseFileName), " -z ", getGWASFile(outputBaseFileName), sep="")

	runOne(gwasCall, verbose, dryRun)

	#EQTL ONLY
	#/broad/mccarroll/nemesh/caviar/CAVIAR-C++/CAVIAR -o 470.eQTL -l 470.LD.txt -z 470.eQTL.Z.txt
	eQTLCall=paste(CAVIARPath, " -o ", paste(outputBaseFileName, ".eQTL -l ", sep=""), getLDScoreMatrixeQTL(outputBaseFileName), " -z ", geteQTLFile(outputBaseFileName), sep="")
	runOne(eQTLCall, verbose, dryRun)

	#BOTH

	#/broad/mccarroll/nemesh/caviar/CAVIAR-C++/eCAVIAR -o 470.both -l 470.LD.txt -z 470.eQTL.Z.txt -z 470.gwas.Z.txt -l 470.LD.txt
	bothCall=paste(eCAVIARPath, " -o ", paste(outputBaseFileName, ".both -l ", sep=""), getLDScoreMatrix(outputBaseFileName), " -z ", getBotheQTLFile(outputBaseFileName),
				   " -l ", getLDScoreMatrix(outputBaseFileName), " -z ", getBothGWASFile(outputBaseFileName), sep="")
	runOne(bothCall, verbose, dryRun)

}


#' Write out the inputs to CAVIAR/eCAVIAR to perform fine mapping on GWAS/EQTL data sets.
#' Outputs have no headers.
#' GWAS/eQTL has 2 columns: snp id and Z score.  Tab separated.
#' LD block has no headers, square matrix of LD scores.  Space separated.
#' @param dfBoth A data frame containing the GWAS/eQTL SNP data for a region
#' @param ldBlockBoth The LD matrix for the intersect of the GWAS / eQTL region, containing the same SNPs as the dfBoth parameter, in the same order.
#' @param dfGWAS A data frame containing the GWAS SNP data for a region
#' @param ldBlockGWAS The LD matrix for the intersect of the GWAS region, containing the same SNPs as the dfGWAS parameter, in the same order.
#' @param dfEQTL A data frame containing the eQTL SNP data for a region.
#' @param ldBlockEQTL The LD matrix for the intersect of the eQTL region, containing the same SNPs as the dfEQTL parameter, in the same order.
#' @param outputBaseFileName The prefix for all files for this region, including the full directory path.
#' @param zScoreAbsValues If true, use the absolute values of the Z-scores in the LD matrix and gwas/eQTL Zscores.
writeECAVIARInputs <-function (dfBoth, ldBlockBoth, dfGWAS, ldBlockGWAS, dfEQTL, ldBlockEQTL, outputBaseFileName, zScoreAbsValues=F) {
    outputGWAS=getGWASFile(outputBaseFileName)
    outputeQTL=geteQTLFile(outputBaseFileName)
    outputBothGWAS=getBothGWASFile(outputBaseFileName)
    outputBothEQTL=getBotheQTLFile(outputBaseFileName)

    outputLDScoreMatrixBoth=getLDScoreMatrix(outputBaseFileName)
    outputLDScoreMatrixGWAS=getLDScoreMatrixGWAS(outputBaseFileName)
    outputLDScoreMatrixEQTL=getLDScoreMatrixeQTL(outputBaseFileName)

    if (zScoreAbsValues) {
        ldBlockBoth=abs(ldBlockBoth)
        ldBlockGWAS=abs(ldBlockGWAS)
        ldBlockEQTL=abs(ldBlockEQTL)

        dfBoth$gwas_Z=abs(dfBoth$gwas_Z)
        dfBoth$eQTL_Z=abs(dfBoth$eQTL_Z)
        dfGWAS$gwas_Z=abs(dfGWAS$gwas_Z)
        dfGWAS$eQTL_Z=abs(dfGWAS$eQTL_Z)
        dfEQTL$gwas_Z=abs(dfEQTL$gwas_Z)
        dfEQTL$eQTL_Z=abs(dfEQTL$eQTL_Z)

    }

    write.table(dfGWAS[,c("snp", "gwas_Z")], outputGWAS, row.names=F, col.names = F, quote=F, sep="\t")
    write.table(dfEQTL[,c("snp", "eQTL_Z")], outputeQTL, row.names=F, col.names = F, quote=F, sep="\t")

    write.table(dfBoth[,c("snp", "gwas_Z")], outputBothGWAS, row.names=F, col.names = F, quote=F, sep="\t")
    write.table(dfBoth[,c("snp", "eQTL_Z")], outputBothEQTL, row.names=F, col.names = F, quote=F, sep="\t")

    write.table(ldBlockBoth, outputLDScoreMatrixBoth, row.names=F, col.names=F, quote=F, sep=" ")
    write.table(ldBlockGWAS, outputLDScoreMatrixGWAS, row.names=F, col.names=F, quote=F, sep=" ")
    write.table(ldBlockEQTL, outputLDScoreMatrixEQTL, row.names=F, col.names=F, quote=F, sep=" ")
}

#' Builds LD blocks for the snps in the input data frame.
#'
#' Because LDBlocks can only be querried out of a VCF by region, we need to select out the intersect of
#' variants that are in both the GWAS/eQTL data as well as the LD matrix.  The data is the ordered so
#' all data sets have their SNPs in the same order for CAVIAR/eCAVIAR.
#'
#'
#' @param df A data frame containing the GWAS/eQTL combined data.
#' @param ldBlock A matrix of LD generated by bcftools for a region
#' @param checkReversedAlleles check the opposite allele orientation in the LD block.  When found, reverse the Z scores in the data frame.
#' @return a list containing the GWAS/eQTL df and ldblock for the intersect of SNPs in both.
intersectGwasAndLDBlock<-function (df, ldBlock, checkReversedAlleles=F) {
    #TODO: if SNP name doesn't match, check reversed allele?
    getReverseName<-function (x) {
        xx=strsplit(x, ":", fixed=T)[[1]]
        return (paste(xx[1], xx[2], xx[4], xx[3], sep=":"))
    }
    getReverseName=Vectorize(getReverseName)

    if (checkReversedAlleles) {
        snp_flip=getReverseName(df$snp)
        snpsBothFlipped=intersect (snp_flip, colnames(ldBlock))
        idx=match(snpsBothFlipped, snp_flip)
        #flip the Z score of the snps that were flipped.
        df[idx,]$gwas_Z=df[idx,]$gwas_Z*-1
        df[idx,]$eQTL_Z=df[idx,]$eQTL_Z*-1
        #replace their snp name.
        df[idx,]$snp=snp_flip[idx]
    }

    snpsBoth=intersect (df$snp, colnames(ldBlock))

    #remaining: setdiff(c(snpsBoth, snpsBothFlipped), colnames(ldBlock))
    if (length(snpsBoth)==0) {
        stop(paste("No SNPs in common between LD block and GWAS/eQTL SNPs for gene", unique (df$gene)))
    }
    #the remaining GWAS/eQTL data:
    idx=match(snpsBoth, df$snp)
    df2=df[idx,]

    #reorder GWAS/eQTL by position
    df2=df2[order(df2$pos, decreasing = F),]

    #the remaining LD block data:
    idx=match(df2$snp, colnames(ldBlock))
    if (length(which(is.na(idx)))>0)
        stop("SNPs requested from LD block that don't exist.  This shouldn't happen!")

    ldBlockFinal=ldBlock[idx,idx]

    result=list(df=df2, ldBlock=ldBlockFinal)
    return (result)
}

#' Generate an output file based on the input file name and a specified output directory.
#'
#' @param inSNPFile the input file name.  Used to generate a prefix name for all output files
#' @param outDir The output directory where files will be generated.
#'
#' @return The output file name prefix, containing the full path + basename of the output file.
getOutputFileNameBase<-function (inSNPFile, outDir) {
    iteration=strsplit (basename(inSNPFile), ".", fixed=T)[[1]][1]
    outBase=paste(outDir, "/", iteration, sep="")
    return (outBase)
}

#
#' Get the chr:start-end string from the eQTL/GWAS data frame.
#'
#' @param df A data frame containing the GWAS/eQTL combined data.
#'
#' @return The region for the data frame.  Always a single chromosome, with a start and end position.
getFormattedGenomicRegion<-function (df) {
    chr=unique(df$chr)
    regionStr=paste(chr, ":", min(df$pos), "-", max(df$pos), sep="")
    return (regionStr)
}

#buil a bcftools call to get the LD matrix, read in the LD matrix, and return.
#' Query a VCF to extract a square matrix of geneotype correlation.
#'
#' @param df A data frame containing the GWAS/eQTL combined data.
#' @param outputBaseFileName The prefix for all files for this region, including the full directory path.
#' @param cleanUpFiles delete the intermediate files generated if true (saves space!)
#' @param dryRun If true, do not run bcftools call, only echo command string if verbose=T.
#' @param verbose Print the command string.
#'
#' @return A square matrix of genotype correlation.
#' @inheritParams runPipeline
buildLDBlock<-function (df, vcfFile, outputBaseFileName, bcfToolsLocation, ldScoreBcftoolsPluginLocation, cleanUpFiles=F, dryRun, verbose=T) {
    outFile=getLDFullRegion(outputBaseFileName)
    regionStr=getFormattedGenomicRegion(df)
	# bcftools filter -Ou -i 'AF>.01 && AF<.99'
    # /broad/mccarroll/nemesh/1kg_VCF/EUR/EUR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.ref_normalized.vcf.gz
    # -r 15:78608329-79136713 | bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT'  |
    # bcftools +/home/unix/nemesh/bin/LDmatrix.so -- -p | gzip > bcftools.ldMatrix.txt.gz

    bcfToolsCall=paste(bcfToolsLocation, " filter -Ou -i 'AF>.01 && AF<.99' ", vcfFile, " -r ", regionStr, "| ",
    				   bcfToolsLocation, " annotate -x ID -I +'%CHROM:%POS:%REF:%ALT'| ",
    				   bcfToolsLocation, " +", ldScoreBcftoolsPluginLocation, " -- -p | gzip > ", outFile, sep="")

    if (verbose) cat (bcfToolsCall, "\n")
    if (dryRun) {
    	return (NULL)
	}
    else {
    	system(bcfToolsCall)
    	ldMatrix=read.table(outFile, header=T, stringsAsFactors = F, check.names = F)
    	if (cleanUpFiles) file.remove(outFile)
    	return (ldMatrix)
    }
}

#' Read the eQTL file.  Restrict to the best gene for each SNP if required
#' @param eQTLFile The eQTL file to read in.
#' @param getBestGenePerSNP Should the set of SNPs be filtered to the best gene per SNP (each SNP appears once),
#' or should a SNP be allowed to appear for each gene it has been tested in?
#' @return The dataframe of eQTL results.  Chromosome/position and ref/alt alleles are parsed out of the SNP name.  Chromosome X is set to chromosome 23.
getEQTLResults<-function (eQTLFile, getBestGenePerSNP=T) {
    b=DropSeq.utilities::fastRead(eQTLFile)
    #x=b[b$SNP=="1:1245368:G:A",]
    bestGene<-function (x) {
        idx=which.min(x$`p-value`)
        x[idx,]
    }

    if (getBestGenePerSNP) {
        bb=b[,bestGene(.SD), by="SNP"]
    } else {
        bb=b
    }
    #parse out the chromosome/position, ref/alt alleles.
    z=strsplit (bb$SNP, ":", fixed=T)
    bb$chr=sapply(z, function (x) x[1])
    bb$pos=as.numeric(sapply(z, function (x) x[2]))
    bb$a1=sapply(z, function (x) x[3])
    bb$a2=sapply(z, function (x) x[4])
    #fix chromosome X
    idxX=which(bb$chr=="X")
    if (length(idxX)>0) bb[idxX,]$chr<-"23"
    return (bb)
}


#' Parse in GWAS results.
#'
#' If the chromosome name contains "chr", that is removed.
#'
#' @param gwasFile A GWAS file to parse.
#' @return A dataframe containing the GWAS results.
getGWasResults<-function (gwasFile) {
    #TODO: validate input doesn't have duplicate records SNP chr/pos/A1/A2 entries.
    a=DropSeq.utilities::fastRead(gwasFile)
    key=paste(a$CHR, a$POS, a$A1, a$A2, sep=":")
    dups=key[which(duplicated(key))]
    if (length(dups)>0) {
        stop(paste ("Duplicate SNPs found in GWAS data: ", paste(head (dups, 5), collapse=",")))
    }
    #colnames(a)[1]="chr"
    #colnames (a)[5]="pos"
    a$CHR=sub ("chr", "", a$CHR)
    return (a)
}

#' Intersect the GWAS and eQTL effects.
#' Find SNPs common bot both data sets genome wide.   Check the reference and alternate
#' alleles for each data set, and make them match in both data sets, flipping the odds ratio
#' of the GWAS data as appropriate so the eQTL and GWAS use the same allele to determine direction
#' of effect.
#' @param gwas The GWAS data frame
#' @param eQTL The eQTL data frame
#'
#' @return A list containing the intersected/modified GWAS and eQTL data frames.
intersectGwasEQTL<-function (gwas, eQTL) {
    gwas$key=paste(gwas$CHR, gwas$POS, sep=":")
    eQTL$key=paste(eQTL$chr, eQTL$pos, sep=":")
    both=intersect(gwas$key, eQTL$key)
    #gwasOnly=setdiff(gwasKey, eQTLKey)

    #gwas SNPs only show up once.
    idxGwas=match(both, gwas$key)
    gwas=gwas[idxGwas,]
    #make distinct eQTL data (each SNP only shows up once) to flip gwas results.
    idxEQTL=match(both, eQTL$key)
    eQTLDistinct=eQTL[idxEQTL,]

    #both should have the same SNPs in the same order.
    idxAllelesMatch=which(gwas$A1==eQTLDistinct$a1 & gwas$A2==eQTLDistinct$a2)
    idxAllelesFlip=which(gwas$A1==eQTLDistinct$a2 & gwas$A2==eQTLDistinct$a1)
    #for alleles that were flipped but match, flip the odds ratio.
    gwas=flipGwasAlleles(gwas, idxAllelesFlip)

    #indels are coded by I?D in the GWA results, and by the actual alleles in the eQTL results...
    #fix those indels in the gwas results.
    gwas=fixIndels(gwas, eQTLDistinct)

    #note in eQTL data the same SNP can show up in multiple genes.
    idxEQTL=which(!is.na(match(eQTL$key, both)))
    eQTL=eQTL[idxEQTL,]


    return (list(gwas=gwas, eQTL=eQTL))

}

#' Union the GWAS and eQTL effects.
#' Find SNPs common bot both data sets genome wide.   Check the reference and alternate
#' alleles for each data set, and make them match in both data sets, flipping the odds ratio
#' of the GWAS data as appropriate so the eQTL and GWAS use the same allele to determine direction
#' of effect.
#' For SNPs present in only 1 data set, retain those "as is".
#' @param gwas The GWAS data frame
#' @param eQTL The eQTL data frame
#'
#' @return A list containing the unioned/modified GWAS and eQTL data frames.
unionGwasEQTL<-function (gwas, eQTL) {
    intersected=intersectGwasEQTL(gwas,eQTL)
    g=intersected$gwas
    q=intersected$eQTL

    #which gwas not in intersected gwas?
    gwas$key=paste(gwas$CHR, gwas$POS, sep=":")
    #eQTL$key=paste(eQTL$chr, eQTL$pos, sep=":")

    #find the gwas SNPs that were not in the intersect and add them to the result.
    #SNPs are only distinct when chr, pos, A1, A2 match.
    #this is confusing as the returned GWAS intersect can have the alleles flipped.
    gwasKeyI=paste(g$CHR, g$POS, g$A1, g$A2, sep=":")
    #if the fixed key is found in either keyset, that SNP is done.
    gwasKeyAll1=paste(gwas$CHR, gwas$POS, gwas$A1, gwas$A2, sep=":")
    gwasKeyAll2=paste(gwas$CHR, gwas$POS, gwas$A2, gwas$A1, sep=":")

    #check if the gwas result is in the intersect with both allele orientations.
    idx=which(!is.na(match(gwasKeyAll1, gwasKeyI)))
    idx1=which(!is.na(match(gwasKeyAll2, gwasKeyI)))

    #what's in the intersect but not in the original results?  Should be 0 entries.
    both=c(gwasKeyAll1[idx],gwasKeyAll2[idx1])
    if (length(both)!=length(unique(both))) {
        stop(paste("Some SNPs occur multiple times in GWAS data.", paste(both[which(duplicated(both))], collapse=", "), sep=" "))
    }
    missing=setdiff(gwasKeyI, both)
    if (length(missing)>0) {
        warning ("Some confusion/bug with intersection with eQTL results and full GWAS result")
    }

    #either one is a hit.  All gwas entries from the intersected results should be found.
    #idxBoth is the set to NOT include.
    idxBoth=union(idx,idx1)

    gwasR=rbind(g, gwas[-idxBoth,])
    if (dim (gwas)[1]!=dim(gwasR)[1]) {
        warning("Number of gwas SNPs changed.")
    }

    #Include gene in the key. This makes it unique if the SNP is used in multiple genes.
    #this one is easier than the GWAS check as alleles don't flip.
    key2All=paste(eQTL$chr, eQTL$pos, eQTL$a1, eQTL$a2, eQTL$gene, sep=":")
    key2I=paste(q$chr, q$pos, q$a1, q$a2, q$gene, sep=":")
    #find eQTL SNPs not in the intersect.
    idxQ=match(setdiff(key2All, key2I), key2All)
    #get rid of intersected eQTL key column so columns are consistent for rbind.
    q[,key:=NULL]
    eQTLR=rbind(q, eQTL[idxQ,])

    if (dim (eQTL)[1]!=dim(eQTLR)[1]) {
        warning("Number of eQTL SNPs changed.")
    }
    #make the final key
    gwasR$key=paste(gwasR$CHR, gwasR$POS, gwasR$A1, gwasR$A2, sep=":")
    eQTLR$key=paste(eQTLR$chr, eQTLR$pos, eQTLR$a1, eQTLR$a2, sep=":")
    return (list(gwas=gwasR, eQTL=eQTLR))

}


#' Exchange the reference and alt as well as the direction of the odds ratio at the given index.
#'
#' Flip the GWAS results alleles at the requested index.
#' Ths also flips the odds ratio (transform by 1/OR).
#' @param gwas The GWAS dataframe
#' @param idxAllelesFlip A vector of positions at which the alleles should be flipped.
#' @return The GWAS data frame with the appropriate alleles/odds ratios flipped.
flipGwasAlleles<-function (gwas,idxAllelesFlip ) {
    if (length(idxAllelesFlip)>0) {
        #flip the odds ratio if it exists
        if ("OR" %in% colnames(gwas)) {
            gwas[idxAllelesFlip,]$OR=1/gwas[idxAllelesFlip,]$OR
        }
        #flip the Z score.
        gwas[idxAllelesFlip]$Z=gwas[idxAllelesFlip]$Z*-1
        a2=gwas[idxAllelesFlip,]$A2
        gwas[idxAllelesFlip,]$A2=gwas[idxAllelesFlip,]$A1
        gwas[idxAllelesFlip,]$A1=a2
    }
    return (gwas)
}


#' Find SNPs are indels in the gwas and eQTL data set and give them the same ref allele / direction of effect.
#'
#' Finds SNPs that have D or I as one of their alleles in the GWAS data set and finds if the eQTL has a similar deletion
#' or insertion by looking at the number of bases in the allele.  Sets the reference/alt allele to match in the two data sets,
#' and flips alleles/odds ratios in the GWAS result as neccesary.
#' @param gwas The GWAS data frame
#' @param eQTL The eQTL data frame
#'
#' @return The GWAS data frame, modified so that the alleles and direction of effect match the eQTL data
fixIndels<-function (gwas, eQTL) {
    #find the indels and slim both data sets down to them.
    ncharAllele1=nchar(eQTL$a1)
    ncharAllele2=nchar(eQTL$a2)
    #gwas is coded as an indel, and one of the eQTL alleles is longer than 1.
    idx=which((gwas$a1=="D" | gwas$a2=="D") & (ncharAllele1>1 | ncharAllele2>1))

    gwas2=gwas[idx,]
    eQTL2=eQTL[idx,]

    ncharAllele1=nchar(eQTL2$a1)
    ncharAllele2=nchar(eQTL2$a2)
    #to fix the indels, match the multi-character allele to the I, since the D is the single character.

    #matching alleles
    idxAllelesMatch=which((gwas2$a1=="D" & ncharAllele1==1) | (gwas2$a2=="D" & ncharAllele2==1))
    idxAllelesFlip=which((gwas2$a1=="D" & ncharAllele2==1) | (gwas2$a2=="D" & ncharAllele1==1))
    idxAllelesDisagree=setdiff(1:dim(gwas2)[1], c(idxAllelesMatch, idxAllelesFlip))
    gwas2=flipGwasAlleles(gwas2, idxAllelesFlip)

    gwas[idx,]=gwas2
    return (gwas)

}

getZ<-function (pval) -1* qnorm(pval/2,F)

#' Find SNPs in both the GWAS and eQTL data sets and intersect
#'
#' For the GWAS and eQTL data, find SNPs that are common to both sets.  Calculate the Z-score for the eQTL
#' data set based on the sign of the beta (>0 = +).
#' @param geneName The gene to operate on
#' @param gwas The GWAS data frame
#' @param eQTL The eQTL data frame
#' @param verbose If true, write out the name of each gene as it is processed.
#'
#' @return A data frame for a gene containing the intersect of SNPs in the eQTL and GWAS data, with Z-scores calculated.
getGwasEQTLOverlappingSNPs<-function (geneName, gwas, eQTL, verbose=F) {
    eQ=eQTL[eQTL$gene==geneName,]
    idx=match(eQ$key, gwas$key)
    if (length(idx)==0) return (NULL)
    if (length(which(is.na(idx)))>0)
        warning("SNPs in eQTL data set not found in GWAS.")
    gwasP=gwas[idx,]

    #eQTL Z scores
    eQTLZScore=sapply(eQ$'p-value', getZ)
    #modify sign by beta
    idxNeg=which(eQ$beta<0)
    #idxNeg=which(eQ$beta<1)
    eQTLZScore[idxNeg]=eQTLZScore[idxNeg]*-1

    if ("OR" %in% colnames(gwasP)) {
        df=data.frame(snp=eQ$SNP, snpid=gwasP$SNP, gene=eQ$gene, chr=unique (eQ$chr), pos=gwasP$POS,
                      gwas_P=gwasP$P,  gwas_or=gwasP$OR, gwas_Z=gwasP$Z,
                      eQTL_P=eQ$`p-value`, eQTL_beta=eQ$beta, eQTL_Z=eQTLZScore, stringsAsFactors = F)
    } else {
        df=data.frame(snp=eQ$SNP, snpid=gwasP$SNP, gene=eQ$gene, chr=unique (eQ$chr), pos=gwasP$POS,
                      gwas_P=gwasP$P, gwas_Z=gwasP$Z,
                      eQTL_P=eQ$`p-value`, eQTL_beta=eQ$beta, eQTL_Z=eQTLZScore, stringsAsFactors = F)
    }
    #plot (df$gwas_Z, df$eQTL_Z, xlab="GWAS Z-score", ylab="eQTL Z-score", main=unique (df$gene))
    #plot (-log10(df$gwas_P), -log10(df$eQTL_P), xlab="GWAS pval [-log10]", ylab="eQTL pval [-log10]", main=unique (df$gene))
    if (verbose) cat (geneName, "\n")
    return (df)
}

#system.time (getGwasEQTLUnionSNPs(geneName, gwas, eQTL, verbose=F))
getGwasEQTLUnionSNPs<-function (geneName, gwas, eQTL, verbose=F) {
    #prep the eQTL data.
    eQ=eQTL[eQTL$gene==geneName,]
    eQTLZScore=sapply(eQ$'p-value', getZ)
    #modify sign by beta
    idxNeg=which(eQ$beta<0)
    #idxNeg=which(eQ$beta<1)
    eQTLZScore[idxNeg]=eQTLZScore[idxNeg]*-1
    eQ$z=eQTLZScore

    #get the GWAS data in the same region
    #this two-step approach is far faster than the single step approach.
    r=range(as.numeric(eQ$pos))
    gwasP=gwas[gwas$CHR==unique(eQ$chr),]
    gwasP=gwasP[gwasP$POS>=r[1] & gwasP$POS<=r[2],]
    #gwasP=gwas[gwas$CHR==unique(eQ$chr) & gwas$POS>=min(eQ$pos) & gwas$POS<=max(eQ$pos),]

    keyBoth=intersect(eQ$key, gwasP$key)
    keyEQTLOnly=setdiff(eQ$key, keyBoth)
    keyGWASOnly=setdiff(gwasP$key, keyBoth)

    #assemble where both data sets have results.

    beQ=eQ[match(keyBoth, eQ$key),]
    bGwasP=gwasP[match(keyBoth, gwasP$key),]
    bGwasP$gene=geneName

    df=data.frame(snp=beQ$SNP, snpid=bGwasP$SNP, gene=beQ$gene, chr=unique (beQ$chr), pos=bGwasP$POS,
                  gwas_P=bGwasP$P, gwas_Z=bGwasP$Z,
                  eQTL_P=beQ$`p-value`, eQTL_beta=beQ$beta, eQTL_Z=beQ$z, stringsAsFactors = F)
    #for the other results, cache a few uniques


    #eQTL only
    beQ=eQ[match(keyEQTLOnly, eQ$key),]
    df2=data.frame(snp=beQ$SNP, snpid=rep(NA, dim(beQ)[1]), gene=beQ$gene, chr=unique (beQ$chr), pos=beQ$pos,
                  gwas_P=rep(NA, dim(beQ)[1]), gwas_Z=rep(NA, dim(beQ)[1]),
                  eQTL_P=beQ$`p-value`, eQTL_beta=beQ$beta, eQTL_Z=beQ$z, stringsAsFactors = F)
    #GWAS only
    bGwasP=gwasP[match(keyGWASOnly, gwasP$key),]
    bGwasP$gene=geneName
    df3=data.frame(snp=bGwasP$key, snpid=bGwasP$SNP, gene=bGwasP$gene, chr=bGwasP$CHR, pos=bGwasP$POS,
                  gwas_P=bGwasP$P, gwas_Z=bGwasP$Z,
                  eQTL_P=rep(NA, dim(bGwasP)[1]), eQTL_beta=rep(NA, dim(bGwasP)[1]), eQTL_Z=rep(NA, dim(bGwasP)[1]), stringsAsFactors = F)

    dfF=rbind (df, df2, df3)
    dfF=dfF[order(dfF$pos),]
}


#inDir="/downloads/GWAS_correlation/CAVIAR_FILES/snpBlocks_daner_NGN2/out"
#manifestFile="/downloads/GWAS_correlation/CAVIAR_FILES/snpBlocks_daner_NGN2/gene_manifest.txt"
#outAllSNPsFile="/downloads/GWAS_correlation/CAVIAR_FILES/snpBlocks_daner_NGN2/ngn2_all_snps.txt"
aggregateSummaryResults<-function (inDir, manifestFile, outAllSNPsFile) {
    m=read.table(manifestFile, header=T, stringsAsFactors=F)
    f=list.files(inDir, pattern=".summary", full.names = T)
    #file=f[1]
    readOne<-function (file, m) {
        x=read.table(file, header=T, stringsAsFactors = F)
        return (x)
    }
    r=lapply(f, readOne, m)
    r=do.call(rbind, r)
    write.table(r, outAllSNPsFile, row.names=F, col.names = T, quote=F, sep="\t")

    #get the max CLPP score per gene.
    r=setDT(r)
    #prevent no visible binding for global variable messages.
    CLPP=gene=NULL;
    maxCLPPPerGene=r[,max(CLPP, na.rm=T), by=gene]
    idx=match(m$gene, maxCLPPPerGene$gene)
    m$MAX_CLPP=maxCLPPPerGene[idx,]$V1

    #interesting genes?
    idx=which(m$MAX_CLPP>=0.01 & m$gwas_P<1e-6 & m$qvalue<0.05)
    m[idx,]

}


# buildLDBlockPlink<-function (df, vcfFile, outputBaseFileName, ldScoreBcftoolsPluginLocation, cleanUpFiles=F) {
# 	#outFileLD=paste(outputBaseFileName, ".plink.ld.gz", sep="")
# 	#outFileSNPs=paste(outputBaseFileName, ".plink.bim", sep="")
# 	outFileLD="/downloads/GWAS_correlation/CAVIAR_FILES/out/plink.keep.ld.gz"
# 	outFileSNPs="/downloads/GWAS_correlation/CAVIAR_FILES/out/plink.keep.bim"
# 	regionStr=getFormattedGenomicRegion(df)
# 	#bcftools filter -Ou -i 'AF>.01 && AF<.99' /broad/mccarroll/nemesh/1kg_VCF/EUR/EUR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -r 15:78608329-79136713 |
# 	#bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT'|
# 	#bcftools +$HOME/bin/LDmatrix.so -- | gzip > test_maf_filter.txt.gz
#
# 	#bcfToolsCall=paste("bcftools filter -Ou -i 'AF>.01 && AF<.99' ", vcfFile, " -r ", regionStr, "| ",
# 	#				   "bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT'| ",
# 	#				   "bcftools +", ldScoreBcftoolsPluginLocation, " -- | gzip > ", outFile, sep="")
#
# 	#system(bcfToolsCall)
#
# 	ldMatrix=read.table(outFileLD, header=F, stringsAsFactors = F, check.names = F)
# 	snps=read.table(outFileSNPs, header=F, stringsAsFactors = F, check.names = F)
# 	colnames(ldMatrix)=snps$V2
# 	if (cleanUpFiles) file.remove(outFile)
# 	return (ldMatrix)
# }

#I'm not quite sure I need this...
# filterStrandFlipSNPs<-function (df) {
# 	s=strsplit (df$snp, ":")
# 	a1=sapply(s, function (x) x[3])
# 	a2=sapply(s, function (x) x[4])
#
# 	idx=which((a1=="A" & a2=="T") | (a2=="A" & a1=="T") | (a1=="C" & a2=="G") | (a2=="C" & a1=="G"))
# 	df=df[-idx,]
# 	return (df)
# }
