# add in a paramter that makes the report emit the effect size of each gene.
# add in a comparison that adds a column when the sign of all experiments agrees.
#
# strip the alleles off the SNP names?

# library (data.table)

# qvalueThreshold=0.05; filterChromosomeList=NULL
# geneDescriptionsFile=c("/downloads/cell_selection/eQTL/hg19.gene_descriptions.txt")
#
# experimentNames=c("coding", "coding_intronic")
# unfilteredQTLFiles=c("/downloads/cell_selection/eQTL/sign_test/d21Ngn2plusGliaLive.whole_cells.maf_0.19_cisDist_10kb.eQTL_results.txt",
#                      "/downloads/cell_selection/eQTL/sign_test/d21Ngn2plusGliaLive.exonic+intronic.nuclei.maf_0.19_cisDist_10kb.eQTL_results.txt")
#
# permutedQTLFiles=c("/downloads/cell_selection/eQTL/d21Ngn2plusGliaLive.whole_cells.maf_0.20_cisDist_10kb.eQTL_permuted_results.txt",
#                    "/downloads/cell_selection/eQTL/d21Ngn2plusGliaLive.exonic+intronic.nuclei.maf_0.20_cisDist_10kb.eQTL_permuted_results.txt")
# outGeneLevelSummaryFile=("/downloads/cell_selection/eQTL/d21NGN2_coding_vs_intronic_aggregate_results.genes.txt")
# outExperimentLevelSummaryFile=("/downloads/cell_selection/eQTL/d21NGN2_wholecell_vs_nuclei_aggregate_results.summary.txt")


# experimentNames=c("VillageC", "Heat1-6")
# unfilteredQTLFiles=c("/downloads/eQTL/aggregate/data/DropulationVillageC_all.unfiltered.eQTL_results_10kb.txt",
#                      "/downloads/eQTL/aggregate/data/heat1-6.unfiltered.eQTL_results_10kb.txt")
#
# permutedQTLFiles=c("/downloads/eQTL/aggregate/data/DropulationVillageC_all.eQTL_permuted_results.txt",
#                    "/downloads/eQTL/aggregate/data/heat1-6.eQTL_permuted_results.txt")
# outGeneLevelSummaryFile=("/downloads/eQTL/aggregate/eQTL_ESC_aggregate_results.genes.txt")
# outExperimentLevelSummaryFile=("/downloads/eQTL/aggregate/eQTL_ESC_aggregate_results.summary.txt")
#
# experimentNames=c("rNPC_Pilot", "rNPC_Alive")
# unfilteredQTLFiles=c("/downloads/eQTL/aggregate/data/rNPC_pilot.unfiltered.eQTL_results_10kb.txt",
#                      "/downloads/eQTL/aggregate/data/rNPC_Alive_AllPools.auto.unfiltered.eQTL_results_10kb.txt")
#
# permutedQTLFiles=c("/downloads/eQTL/aggregate/data/rNPC_pilot.eQTL_permuted_results.txt",
#                    "/downloads/eQTL/aggregate/data/rNPC_Alive_AllPools.auto.eQTL_permuted_results.txt")
# outGeneLevelSummaryFile=("/downloads/eQTL/aggregate/eQTL_rNPC_aggregate_results.genes.txt")
# outExperimentLevelSummaryFile=("/downloads/eQTL/aggregate/eQTL_rNPC_aggregate_results.summary.txt")

# experimentNames=c("VillageC", "rNPC_Alive")
# unfilteredQTLFiles=c("/downloads/eQTL/aggregate/data/DropulationVillageC_all.unfiltered.eQTL_results_10kb.txt",
#                      "/downloads/eQTL/aggregate/data/rNPC_Alive_AllPools.auto.unfiltered.eQTL_results_10kb.txt")
#
# permutedQTLFiles=c("/downloads/eQTL/aggregate/data/DropulationVillageC_all.eQTL_permuted_results.txt",
#                    "/downloads/eQTL/aggregate/data/rNPC_Alive_AllPools.auto.eQTL_permuted_results.txt")
# outGeneLevelSummaryFile=("/downloads/eQTL/aggregate/tissue_specific_eQTL_aggregate_results.genes.txt")
# outExperimentLevelSummaryFile=("/downloads/eQTL/aggregate/tissue_specific_eQTL_aggregate_results.summary.txt")

# experimentNames=c("VillageC", "rNPC_Alive")
# unfilteredQTLFiles=c("/downloads/Gender_Expression/DropulationVillageC_all.unfiltered.gender_biased_genes.txt",
#                      "/downloads/Gender_Expression/rNPC_Alive_AllPools.auto.unfiltered.gender_biased_genes.txt")
#
# permutedQTLFiles=c("/downloads/Gender_Expression/DropulationVillageC_all.gender_biased_genes.txt",
#                    "/downloads/Gender_Expression/rNPC_Alive_AllPools.auto.gender_biased_genes.txt")
# outGeneLevelSummaryFile=("/downloads/Gender_Expression/tissue_specific_eQTL_aggregate_results.genes.txt")
# outExperimentLevelSummaryFile=("/downloads/Gender_Expression/tissue_specific_eQTL_aggregate_results.summary.txt")
# expNames=c("VillageC", "rNPC")
# filterChromosomeList=NULL

#' Summarize results across many eQTL experiments.
#'
#' Captures the set of genes with a qvalue below the threshold in any experiment, and the
#' q-values for each experiment.  Runs the sign test on the appropriate SNP for each experiment.
#'
#' @param experimentNames The names of the experiment to appear in the columns of the outputs.
#' @param permutedQTLFiles A vector of eQTL result file that has gone through permutation (or other gene-level correction) and has a q-value.
#' @param unfilteredQTLFiles A vector of eQTL result file run on all SNPs on the data set.  This data set should contain all SNP/Gene pairs in the permutedQTLFile,
#' even if they wouldn't normally be included in the data set.
#' @param qvalueThreshold A minimum q-value threshold for a gene from the permutedQTLFile to be included in the result set.
#' @param geneDescriptionsFile A file containing gene names in column 2 and gene descriptions in column 3.
#' @param outGeneLevelSummaryFile Output file containing the gene-level results
#' @param outExperimentLevelSummaryFile Output file containing the experiment level sign test results.
#' @param filterChromosomeList remove chromosomes in this vector.  Useful when you want to filter out the X chromosome.
#' @export

summarizeSignTestManyExperiments<-function (experimentNames, permutedQTLFiles, unfilteredQTLFiles, geneDescriptionsFile=NULL, qvalueThreshold=0.05, outGeneLevelSummaryFile=NULL, outExperimentLevelSummaryFile=NULL, filterChromosomeList=NULL) {
    #the union of genes and the q-values per gene/experiment
    addEffectSizes=F
    r=getQValuesAcrossExperiments(experimentNames, permutedQTLFiles, addEffectSizes=addEffectSizes, qvalueThreshold, geneDescriptionsFile, filterChromosomeList)

    #the sign tests.  Experiment level and gene level.
    d=lapply(permutedQTLFiles, testEffectsSignsMatchMany, unfilteredQTLFiles,qvalueThreshold,verbose=F, filterChromosomeList=filterChromosomeList)

    #set up the summary.
    summaryDF=do.call(rbind, lapply (d, function (x) x$summaryDF))
    rownames(summaryDF)=experimentNames
    colnames(summaryDF)=experimentNames
    if (!is.null(outExperimentLevelSummaryFile)) write.table(summaryDF, outExperimentLevelSummaryFile, row.names=T, col.names=T, quote=F, sep="\t")

    #set up the gene level
    dd=lapply(d, function (x) x$signTestDF)
    geneSignDF<-t(do.call(rbind, lapply(dd, function (x) x[match(r$gene, x$gene),]$signTestSummary)))
    colnames(geneSignDF)=paste(experimentNames, "sign test", sep=" ")

    #expression is the uniform, so get the global expression stats.
    exp=do.call(rbind, dd)
    exp=unique(exp[,c("gene", "expressionLevel")])
    idx=match(r$gene, exp$gene)
    exp=exp[idx,]

    #do all the signs match across experiments?
    dd=lapply(d, function (x) x$effectDirectionSame)
    edSame=t(do.call(rbind, lapply(dd, function (x) x[match(r$gene, x$gene),]$effectDirectionSame)))
    signsAgree=apply (edSame, 1, all, na.rm=T)

    #bind it all up.
    result=cbind (r, geneSignDF, medianExpression=exp$expressionLevel, "effect direction concordant"=signsAgree)

    #write it out.
    if (!is.null(outGeneLevelSummaryFile)) write.table(result, outGeneLevelSummaryFile, row.names=F, col.names=T, quote=F, sep="\t")
    return (list(experimentLevel=summaryDF, geneLevel=result))

}

#' Create the union of genes reported as significant in any eQTL experiment
#'
#' Report the union of genes that were significant (<= qvalueThreshold) in any experiment, and the q-value in each experiment.
#' @param addEffectSizes Emit the effect sizes of each significant eQTL.
#' @param filterChromosomeList remove chromosomes in this vector.  Useful when you want to filter out the X chromosome.
#' @inherit summarizeSignTestManyExperiments
#' @return A data frame with the gene, qvalues for experiments, and an optional gene description field.
getQValuesAcrossExperiments<-function (experimentNames, permutedQTLFiles, addEffectSizes=F, qvalueThreshold=0.05, geneDescriptionsFile=NULL, filterChromosomeList=NULL) {
    #gather up the union of genes and build the q-values part of the result.
    b=lapply (permutedQTLFiles, fread)

    if (!is.null(filterChromosomeList)) {
        b=lapply(b, function (x) xx=x[which(is.na(match(x$chr, filterChromosomeList))),])
    }

    allGenes=sort(unique(unlist(sapply(b, function (x) unique (x[x$qvalue<=qvalueThreshold,]$gene)))))
    qvalues<-t(do.call(rbind, lapply(b, function (x) x[match(allGenes, x$gene),]$qvalue)))
    colnames (qvalues)=paste(experimentNames, "q-value", sep=" ")

    #sort genes and qvalues by the min qvalue in the row.
    idx=order(apply (qvalues, 1, min, na.rm=T))
    qvalues=format(qvalues[idx,], digits=3, scientific=T)
    allGenes=allGenes[idx]

    #get the chromosome for each gene.
    z=do.call(rbind, b)
    z=z[match(allGenes, z$gene),]
    chromosomes=sapply (strsplit (z$SNP, ":"), function (x) x[1])

    if (addEffectSizes) {
        effectSizes<-t(do.call(rbind, lapply(b, function (x) x[match(allGenes, x$gene),]$beta)))
        effectSizes=round(effectSizes, digits=4)
        colnames (effectSizes)=paste(experimentNames, "effect size", sep=" ")
        #bind the effect sizes to the qvalues.
        qvalues=cbind (qvalues, effectSizes)
    }

    #get gene descriptions if appropriate.
    if (is.null(geneDescriptionsFile)) {
        df=data.frame(gene=allGenes, qvalues, check.names = F, stringsAsFactors = F)
        return (df)
    }

    z=parseGeneDescriptionFile(geneDescriptionsFile)
    idx=match(allGenes, z[,3])
    df=data.frame(gene=allGenes, chr=chromosomes, description=z[idx,2], qvalues, check.names = F, stringsAsFactors = F)
    idx=which(is.na(df$description))
    if (length(idx)>0) {
        df[idx,]$description=""
    }


    return (df)
}


parseGeneDescriptionFile<-function (geneDescriptionsFile) {
    z=fread(geneDescriptionsFile, data.table=F)
    q=strsplit (z[,2], split=" [", fixed=T)
    z[,2]=unlist(lapply(q, function (x) x[1]))
    return (z)
}


#unfilteredQTLFile="/downloads/emi/Adipose_Subcutaneous_Analysis_cis-eQTLs.txt.gz"
#unfilteredQTLFile="/downloads/emi/test.txt"

#permutedQTLFile="/downloads/emi/hg19_no_C4_Cumulative_NucleiDropulation_BA46_2019-12-12_r2_0.5_top50_glutamatergic_neurons.maf_0.05_cisDist_10kb.eQTL_eigenMT_results.txt"
#outMetricsFile="/downloads/emi/out.metrics.txt"
#qvalueThreshold=0.05; verbose=T; flipAlleles=T; stripENSGSuffix=T; chunk_size=1000;filterChromosomeList=NULL

#' Performs a sign test on two data sets to see if SNP/Gene combinations show the same direction of effect.
#'
#' @param permutedQTLFile An eQTL result file that has gone through permutation (or other gene-level correction) and has a q-value.
#' @param unfilteredQTLFile An eQTL result file run on all SNPs on the data set.  This data set should contain all SNP/Gene pairs in the permutedQTLFile,
#' even if they wouldn't normally be included in the data set.
#' @param qvalueThreshold A minimum q-value threshold for a gene from the permutedQTLFile to be included in the result set.
#' @param filterChromosomeList remove chromosomes in this vector.  Useful when you want to filter out the X chromosome.
#' @param verbose Should "helpful" information be printed as function is run?
#' @param chunk_size how much of the unfiltered eQTL file should be loaded at a time.  These files can be quite large.
#' @param outMetricsFile If not null, writes out a file with statistics on the run (number/fraction eGenes passing sign test)
#' @param flipAlleles If true, looks for instances where the reference and alternate alleles in the two data sets are reversed.  Sets
#' the reference and alternate allele in the unfilteredQTLFile to be the same as the permutedQTLFile, and where the alleles are flipped
#' the beta is multiplied by -1.
#' @param stripENSGSuffix Should the suffix be removed from ENSG IDs if they are used as the gene names?  For example
#' an ENSG can have a version number: ENSG00000000460.12 where 12 is the version.  If this is set to be true the ".12" is removed.
#' @return A list containing 2 elements: 1) A data frame containing a list of genes/SNPs that passed the FDR threshold, their q-value in permutedQTLFile,
#' the direction of the effect in the permutedQTLFile, and the direction in the unfilteredQTLFile.
#' 2) A calculation of the fraction of snp/gene pairs where the sign of the two results agree.
#' @import data.table readr
#' @export
testEffectsSignsMatch<-function (unfilteredQTLFile, permutedQTLFile, qvalueThreshold=0.05, verbose=T, filterChromosomeList=NULL,
                                 chunk_size=1000000, outMetricsFile=NULL, flipAlleles=T, stripENSGSuffix=T) {
    #unfilteredQTLFile=unfilteredQTLFiles[2];permutedQTLFile=permutedQTLFiles[1]
    if (!file.exists(unfilteredQTLFile)) stop(paste("Can't read file", unfilteredQTLFile))
    if (!file.exists(permutedQTLFile)) stop(paste("Can't read file", permutedQTLFile))

    a=suppressWarnings(fread(permutedQTLFile))
    #b=suppressWarnings(fread(unfilteredQTLFile))

    if (!is.null(filterChromosomeList)) a=a[which(is.na(match(a$chr, filterChromosomeList))),]

    a=a[a$qvalue<=qvalueThreshold,]
    numEQTLsOriginal=dim (a)[1]
    if (verbose) cat ("Number of eQTL SNPs", dim (a)[1], "\n")

    #readr filter generation and filter.
    if (length(grep(".gz", unfilteredQTLFile))==1) {
        cat (paste("Temp dir in use: ", tempdir(), "\n"))
    }

    addSNPLocationAndAllelesPrivate<-function (z) {
        s=strsplit (z$SNP, ":")
        #if there are encodings that miss the last allele, then let's add that allele back as a "-".
        #it'll probably be ignored, but at least we don't crash.
        idxMissingAllele=which(sapply(s, length)<4)
        for (i in idxMissingAllele) {
            #expected 4 places of info.
            pad=rep("-", 4-(length(s[[i]])))
            s[[i]]=c(s[[i]], pad)
        }
        d=do.call(rbind, s)
        colnames (d)=c("chr", "start", "ref", "alt")
        cbind (z, d)
    }


    filterChunkFullID<-function (x, pos, snpList) x[!is.na(match(x$SNP, snpList)),]
    filterChunkPosition<-function (x, pos, snpList) {
        x=addSNPLocationAndAllelesPrivate(x)
        snpIDs=paste(x$chr, x$start, sep=":")
        z=x[!is.na(match(snpIDs, snpList)),]
        return (z)
    }

    filterChunkFuncGen<-function (snpList) {
        function(x, pos) filterChunkPosition(x, pos, snpList)
    }


    if (!flipAlleles) {
        #Case: don't flip alleles.
        b=readr::read_delim_chunked(unfilteredQTLFile, readr::DataFrameCallback$new(purrr::partial(filterChunkFullID, snpList=a$SNP)), delim="\t", chunk_size =chunk_size, progress = T)
    } else {
        #flip alleles where neccesary
        a=addSNPLocationAndAllelesPrivate(a)
        snpList=paste(a$chr, a$start, sep=":")
        b=readr::read_delim_chunked(unfilteredQTLFile, readr::DataFrameCallback$new(purrr::partial(filterChunkPosition, snpList=snpList)), delim="\t", chunk_size =chunk_size, progress = T)
        b=flipRefAltAlleles(df=b,sites=a)
    }

    if (stripENSGSuffix) {
        a$gene=stripENSGSuffix(a$gene)
        b$gene=stripENSGSuffix(b$gene)
    }

    #figure out what's missing.
    missingGenes=setdiff(a$gene, b$gene)
    cat (paste("Untested Genes [", length(missingGenes), "]\n"))

    #key on the gene+snp.
    aKey=paste(a$gene, a$SNP, sep=":")
    bKey=paste(b$gene, b$SNP, sep=":")
    both=intersect(aKey, bKey)
    numComparedEqtls=length(both)

    #which SNPs are absent from the non-permuted data set?
    missingEqtls=as.vector(setdiff(aKey, both))
    names(missingEqtls)=NULL
    cat (paste("Untested eQTLs [", length(missingEqtls), "]\n"))

    #test if any of the SNP/Genes from permuted results aren't in the all-inclusive results.
    #aNotB=setdiff (aKey, bKey)
    #if (length(aNotB)>0)
    #    warning(paste("Not all of the permuted SNP/Genes appear in the unfiltered set: Num disagree [", length(aNotB),
    #    "] \nexamples [", paste(head (aNotB, 5), collapse=" "), "]"))

    aa=a[match(both, aKey),]
    bb=b[match(both, bKey),]

    #encode +/- direction of effect.
    getSignEncoding<-Vectorize(function (x) {
        if (is.na(x)) return (NA)
        if (x>0) return ("+")
        return ("-")
    })


    df=data.frame(gene=aa$gene, SNP=aa$SNP, qvalue=aa$qvalue, maf=aa$MAF,
                  betaPermuted=aa$beta, betaUnfiltered=bb$beta,
                  effectSizePermuted=aa$effect_size,
                  effectDirectionPermuted=getSignEncoding(aa$beta),
                  effectDirectionUnfiltered=getSignEncoding(bb$beta),
                  pvaluePermuted=aa$`p-value`,pvalUnfiltered=bb$`p-value`,
                  expressionPermuted=aa$median_expression, expressionUnfiltered=bb$median_expression,
                  stringsAsFactors=F, check.names=F)

    df$sign_test_result=df$effectDirectionPermuted==df$effectDirectionUnfiltered

    numSignsMatch=length(which(df$sign_test_result))
    fracSignsMatch=round (numSignsMatch/length(which(!is.na(df$effectDirectionUnfiltered))),4)

    dff=df[which(df$expressionPermuted>0 & df$expressionUnfiltered>0),]
    numSignsMatch2=length(which(dff$effectDirectionPermuted==dff$effectDirectionUnfiltered))
    fracSignsMatch2=round (numSignsMatch2/length(which(!is.na(dff$effectDirectionUnfiltered))),3)

    #do any of the SNPs match?
    testMullipleSNPsPerGene<-function (gene_name, df) {
        x=df[df$gene==gene_name,]
        matchIdx=which(x$effectDirectionPermuted==x$effectDirectionUnfiltered)
        anyMatch=length(matchIdx)>0
        allMatch=length(matchIdx)==dim (x)[1]
        o=data.frame(gene=gene_name, anyMatch=anyMatch, allMatch=allMatch, stringsAsFactors=F)
        return (o)
    }

    z=do.call(rbind, lapply(unique (df$gene), testMullipleSNPsPerGene, df))
    eGenesAnyMatchFrac=length(which(z$anyMatch==T))/dim (z)[1]
    eGenesAllMatchFrac=length(which(z$allMatch==T))/dim (z)[1]

    result=list(results=df, fracSignsMatch=fracSignsMatch, eGenesAnyMatchFrac=eGenesAnyMatchFrac,
                eGenesAllMatchFrac=eGenesAllMatchFrac, numEQTLsOriginal=numEQTLsOriginal, numComparedEqtls=numComparedEqtls)

    if (!is.null(outMetricsFile)) {
        dfStats=data.frame(unfilteredQTLFile=unfilteredQTLFile, permutedQTLFile=permutedQTLFile, qvalueThreshold=qvalueThreshold,
                           fracSignsMatch=fracSignsMatch, eGenesAnyMatchFrac=eGenesAnyMatchFrac, eGenesAllMatchFrac=eGenesAllMatchFrac,
                           numEQTLsOriginal=numEQTLsOriginal, numComparedEqtls=numComparedEqtls, stringsAsFactors = F)
        write.table(dfStats, outMetricsFile, row.names=F, col.names=T, quote=F, sep="\t")
    }
    return (result)

}

#' Generate sign test results for an eQTL data set as compared to a reference eQTL data set, binned by index SNP MAF.
#'
#' Bin sign test results by minor allele frequency and generate useful plots.
#'
#' @param permutedQTLFile An eQTL result file that has gone through permutation (or other gene-level correction) and has a q-value.
#' @param unfilteredQTLFile An eQTL result file run on all SNPs on the data set.  This data set should contain all SNP/Gene pairs in the permutedQTLFile,
#' even if they wouldn't normally be included in the data set.
#' @param qvalueThreshold A minimum q-value threshold for a gene from the permutedQTLFile to be included in the result set.
#' @param binSize The size of each bin of minor allele frequency.  This builds a sequence of bins from 0 to 0.5 by this size.
#' The sign test result is estimated for index SNPs within each bin.
#' @return A data frame containing the MAF bins and sign test results at each tested bin.
#' @seealso DropSeq.eqtl::testEffectsSignsMatch
#' @export
binSignTestByMAF<-function (permutedQTLFile, unfilteredQTLFile, qvalueThreshold=0.05, binSize=0.025) {

    getSignTestBin<-function (index, bins, d) {
        s=bins[index]
        e=bins[index+1]
        x=d[d$maf>s & d$maf<=e,]
        num=dim(x)[1]
        fracSignTest=length(which(x$sign_test_result==T))/dim(x)[1]
        if (is.nan(fracSignTest)) fracSignTest=NA
        df=data.frame(mafLow=s, mafHigh=e, numEqtls=num, fracSignTest=fracSignTest, stringsAsFactors = F)
        return (df)
    }

    #debug (DropSeq.eqtl::testEffectsSignsMatch)
    r=DropSeq.eqtl::testEffectsSignsMatch(unfilteredQTLFile, permutedQTLFile, qvalueThreshold=qvalueThreshold)
    d=r$results

    #binSize=0.025
    bins=seq(0, 0.5, by=binSize)

    df=do.call(rbind, lapply(1:(length(bins)-1), getSignTestBin, bins, d))
    labels=paste(df$mafLow, df$mafHigh, sep="-")

    old.par <- par(no.readonly = TRUE)

    par(mar=c(6, 4, 2, 2))
    par(mfrow=c(2,1))

    plot (1:length(labels), df$fracSignTest, xlab="", ylab="sign test",
          las=2, axes=F, ylim=c(0.5, 1))

    axis(side=2)
    axis(side=1, at=1:length(labels), labels=labels, las=2)
    mtext("minor allele frequency bin", side=1, line=6)
    title ("Sign test results binned by SNP allele frequency")

    barplot (df$numEqtls, xlab="", ylab="number of eQTLs", las=2, names.arg = labels)

    mtext("minor allele frequency bin", side=1, line=6)
    par(old.par)

    return (df)

}


#' Generate sign test results for an eQTL data set as compared to a reference eQTL data set, binned by q-value of the results.
#'
#' Bin sign test results by q-value and generate useful plots.
#'
#' @param permutedQTLFile An eQTL result file that has gone through permutation (or other gene-level correction) and has a q-value.
#' @param unfilteredQTLFile An eQTL result file run on all SNPs on the data set.  This data set should contain all SNP/Gene pairs in the permutedQTLFile,
#' even if they wouldn't normally be included in the data set.
#' @param qvalueThreshold A minimum q-value threshold for a gene from the permutedQTLFile to be included in the result set.
#' @param numBins The number of quantiles to partition data into.  The sign test result is estimated for index SNPs within each bin.
#' @return A data frame containing the MAF bins and sign test results at each tested bin.
#' @export
#' @seealso DropSeq.eqtl::testEffectsSignsMatch
binSignTestByQvalueEffectSize<-function (permutedQTLFile, unfilteredQTLFile, qvalueThreshold=0.05, numBins=5) {
	#for R CMD CHECK
	qValueQuantile=NULL;effectSizeQuantile=NULL;

    getSignTestBinByBin<-function (x) {
        num=dim(x)[1]
        fracSignTest=length(which(x$sign_test_result==T))/dim(x)[1]
        if (is.nan(fracSignTest)) fracSignTest=NA
        df=data.frame(numEqtls=num, fracSignTest=fracSignTest, stringsAsFactors = F)
        return (df)
    }

    #debug (DropSeq.eqtl::testEffectsSignsMatch)
    r=DropSeq.eqtl::testEffectsSignsMatch(unfilteredQTLFile, permutedQTLFile, qvalueThreshold=qvalueThreshold)
    d=r$results
    data.table::setDT(d)

    #breaks +1.  Each bin is breaks[i] to breaks[i+1]
    breaks=quantile(-log10(d$qvalue), probs = seq(0,1,length.out = (numBins+1)))
    d[,qValueQuantile := as.numeric(cut(-log10(qvalue), breaks = breaks, include.lowest = T)),]

    #by q-value
    z=d[,getSignTestBinByBin(.SD), by=qValueQuantile]
    z=z[order(z$qValueQuantile),]
    labels=sapply(1:(length(breaks)-1), function (idx, breaks) paste(round (breaks[idx],2), round (breaks[idx+1],2), sep="-"), breaks)
    z$qValueRange_log10=labels
    z=z[,c("fracSignTest", "numEqtls", "qValueRange_log10")]
    plot (1:dim(z)[1], z$fracSignTest, axes=F, xlab="q-value bin [-log10]", ylab="sign test result in bin", ylim=c(0, 1), pch=16, type='b')
    axis(2)
    axis(1, at=1:dim(z)[1], labels=z$qValueRange_log10)
    title ("sign test results by q-value group")

    #by effect size
    breaks=quantile(abs(d$effectSizePermuted), probs = seq(0,1,length.out = (numBins+1)))
    d[,effectSizeQuantile := as.numeric(cut(abs(d$effectSizePermuted), breaks = breaks, include.lowest = T)),]


    z2=d[,getSignTestBinByBin(.SD), by=effectSizeQuantile]
    z2=z2[order(z2$effectSizeQuantile),]
    labels=sapply(1:(length(breaks)-1), function (idx, breaks) paste(round (breaks[idx],2), round (breaks[idx+1],2), sep="-"), breaks)
    z2$effect_size_range=labels

    plot (1:dim(z2)[1], z2$fracSignTest, axes=F, xlab="effect size bin", ylab="sign test result in bin", ylim=c(0, 1), pch=16, type='b')
    axis(2)
    axis(1, at=1:dim(z2)[1], labels=z2$effect_size_range)
    title ("sign test results by effect size group")

    return (list(qvalue=z, effect_size=z2))

}

#ensgIDs=a$gene
stripENSGSuffix<-function (ensgIDs) {
    sapply(strsplit (ensgIDs, split=".", fixed=T), function (x) x[1])
}


#' Align the reference and alternate alleles of two sets of eQTL results so they match up
#'
#' Normally one would use the permuted data for the sites, and the unpermuted data for the df.
#  This reorients the unpermuted site alleles to match the permuted results
#'
#' @param df A set of eQTL results where some alleles may need to be flipped.
#' @param sites A set of eQTL results that act as a reference for what allele should be considered
#' the reference and alternate at a particular SNP.
#' @param verbose If true print extra debug info on the number of alleles flipped.
#'
#' @return The df dataframe with the ref/alt alleles oriented to match the sites dataframe.
#' SNPs in the df that don't have the same alleles as the sites SNPs are filtered.  Flipped SNPs have their
#' beta multiplied by -1 as well to reverse the direction of effect of the eQTL.
flipRefAltAlleles<-function (df, sites, verbose=T) {
    #df=b;sites=a
    df$key=paste(df$chr, df$start, df$ref, df$alt, sep=":")
    df$key_r=paste(df$chr, df$start, df$alt, df$ref, sep=":")

    sites$key=paste(sites$chr, sites$start, sites$ref, sites$alt, sep=":")

    sameOrientation=intersect(df$key, sites$key)
    oppositeOrientation=intersect(df$key_r, sites$key)

    z=df[!is.na(match(df$key, sameOrientation)),]
    z2=df[!is.na(match(df$key_r, oppositeOrientation)),]
    #flip the beta and the keys!
    if (dim (z2)[1]>0) {
        z2$beta=z2$beta*-1
        ref=z2$ref
        z2$ref=z2$alt
        z2$alt=ref
    }

    numMatch=dim(z)[1]
    numFlip=dim (z2)[1]
    numUnmatched= dim (df)[1]-(dim (z)[1]+dim (z2)[1])

    if (verbose) cat(paste("Sites match alleles [", numMatch, " ->",  round (numMatch/dim(df)[1]*100, 2),"%] Sites flipped alleles [", numFlip, " ->",  round (numFlip/dim(df)[1]*100, 2), "%] Sites Other [", numUnmatched, " ->",  round (numUnmatched/dim(df)[1]*100, 2), "%]\n"))

    df2=rbind (z, z2)
    setDT(df2)
    df2[,c("key", "key_r"):=NULL]
    setDF(df2)
    return (df2)
}



#' Run the sign test on multiple data sets and aggregate the results.
#'
#' Given one set of genes with permuted p-values, scan other data sets to see if their sign agrees.  Aggregate results of multiple unfiltered data sets into a single result.
#' @param unfilteredQTLFiles A list of unfiltered eQTL files
#' @param filterChromosomeList remove chromosomes in this vector.  Useful when you want to filter out the X chromosome.
#' @return A data frame containing the gene and a column containing the sign test data across experiments.  This is
#' the snpName : sign of experiment.  For example SNP123:+/+/+/+ would indicate 4 experiments where the direction is the same for all 4.
#' Also returns a single row data frame with the fraction of the genes from the permutedQTLFile where the sign matches in the unfiltered experiment.
#' @export
#' @inherit testEffectsSignsMatch
testEffectsSignsMatchMany<-function (permutedQTLFile, unfilteredQTLFiles, qvalueThreshold=0.05, verbose=T, filterChromosomeList=NULL) {
    #permutedQTLFile=permutedQTLFiles[1]
    d=lapply (unfilteredQTLFiles, testEffectsSignsMatch, permutedQTLFile, qvalueThreshold, verbose=T, filterChromosomeList)
    #the fraction of sign tests that match on each experiment.
    summaryDF=rbind (lapply(d, function (x) x$fracSignsMatch))
    dd=lapply(d, function (x) x$results)

    #gets all the unfiltered data effect directions.
    effectDirectionUnfiltered=t(do.call(rbind,(lapply (dd, function (x) x$effectDirectionUnfiltered))))
    #summarizes those effects.
    geneSummary=paste(dd[[1]]$SNP, "[", apply(effectDirectionUnfiltered, 1, function (x) paste(x, collapse="/")), "]", sep="")

    #get the expression levels per experiment from the unfiltered data.
    expressionLevelUnfiltered=t(do.call(rbind,(lapply (dd, function (x) round (x$expressionUnfiltered,2)))))
    expressionLevelUnfiltered=apply(expressionLevelUnfiltered, 1, function (x) paste(x, collapse=" / "))

    #generate an index when the sign direction matches across all experiments.
    effectDirectionSame=apply(effectDirectionUnfiltered, 1, function (x) length(unique(x))==1)
    edsDF=data.frame(gene=dd[[1]]$gene, effectDirectionSame, stringsAsFactors=F)

    df=data.frame(gene=dd[[1]]$gene, signTestSummary=geneSummary, expressionLevel=expressionLevelUnfiltered, stringsAsFactors = F)
    result=list(summaryDF=summaryDF, signTestDF=df, effectDirectionSame=edsDF)
    return (result)
}



# plot (aBest$beta, aBest$betaNew, xlab=paste(expNames[idxFiles[1]], "beta"), ylab=paste(expNames[idxFiles[2]], "beta"),  main=paste("% effect signs match=", round(fracSignsMatch*100,1)))
# points(aBest[idxBad,]$beta, aBest[idxBad,]$betaNew, col="red")
# abline (v=0)
# abline (h=0)
#
# plot (aBest$beta, aBest$betaNew, xlab=paste(expNames[idxFiles[1]], "beta"), ylab=paste(expNames[idxFiles[2]], "beta"),  main=paste("% effect signs match=", round(fracSignsMatch*100,1)), xlim=c(-1, 1), ylim=c(-1, 1))
# points(aBest[idxBad,]$beta, aBest[idxBad,]$betaNew, col="red")
# abline (v=0)
# abline (h=0)
#
# plot (aBest$beta, aBest$betaNew, xlab=paste(expNames[idxFiles[1]], "beta"), ylab=paste(expNames[idxFiles[2]], "beta"),  main=paste("% effect signs match=", round(fracSignsMatch*100,1)), xlim=c(-0.3, 0.3), ylim=c(-0.3, 0.3))
# points(aBest[idxBad,]$beta, aBest[idxBad,]$betaNew, col="red")
# abline (v=0)
# abline (h=0)
