# library (data.table)
# library(qvalue)
# library (MatrixEQTL)
# source ("~/dropseqrna/transcriptome/R/packages/DropSeq.eqtl/R/Matrix_eQTL_CommonFunctions.R")
# source ("~/dropseqrna/transcriptome/R/packages/DropSeq.eqtl/R/Matrix_eQTL_PermutationTesting.R")

# SNP_file_name="/downloads/eQTL/SMN/SMN_snp_genotypes.txt"
# snps_location_file_name="/downloads/eQTL/SMN/SMN_snp_locations.txt"
# phenotypeFileName="/downloads/eQTL/SMN/SMN_phenotype.txt"
# CNV_file_name="/downloads/eQTL/SMN/SMN_CN.txt"
# CNV_location_file_name="/downloads/eQTL/SMN/SMN_locations.txt"
# output_file_name_cis="/downloads/eQTL/SMN/SMN_sort_GWAS.linear"

# SNP_file_name="/broad/mccarroll/dmeyer/gwas/zikv2/genotypes.txt.gz"
# snps_location_file_name="/broad/mccarroll/dmeyer/gwas/zikv2/locations.txt"
# phenotypeFileName="/broad/mccarroll/dmeyer/gwas/zikv2/zikv.ug.wide.txt"
# CNV_file_name=NULL
# CNV_location_file_name=NULL
# output_file_name_cis="/broad/mccarroll/dmeyer/gwas/zikv2/zikv.ug.GWAS.permuted.linear"
# numPermutations=100;minNumPassingPermuations=15

#testSinglePhenotype(SNP_file_name, snps_location_file_name, phenotypeFileName, CNV_file_name, CNV_location_file_name, output_file_name_cis)

#' Test a single phenotype against all genotypes to find associations
#'
#' Take in a matrix of genotypes, snp locations, and a single phenotype formatted like an expression matrix input.
#' Optionally, supply a matrix of CNV genotypes and CNV locations.
#' Runs all genotypes/CNVs against the single expression vector.
#'
#' @param SNP_file_name A matrix of genotypes, 1 row per SNP, 1 column per sample.  Encoded as 0,1,2 copies of the alternate allele.
#' @param snps_location_file_name The location of each SNP.  Same number of lines as the SNP_file_name, encodes the snp name, chromosome, position.
#' @param phenotypeFileName A matrix of phenotypes, 1 row per phenotype, 1 column per sample.
#' @param CNV_file_name A matrix of CNV genotypes, 1 row per SNP, 1 column per sample.  Encoded as whole or dosage copy number.  Optional.
#' @param CNV_location_file_name The location of each CNV  Same number of lines as the SNP_file_name, encodes the snp name, chromosome, start position, end position.
#' @param output_file_name_cis The eQTL result file.
#' @param numPermutations The maximum number of permutations run on the result set.  Permutations are adaptive - for SNPs
#' where the permuted results quickly exceed the empiric pvalues, a smaller number of permutations are run to conserve time.
#' @param minNumPassingPermuations How many times do we need to see random data exceed the empiric pvalue to be satisfied with the permuted p?
#' @param chromosome If provided limit analysis to a single chromosome.
#' @export
#' @import MatrixEQTL data.table utils qvalue
testSinglePhenotype<-function (SNP_file_name, snps_location_file_name, phenotypeFileName,
                               CNV_file_name=NULL, CNV_location_file_name=NULL, output_file_name_cis,
                               numPermutations=10000, minNumPassingPermuations=15, chromosome=NULL) {

    validateSingleFileExists(SNP_file_name)
    validateSingleFileExists(snps_location_file_name)
    validateSingleFileExists(phenotypeFileName)
    if (!is.null(CNV_file_name)) validateSingleFileExists(CNV_file_name)
    if (!is.null(CNV_location_file_name)) validateSingleFileExists(CNV_location_file_name)

    #load up the data matrixes.
    #SNPs
    snps=DropSeq.utilities::fastRead(SNP_file_name, check.names=F)
    snpspos = data.table::fread(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);

    #phenotypes
    p=read.table(phenotypeFileName, header=T, stringsAsFactors=F, sep="\t", check.names = F)
    p=validatePhenotypeGenotypeInputs(snps, p)
    #optional CNV data.
    if (!is.null(CNV_file_name) & !is.null(CNV_location_file_name)) {
        cnv=  fread  (CNV_file_name)
        cnvspos=fread(CNV_location_file_name)
        idx=match(colnames (cnv), colnames (snps))
        cnv=cnv[,idx, with=F]
        snps=rbind(snps, cnv)
        cnvspos2=data.frame(snp=cnvspos$snp, chr=cnvspos$chr, pos=(cnvspos$e+cnvspos$s)/2, stringsAsFactors = F)
        snpspos=rbind(snpspos, cnvspos2)
    }

    #optionally filter on chromosome.
    if (!is.null(chromosome)) {
        cat (paste("Chromosome", chromosome, "selected, limiting set of eQTLs to permute", "\n"))
        idx=which(snpspos$chr==chromosome)
        snps=snps[idx,]
        snpspos=snpspos[idx,]
    }

    #convert to matrix eQTL format.
    p=convertToSlicedData(p)
    snps=convertToSlicedData(snps)

    #TODO: If we want to capture the SNPs with an empiric pvalue > 0.1, then break this work
    #up by chromosome and knit the results back together.  They aren't significant, so don't
    #need to be permuted later, which will solve the too many SNPs problem.
    invisible(capture.output(z<-MatrixEQTL::Matrix_eQTL_engine(snps=snps, gene=p, output_file_name=NULL, pvOutputThreshold=0.1, verbose=F), type="message"))
    zz=data.frame(z$all$eqtls, stringsAsFactors = F)
    colnames(zz)[match(c("snps", "statistic", "pvalue"), colnames(zz))]=c("SNP", "t-stat","P")
    idx=match(zz$SNP, snpspos$snp)
    zz$CHR=snpspos[idx,]$chr
    zz$BP=snpspos[idx,]$pos
    zz$CHR=sub("chr", "", zz$CHR)

    #output must be sorted in genomic order to be read by IGV.
    idx=which(zz$CHR=="X")
    if (length(idx)>0) zz[idx,]$CHR<-"23"

    idx=order(as.numeric(zz$CHR), as.numeric(zz$BP))
    zz=zz[idx,]

    idx=which(zz$CHR==23)
    if (length(idx)>0) zz[idx,]$CHR<-"X"

    empiricResults=data.frame(CHR=zz$CHR, BP=floor(zz$BP), SNP=as.character(zz$SNP), P=as.numeric(zz$P), stringsAsFactors=F)
    #run adaptive permutations.
    system.time (empiricResults<-runSinglePhenotypePermutationsAdaptive (snps, snpspos, p, empiricResults, numPermutations,minNumPassingPermuations, output_file_name_cis=output_file_name_cis))
    write.table(empiricResults, output_file_name_cis, row.names=F, col.names=T, quote=F, sep=" ")

}



#' Run permutation testing on the empiric results.
#'
#' Use adaptive permutation testing - test each SNP 100 times, then eliminate all results where the permuted result was
#' significant at least minNumPassingPermuations times.
#' @param snps Matrix eQTL package SlicedData object containing the genotype matrix (see: SNP_file_name)
#' @param snpspos data table containingthe SNP locations (see: snps_location_file_name)
#' @param p Matrix eQTL package SlicedData object containing the phenotype to test
#' @param empiricResults Dataframe containing the matrix eQTL empiric results for all SNPs in the snps object.
#' @param minNumPassingPermuations How many times do we need to see random data exceed the empiric pvalue to be satisfied with the permuted p?
#' @param output_file_name_cis If not null, write intermediate outputs to this file.
#' @inheritParams runSinglePhenotypePermutationsSimple
runSinglePhenotypePermutationsAdaptive<-function (snps, snpspos, p, empiricResults, numPermutations=10000, minNumPassingPermuations=15, output_file_name_cis) {
    #short circuit if no permutations requested
    if (numPermutations==0) return (empiricResults)

    #how many order of magnitude do we execute on?
    #maxOrderOfMagnitude=ceiling(log10(numPermutations))
    #hard codee max order of magnitude, then go to snp/gene pair individual testing.
    maxOrderOfMagnitude=4
    startOOM=2

    #set up starting conditions for data.
    empiricResults$excceds_empiric_p=0L
    empiricResults$total_permutations=0L

    for (oom in startOOM:maxOrderOfMagnitude) {
        #how many permutations should be run?
        totalPermutationsRun=max(empiricResults$total_permutations)
        maxPermutations=(10^oom)

        snpsToTest=unique(empiricResults[empiricResults$excceds_empiric_p<=minNumPassingPermuations,]$SNP)
        cat ("Adaptive permutation testing for log10 [", oom, "] Testing SNPs [", length(snpsToTest), "]\n")

        if (length(snpsToTest)==0) {
            cat ("No SNPs left to test.")
            break # break out of loop if no data left to test.
        }

        #Clone original genotype data and filter to SNPs to test this cycle.
        idx=which(!is.na(match(snps$GetAllRowNames(), snpsToTest)))
        snpsThis=snps$Clone()
        snpsThis$RowReorder(idx)

        #work on the subset of data.
        empiricResultsTested=empiricResults[match(snpsToTest, empiricResults$SNP),]

        for (i in (totalPermutationsRun+1):maxPermutations) {
            cat ("Permutation [", i, "/", maxPermutations,"]\n", sep="" )
            #these seem to take equivilent amounts of time
            empiricResultsTested=singlePhenotypePermutation(snpsThis, p, empiricResultsTested, verbose = F)
            #empiricResultsTested=singlePhenotypePermutationFaster(snpsThis, p, empiricResultsTested, verbose = F)
        }

        #update the full data with the permutation results from the subset
        idx=match(empiricResultsTested$SNP, empiricResults$SNP)
        if (any(is.na(idx)))
            stop ("NA index to match up permuted results with original!")
        empiricResults[idx,c("excceds_empiric_p", "total_permutations")]=empiricResultsTested[,c("excceds_empiric_p", "total_permutations")]

        #All permutations finished this round, update permuted pvalue.
        empiricResults$snp_permuted_pvalue=(empiricResults$excceds_empiric_p+1)/(empiricResults$total_permutations)
        #write.table(empiricResults, output_file_name_cis, row.names=F, col.names=T, quote=F, sep=" ")

    }

    #All permutations finished, update permuted pvalue.
    empiricResults$snp_permuted_pvalue=(empiricResults$excceds_empiric_p+1)/(empiricResults$total_permutations)
    #which SNPs have less than minNumPassingPermuations hits?
    idxPermuteMore=which(empiricResults$excceds_empiric_p<minNumPassingPermuations)

    if (length(idxPermuteMore)>0) {
        snpNames=empiricResults[idxPermuteMore,]$SNP
        phenotype=as.numeric (p$getSlice(1))
        empiricResults<-geneSNPPairPermutations (snpNames, empiricResults, phenotype, snps, numPermutations, minNumPassingPermuations, batchSize=1e6)
    }
    #write.table(empiricResults, output_file_name_cis, row.names=F, col.names=T, quote=F, sep=" ")
    return (empiricResults)

}






#' Permute the genotype data, then retest the genotype-phenotype relationship
#'
#' For each SNP, increment the number of tests by 1.
#' For each SNP where the permuted test P < empiric test P, increment excceds_empiric_p by 1.
#' @inheritParams runSinglePhenotypePermutationsAdaptive
#' @param empiricResultsTested The set of eQTL results to permute.
singlePhenotypePermutation<-function (snps, p, empiricResultsTested, verbose=T) {

    #this might be more efficient if the phenotype was the permuted instead of the genotype matrix.
    #this would prevent the copy by value R semantics when passing in the ojbect and modifying it.
    originalColNames=snps$columnNames

    snps$ColumnSubsample(sample (1:snps$nCols()))

    if (verbose) {
        result<-MatrixEQTL::Matrix_eQTL_engine(snps=snps, gene=p, output_file_name=NULL,
                                               pvOutputThreshold=1, verbose=F)
    } else {
        invisible(capture.output(result<-MatrixEQTL::Matrix_eQTL_engine(snps=snps, gene=p, output_file_name=NULL,
                                                                        pvOutputThreshold=1, verbose=F), type="message"))
    }

    #put the SNP back into the original column ordering
    idx=match(originalColNames, snps$columnNames)
    snps$ColumnSubsample(idx)

    #order the new eQTLs in the same order as the known best results, then count times the permuted result exceeds the old
    newEQTLs=result$all$eqtls
    newEQTLs=newEQTLs[match(empiricResultsTested$SNP, newEQTLs$snps),]

    #increment the number of tests, and capture the permuted tests that exceed the empiric pvalue
    empiricResultsTested$total_permutations=empiricResultsTested$total_permutations+1
    idxExceeds=which(newEQTLs$pvalue<empiricResultsTested$P)
    if (length(idxExceeds)>0)
        empiricResultsTested[idxExceeds,]$excceds_empiric_p=empiricResultsTested[idxExceeds,]$excceds_empiric_p+1

    return (empiricResultsTested)
}

#' Permute the phenotype data, then retest the genotype-phenotype relationship
#'
#' It's possible the single permutation is being slowed down by the pass by value semantics
#' of R when the genotype matrix is sampled and thus changed.  Since there's only
#' a single phenotype, permute that instead.
#'
#' For each SNP, increment the number of tests by 1.
#' For each SNP where the permuted test P < empiric test P, increment excceds_empiric_p by 1.
#' @inheritParams testSinglePhenotype
singlePhenotypePermutationFaster<-function (snps, p, empiricResultsTested, verbose=T) {

    #this might be more efficient if the phenotype was the permuted instead of the genotype matrix.
    #this would prevent the copy by value R semantics when passing in the ojbect and modifying it.
    originalColNames=p$columnNames

    p$ColumnSubsample(sample (1:p$nCols()))

    if (verbose) {
        result<-MatrixEQTL::Matrix_eQTL_engine(snps=snps, gene=p, output_file_name=NULL,
                                               pvOutputThreshold=1, verbose=F)
    } else {
        invisible(capture.output(result<-MatrixEQTL::Matrix_eQTL_engine(snps=snps, gene=p, output_file_name=NULL,
                                                                        pvOutputThreshold=1, verbose=F), type="message"))
    }

    #put the SNP back into the original column ordering
    idx=match(originalColNames, p$columnNames)
    p$ColumnSubsample(idx)

    #order the new eQTLs in the same order as the known best results, then count times the permuted result exceeds the old
    newEQTLs=result$all$eqtls
    newEQTLs=newEQTLs[match(empiricResultsTested$SNP, newEQTLs$snps),]

    #increment the number of tests, and capture the permuted tests that exceed the empiric pvalue
    empiricResultsTested$total_permutations=empiricResultsTested$total_permutations+1
    idxExceeds=which(newEQTLs$pvalue<empiricResultsTested$P)
    if (length(idxExceeds)>0)
        empiricResultsTested[idxExceeds,]$excceds_empiric_p=empiricResultsTested[idxExceeds,]$excceds_empiric_p+1

    return (empiricResultsTested)
}


#' Checks to see that all genotypes have a phenotype
#'
#' Reorders phenotype data if necessary to be in the same order as genotype data.
#'
#' @param snps A matrix of SNP genotypes (or any other data that has donors in columns)
#' @param p A matrix of phenotypes
#'
#' @return the phenotype matrix reordered to be have the same ordering as the genotype matrix.
#' @export
validatePhenotypeGenotypeInputs<-function (snps, p) {

    missingDonors=setdiff(colnames(snps), colnames(p))

    if (length(missingDonors)>0) {
        stop (paste("Not all genotyped donors have a matching phenotype column:",
                    paste(missingDonors, collapse=" ")))
    }

    #force phenotype and snp data into the same order by reordering the phenotype data.
    idx=match(colnames(snps), colnames(p))
    #add 1 to the start to retain the ID column.
    p=p[,c(idx),drop=F]
    #final dumb assertion.
    if (any(colnames (p)!=colnames (snps)))
        stop ("Could not resort phenotypes to match up to genotypes.")
    return (p)
}

#' Run permutation testing on the empiric results for the subset of SNPs with marginally significant p-values
#'
#' This is a pretty brute force approach without adaptive permutation testing.
#' @inheritParams testSinglePhenotype
runSinglePhenotypePermutationsSimple<-function (snps, snpspos, p, empiricResults,
                                                permutationPvalueThreshold=1e-3, numPermutations=10000) {

    snpsToTest=empiricResults[empiricResults$P<=pvalueThreshold,]$SNP

    r1=filterSNPDataSlice(snps, snpspos, snpsToTest)
    snps=r1$snps

    #to simplify merging in peremutation results, slim down empiric results to the subset
    #that is tested.
    empiricResultsTested=empiricResults[match(snpsToTest, empiricResults$SNP),]

    empiricResultsTested$excceds_empiric_p=0L
    empiricResultsTested$total_permutations=0L

    for (i in 1:numPermutations) {
        cat ("Permutation [", i, "/", numPermutations,"]\n", sep="" )
        empiricResultsTested=singlePhenotypePermutation(snps, p, empiricResultsTested)
    }
    empiricResultsTested$snp_permuted_pvalue=(empiricResultsTested$excceds_empiric_p+1)/(empiricResultsTested$total_permutations)
    #update the data with the permutation results.
    empiricResults$permuted_pvalue=NA
    empiricResults$total_permutations=NA
    idx=match(empiricResultsTested$SNP, empiricResults$SNP)
    empiricResults[idx,c("total_permutations", "permuted_pvalue")]=empiricResultsTested[,c("total_permutations", "snp_permuted_pvalue")]
    return (empiricResults)
}


##################################################################
# Single phenotype / snp permutation
##################################################################

#' Run adaptive permutation teeting on a list of SNPs with a single phenotype
#'
#' This permutes a single genotype / phenotype pair.  The pair is permuted a maximum of
#' numPermutations times.  Permutation can end early if the number of permuted results that exceed
#' the empiric p-value is >= minNumPassingPermuations, at which point the permuted pvalue is
#' considered to be well estimated.  This can significantly short circuit permutation of a single SNP.
#'
#' @param snpNames The names of the SNPs to test.
#' @param empiricResults The partially permuted results from matrix eQTL
#' @param phenotype A numeric vector of phenotype or expression data to test against all SNP genotypes
#' @param snps The SNP genotype matrix.
#' @param numPermutations The maximum number of permutations that should be run
#' @param minNumPassingPermuations The number of permutation results > empiric result required to stop permutation before reaeching numPermutations
#' @param batchSize how many permutations should be run each iteration.  This defines the tradeoff of memory vs CPU time.
geneSNPPairPermutations<-function (snpNames, empiricResults, phenotype, snps, numPermutations, minNumPassingPermuations=15, batchSize=1e6) {
    #snpName=snpNames[1]
    getOne<-function (index, snpNames, empiricResults, snps, phenotype, startOOM) {
        snpName=snpNames[index]
        g=as.numeric (snps$FindRow(snpName)$row)
        e=phenotype
        cat (paste(format(Sys.time()), " Adaptive permutation testing for SNP [", snpName, "] [" , index, "] of [", length(snpNames),"]\n", sep=""))
        permutedResult=geneSnpPairAdaptivePermutation(g, e, numPermutations, minNumPassingPermuations, batchSize, startOOM)
        z=empiricResults[empiricResults$SNP==snpName,]
        z[,c("excceds_empiric_p", "total_permutations", "snp_permuted_pvalue")]<-permutedResult[c("numSucceess", "numTrials", "permuted_pvalue")]
        return (z)
    }
    #determine the number of permutations already run on the data
    #if the max is 10^4, then start with 10^5.
    startOOM=floor(log10(max(empiricResults$total_permutations)))+1
    #run the rest of the permutations.
    #snpNames=snpNames[1:10]
    result=do.call(rbind, lapply(1:length(snpNames), getOne, snpNames, empiricResults, snps, phenotype, startOOM))
    #update the permuteed SNPs results with the new permutations
    idx=match(result$SNP, empiricResults$SNP)
    empiricResults[idx,]=result
    return (empiricResults)
}




#' Run adaptive permutation on one SNP Gene (or phenotype) pair.
#'
#' @param g The genotypes to test as the independent variable
#' @param e The expression value (or other response like a phenotype) to test.
#' @param numPermutations The maximum number of permutations
#' @param minNumPassingPermuations If an order of magnitude of permutations exceeed this number, stop permutation
#' @param batchSize how many permutations to run in a block.  Larger numbers use more memory
#' @return data frame of results.
geneSnpPairAdaptivePermutation<-function (g, e, numPermutations=1e8, minNumPassingPermuations=15, batchSize=1e6, startOOM=5) {

    numSucceess=0
    numTrials=0
    maxOrderOfMagnitude=ceiling(log10(numPermutations))

    #Find the startOOM from the maximum peremuted results empiric.

    for (oom in startOOM:maxOrderOfMagnitude) {
        #TODO eventually work the SNP and gene name in here?
        cat (paste(format(Sys.time()), " Adaptive permutation testing for log10 [", oom, "]" , "\n", sep=""))
        permutationsThisCycle=(10^oom)-numTrials

        result=geneSNPPairPermutation (g, e, permutationsThisCycle, batchSize)
        numTrials=numTrials+result$numTrials
        numSucceess=numSucceess + result$numSucceess
        #stop permutations if you found enough successes to estimate the pvalue well.
        if (numSucceess>=minNumPassingPermuations) break
    }

    lmResult=summary (lm (e ~ g))
    empiricP=lmResult$coefficients[2,4]

    if (numSucceess==0) numSucceess=1
    permuted_pvalue=numSucceess/numTrials
    df=data.frame(default_lm_pvalue=lmResult$coefficients[2,4], permuted_pvalue=permuted_pvalue,
                  numSucceess=numSucceess, numTrials=numTrials, stringsAsFactors = F)
    return (df)

}

#' Run the requested number of permutations distributed across batches
#' @return  Report the number of trials and successees.
#' @import permuco
#'
geneSNPPairPermutation<-function (g, e, nPermutations, batchSize=1e6) {
    # if the number of permutations is smaller than a batch...
    if (nPermutations<batchSize) batchSize=nPermutations
    numSucceess=0
    numTrials=0
    #cat ("Adaptive permutation runnning additional [",nPermutations, "] tests\n")
    while (numTrials<nPermutations) {
        # each iteration needs +1 permutation since the first result of the batch is the empiric.
        r=permuco::lmperm(e ~ g, np=batchSize+1, method="manly")
        #get rid of the empiric result from the permutations.
        d=r$distribution[-1]
        empiricT=r$table[2,3]
        numSucceess=numSucceess+length(which(abs(d)>abs(empiricT)))
        numTrials=numTrials+batchSize
        cat (paste(format(Sys.time()), " Permutation [" ,numTrials, "\\", nPermutations, "]\n", sep=""))
    }
    df=data.frame(numSucceess=numSucceess, numTrials=numTrials, stringsAsFactors = F)
    return (df)
}
