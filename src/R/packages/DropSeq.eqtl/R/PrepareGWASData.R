#normalize and merge GWAS data sources to a standard format.
#highly influenced by ldscore munge_sumstats, but I want more columns in my output,
#and want to merge outputs and have genomic sort on output.

# Target format should have columns (prefer CHR, POS, SNPID, A1, A2, P, Z) for final column names/order.
#     chr        snpid        a1      a2      pos     Z
#     1       rs4951859       C       G       729679  -1
#     1       rs142557973     T       C       731718  1
#     1       rs141242758     T       C       734349  1.5


# library (data.table)
#may need to lift data over after other format modifications
#snpLocs=liftOverCoordinates(snpLocs, inEqtlFileSequenceDictionaryFile, outputSequenceDictionaryFile, liftOverChainFile, picardPath)


# inFiles=c("/downloads/GWAS_correlation/daner_PGC_SCZ_1016.f.b.GCadjusted.c", "/downloads/GWAS_correlation/SCZ_Index_SNPs.txt")
# outFile="/downloads/GWAS_correlation/daner_PGC_SCZ_1016_plus_INDEX_SNPS.txt"

# inFiles=c("/downloads/GWAS_correlation/ckqny.scz2snpres.gz")
# outFile="/downloads/GWAS_correlation/ckqny.scz2snpres.prepped.txt"

# inFile=c("/downloads/GWAS_correlation/ckqny.scz2snpres.gz")
# outFile="/downloads/GWAS_correlation/ckqny.scz2snpres.hg19.imputation_0.9.txt"
# imputationThreshold=0.9
# inSequenceDictionaryFile=NULL;outputSequenceDictionaryFile=NULL;liftOverChainFile=NULL;picardPath=NULL

#note, because liftover expects hg19 contigs to have "chr" in them, the input seqeunce dictionary needs that too for hg19->hg38

# inFile=c("/downloads/PGC3_SCZ_wave3_public_without_frequencies.v1.tsv")
# outFile="/broad/mccarroll/dropulation_census_metadata/GWAS/PGC3_SCZ_wave3_public_without_frequencies.v1.hg38.txt"
# imputationThreshold=0.9
# inSequenceDictionaryFile="/downloads/Frazer/hg19_UCSC.dict"
# outputSequenceDictionaryFile="/downloads/Frazer/GRCh38_maskedAlt.dict"
# liftOverChainFile="/downloads/Frazer/hg19ToHg38.over.chain.gz"
# picardPath="/broad/mccarroll/software/dropseq/priv/3rdParty/picard/picard.jar"



#' Convert GWAS data to a standard format
#'
#' Convert GWAS data to a standard format that can be read in by LDScore, used in GWAS correlation analysis,
#' and is suitable to be indexed by Tabix.
#'
#' The inSequenceDictionaryFile, outputSequenceDictionaryFile, liftOverChainFile, and picardPath can be left null
#' unless you desire to lift over your GWAS results to a different build.
#'
#' @param inFile The input GWAS file - this accepts a variety of standard formats.
#' @param outFile The GWAS data, reformatted, filtered, lifted over, and sorted for indexing.
#' @param imputationThreshold The minimum imputation threshold to retain a marker
#' @param inSequenceDictionaryFile The sequence dictionary file containing the contigs for input file.
#' @param outputSequenceDictionaryFile The sequence dictionary file containing the contigs for the output file
#' @param liftOverChainFile The chain file that maps position mapping from the old to new genome builds
#' @param picardPath Where the picard jar file is located.
#' @export
prepareGWASData<-function (inFile, outFile, imputationThreshold=0.9, inSequenceDictionaryFile=NULL,
                           outputSequenceDictionaryFile=NULL, liftOverChainFile=NULL, picardPath=NULL) {

    #d=lapply(inFile, prepareGWASDataSingle, imputationThreshold)
    #d=do.call(rbind, d)
    d=prepareGWASDataSingle(inFile, imputationThreshold)
    d=removeDuplicateSNPs(d)

    d=sortGWASData(d)
    #capitalize any columns that aren't already.
    colnames(d)=toupper(colnames(d))

    #lift over if all the proper files are in place.
    if (!is.null(inSequenceDictionaryFile) & !is.null(outputSequenceDictionaryFile) & !is.null(liftOverChainFile) & !is.null(picardPath)) {
        d=liftoverGWAS(d, inSequenceDictionaryFile, outputSequenceDictionaryFile, liftOverChainFile, picardPath)
        #need to resort if coordinates lifted over.
        d=sortGWASData(d)
    }

    write.table(d, outFile, row.names = F, col.names = T, quote=F, sep="\t")
}



liftoverGWAS<-function (d, inSequenceDictionaryFile, outputSequenceDictionaryFile, liftOverChainFile, picardPath) {
    #get the coordinates for liftover
    snpLocs=unique (data.frame(chr=d$CHR, start=d$POS, end=d$POS,marker_id=d$SNP, stringsAsFactors = F))
    snpLocs=liftOverCoordinates(snpLocs, inSequenceDictionaryFile, outputSequenceDictionaryFile, liftOverChainFile, picardPath)
    idx=match(d$SNP, snpLocs$marker_id)
    d$CHR=snpLocs[idx,]$new_chr
    d$POS=snpLocs[idx,]$new_start
    #purge any NA results where the SNP didn't map.
    d=d[!is.na(d$CHR),]
    return (d)
}


#rs2220276
removeDuplicateSNPs<-function (d) {
    #If merged studies have the same SNP with same allele and different p-values, need to select the minimal P.
    #check alleles in both orientations in case a merged study uses the other orientation.
    key=paste(d$CHR, d$POS, d$A1, d$A2, sep=":")
    key2=paste(d$CHR, d$POS, d$A2, d$A1, sep=":")
    keys=c(key,key2)
    dupes=keys[which(duplicated(keys))]
    idxDupes=which(!is.na(match(key, dupes)))
    if (length(idxDupes)==0) return (d)

    #partition out the non-dupe records.
    dOK=d[-idxDupes,]
    dupeRecs=d[idxDupes,]
    #only work off the chromosome and position, the allele flipping was checked above.
    dupeRecs$key=paste(dupeRecs$CHR, dupeRecs$POS, sep=":")
    setDT(dupeRecs)

    getBest<-function (x) {
        idx=which.min(x$P)
        x[idx,]
    }

    r=dupeRecs[,getBest(.SD),by=key]
    r[,key:=NULL]
    setcolorder(r, colnames(dOK))
    r=rbind (dOK, r)
    return (r)
}

sortGWASData<-function (d) {
    #this is annoying primarily because of chrX, chrY, chrMT
    chr=d$CHR
    chr=sub("chr", "", chr)
    chr[chr=="X" | chr=="chrX"]<-23
    chr[chr=="Y" | chr=="chrY"]<-24
    chr[chr=="MT" | chr=="chrMT" | chr=="M" | chr=="chrM"]<-25

    idx=order(chr, d$POS)
    return (d[idx,])
}

#inFile=inFiles[1]
prepareGWASDataSingle<-function (inFile, imputationThreshold=0.9) {
    a=read.table(inFile, header=T, stringsAsFactors = F, sep="\t")
    dataColNames=toupper(colnames(a))

    #pre-process alleles if both alleles are in one column.
    a=splitAlleles(a)

    #filter data by info field
    a=filterByInfo(a, dataColNames, imputationThreshold)

    chrIdx=getChromosomeIndex(dataColNames)
    posIdx=getPositionIndex(dataColNames)
    snpIdx=getSNPIdx(dataColNames)
    a1Idx=getA1Index(dataColNames)
    a2Idx=getA2Index(dataColNames)
    pvalIdx=getPValueIdx(dataColNames)

    signedStat=getSignedStatIndexAndName (dataColNames)
    zScore=getZscore(a, signedStat, pvalIdx)

    df=data.frame(a[,c(chrIdx, posIdx, snpIdx, a1Idx, a2Idx, pvalIdx, signedStat$idx)])
    #need to fix colnames.
    colnames(df)[1:6]=c("CHR", "POS", "SNP", "A1", "A2", "P")
    df=cbind(df, Z=zScore)
    return (df)
}

getZscore<-function (a, signedStat, pvalIdx) {
    p=a[,pvalIdx]
    ss=a[,signedStat$idx]
    getThresholdForSignedStat<-function (signedStatName) {
        if (signedStatName=="OR") return (1)
        if (signedStatName=="LOG_ODDS" | signedStatName=="BETA" | signedStatName=="Z") return (0)
        else stop ("No recognized signed stat found!")
    }
    ssThreshold=getThresholdForSignedStat(signedStat$name)
    zScores=round (sapply(p, getZ),6)
    idxNegative=which(ss<ssThreshold)
    if (length(idxNegative)>0) zScores[idxNegative]=zScores[idxNegative]*-1
    return (zScores)
}

#get Z-score.
getZ<-function (pval) -1* qnorm(pval/2,F)


##############################
# PARSE OUT THE REQUIRED COLUMNS.
###############################

getGenericIndex<-function (dataColNames, possibleNames, colName, required=F) {

    idx=which(!is.na(match(dataColNames, possibleNames)))
    if (length(idx)>1)
        warning (paste ("Multiple ", colName, "possible columns found", paste(dataColNames[idx], collapse=", ")))
    if (length(idx)==0) {
        strMsg=paste("Searching for", colName, "could not find in [", paste(dataColNames, collapse = " "),"]\n")
        if (required)
            stop(strMsg)
        else
            cat (strMsg)
    }
    else {
        strMsg=paste("Searching for", colName, "found as [", dataColNames[idx], "]\n")
        cat (strMsg)
    }

    return (idx)
}

getChromosomeIndex=function (dataColNames) {
    chrNames=c("CHR", "HG19CHRC")
    getGenericIndex(dataColNames, possibleNames = chrNames, colName="Chromosome", required=T)
}

#what column has the position information?
getPositionIndex<-function (dataColNames) {
    posNames=c("BP", "POS")
    getGenericIndex(dataColNames, possibleNames = posNames, colName="Position", required=T)
}

getA1Index<-function (dataColNames) {
    a1Names=c("A1", "ALLELE1", "ALLELE_1", "EFFECT_ALLELE", "REFERENCE_ALLELE", "INC_ALLELE", "EA")
    getGenericIndex(dataColNames, possibleNames = a1Names, colName="Allele 1", required=T)
}

getA2Index<-function (dataColNames) {
    a1Names=c("A2", "ALLELE2", "ALLELE_2", "OTHER_ALLELE", "NON_EFFECT_ALLELE", "DEC_ALLELE", "NEA")
    getGenericIndex(dataColNames, possibleNames = a1Names, colName="Allele 2", required=T)
}

getSNPIdx<-function (dataColNames) {
    snpNames=c("SNP", "MARKERNAME", "SNPID", "RS", "RSID", "RS_NUMBER", "RS_NUMBERS")
    getGenericIndex(dataColNames, possibleNames = snpNames, colName="SNP ID", required=T)
}

getPValueIdx<-function (dataColNames) {
    pvalNames=c("P", "PVALUE", "P_VALUE", "PVAL", "P_VAL", "GC_PVALUE")
    getGenericIndex(dataColNames, possibleNames = pvalNames, colName="SNP PVALUE", required=T)
}

getSignedStatIndexAndName<-function (dataColNames) {
    zScore=c("ZSCORE", "Z-SCORE", "GC_ZSCORE", "Z")
    or=c("OR")
    beta=c("B", "BETA", "EFFECTS", "EFFECT")
    logOdds=c("LOG_ODDS")

    zScoreIdx=getGenericIndex(dataColNames, possibleNames = zScore, colName="Z score", required=F)
    orIdx=getGenericIndex(dataColNames, possibleNames = or, colName="Odds Ratio", required=F)
    betaIdx=getGenericIndex(dataColNames, possibleNames = beta, colName="Beta", required=F)
    logOddsIdx=getGenericIndex(dataColNames, possibleNames = logOdds, colName="Log Odds", required=F)

    #check for multiple signed stats.
    allSignedStatsIdx=c(zScoreIdx, orIdx, betaIdx, logOddsIdx)
    if (length(allSignedStatsIdx)>1)
        warning (paste("Multiple signed stats found", paste(dataColNames[allSignedStatsIdx], collapse=", "), collapse=""))

    if (length(allSignedStatsIdx)==0)
        stop ("No signed stats found!")

    if (length(zScoreIdx)==1) return (list(name="Z", idx=zScoreIdx))
    if (length(orIdx)==1) return (list(name="OR", idx=orIdx))
    if (length(betaIdx)==1) return (list(name="BETA", idx=betaIdx))
    if (length(logOddsIdx)==1) return (list(name="LOG_ODDS", idx=logOddsIdx))

}


######################
# SPECIAL PURPOSE FILTER
############################
#need check for a column named "A1A2", which has both alleles, need to split them up into two columns.
splitAlleles<-function (a) {
    idx=getGenericIndex(colnames(a), possibleNames = c("A1A2"), colName="A1A2", required=F)
    if (length(idx)==0) {
        cat ("A1A2 column not found, skipping split\n")
        return (a) #no-op
    }
    z=strsplit (a[,idx], split="/", fixed=T)
    a1=sapply(z, function (x) x[1])
    a2=sapply(z, function (x) x[2])
    a$A1=a1
    a$A2=a2
    return (a)
}

filterByInfo<-function (a, dataColNames, infoThreshold=0.9) {
    infoNames=c("INFO")
    idxCol=getGenericIndex(dataColNames, possibleNames = infoNames, colName="INFO")
    idx=which(a[,idxCol]>=infoThreshold)
    a=a[idx,]
    return (a)
}

liftOverCoordinates<-function (snpLocs, inSequenceDictionaryFile, outputSequenceDictionaryFile, liftOverChainFile, picardPath) {
    #you need to append "chr" to the front of the contig names because someone hates me.
    #only if chr isn't already there.
    if (length(grep ("chr", snpLocs$CHR))==0) {
        snpLocs$chr=paste("chr", snpLocs$chr, sep="")
    }
    outIntervalFile=tempfile(pattern="output_", fileext=".interval")
    outLifted=tempfile(pattern="output_snps_", fileext=".lifted.txt")
    outRejected=tempfile(pattern="output_snps_", fileext=".rejected.txt")
    if (length(unique (snpLocs$marker_id))!=length(snpLocs$marker_id)) {
        stop("Marker IDs are not unique, can't lift over unambiguously")
    }

    writeIntervalFile(snpLocs, inSequenceDictionaryFile, outIntervalFile)

    #java -jar /broad/mccarroll/software/dropseq/priv/3rdParty/picard/picard.jar LiftOverIntervalList I=/var/folders/4s/c0wshbjd24s6771krpk5q2shl5h3l6/T//RtmpiVyjzg/hg19_7b267ed9adf3.interval CHAIN=/downloads/Frazer/hg19ToHg38.over.chain.gz SD=/downloads/Frazer/GRCh38_maskedAlt.dict O=lifted.txt REJECT=rejected.txt

    #O=lifted.txt REJECT=rejected.txt
    cmd=paste("java -jar ", picardPath, " LiftOverIntervalList I=", outIntervalFile,
              " CHAIN=", liftOverChainFile, " SD=", outputSequenceDictionaryFile,
              " OUTPUT=", outLifted, " REJECT=", outRejected, sep="")
    system(cmd, ignore.stderr = T, ignore.stdout = T)
    #system(cmd)
    lifted=read.table(outLifted, comment.char = "@", stringsAsFactors = F)
    #rejected=read.table(outRejected, comment="@")

    #the interval of the lifted file has the marker ID, just need to look that up.
    idx=match(snpLocs$marker_id, lifted$V5)

    snpLocs$new_chr=as.character(lifted[idx,]$V1)
    snpLocs$new_start=lifted[idx,]$V2
    snpLocs$new_end=lifted[idx,]$V3

    #get rid of NA entries
    idxNA=which(is.na(snpLocs$new_chr))
    if (length(idxNA)>0) {
        cat ("Dropped [", length(idxNA), "] SNP sites due to liftover failure.\n")
        snpLocs=snpLocs[-idxNA,]
    }
    file.remove(c(outIntervalFile, outLifted, outRejected))
    #remove the "chr" added to the contig names for hg19 data, which shouldn't have it.
    snpLocs$chr=sub("chr", "", snpLocs$chr)
    return (snpLocs)
}

writeIntervalFile<-function (dfLoc, dictFile, outIntervalFile, addChr=F) {
    a=read.table(dictFile, header=F, stringsAsFactors = F, sep="\t", fill=T)
    if (addChr) {
        idx=grep ("SN", a$V2)
        a[idx,]$V2=sub(":", ":chr", a[idx,]$V2)
    }

    write.table(a, outIntervalFile, row.names=F, col.names = F, quote=F, sep="\t")
    z=data.frame(chr=dfLoc$chr, start=as.integer(dfLoc$start), end=as.integer(dfLoc$end), strand="+", name=dfLoc$marker_id)
    write.table(z, outIntervalFile, row.names=F, col.names = F, quote=F, sep="\t", append=T)
}


#perhaps converting to a Z-score up front would be better, because I can't guarantee if I'll see OR, beta, etc.
# filter on info column, MAF.  Thresholds: INFO > 0.9, MAF > 0.01


#how to map from signed statistic to Z-score.
# null_values = {
#
#     'LOG_ODDS': 0,
#     'BETA': 0,
#     'OR': 1,
#     'Z': 0
# }


#Mapping from different column names to cannonical name.

# MAF
# 'EAF': 'FRQ',
# 'FRQ': 'FRQ',
# 'MAF': 'FRQ',
# 'FRQ_U': 'FRQ',
# 'F_U': 'FRQ',



