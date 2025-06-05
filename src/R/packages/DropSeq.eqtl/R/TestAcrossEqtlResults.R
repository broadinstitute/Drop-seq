#' Run a function across eQTL results
#'
#' @param rootDirs character vector of paths to root directories containing eQTL results
#' @param queryFile A tab delimited file with 2 columns.  SNP: the rsID of the SNP of interest.
#' @param FUN function called across query pairs for each rootDirs[i] i=1:length(rootDirs) 
#'            with signature function(gene = character(), geneExpression=numeric(), snp=character(), rsid=character(), 
#'            genotype=integer(), files=fetchEqtlFiles(rootDirs[i]))
#' @param usePeer boolean flag to indicate whether or not to use peer-adjusted expression
#'
#' @importFrom data.table fread
#' @importFrom IRanges IRanges
#' @export
doAcrossResults<-function (rootDirs, queryFile, FUN, usePeer) {
    query=fread(queryFile)

    # then for each rootDir, grab data for each gene/snp pair in queryFile do FUN
    for (i in 1:length(rootDirs)) {
        files=fetchEqtlFiles(rootDirs[i])
        # index genotype matrix if needed, otherwise return indexed bgzipped files
        gtMatrix=indexGenotypeMatrix(files$genotypeMatrixFile, forceIndex=T)
        if (usePeer) {
            if (!file.exists(files$peerExpressionFile)) {
                stop("Trying to use PEER-adjusted expression but PEER expression file not found")
            }
            dge=fread(files$peerExpressionFile)
        } else {
            dge=fread(files$expressionFile)
        }
        varLoc=fread(files$variantLocationsFile)
        for (j in 1:length(query)) {
            # grab out genotype from genotype_matrix.txt.bgz
            gt.res=with(varLoc[match(query$SNP[j], varLoc$id),],
                        queryTabixFile(gtMatrix, chr, IRanges(pos, width=1)))
            g=as.integer(gt.res[1,-c(1:3)])
            names(g)=names(gt.res)[-(1:3)]

            # next grab out expression from dge
            e=unlist(dge[id == query$GENE[j],-1])
            names(e)=names(dge)[-1]
            e=e[names(g)]

            # then run function with filled signature
            FUN(query$GENE[j], e, gt.res[1,3], query$SNP[j], g, files)
        }
    }
    invisible()
}

#' Reformat, bgzip, and index genotype matrix
#'
#' Writes out header as a single-line, tab-separated comment using the comment character '#'
#'
#' @param genotypeMatrixFile file path for genotype matrix in txt.gz format
#' @param forceBgzip forces block-gzip compression of genotype matrix if true, otherwise uses existing file
#' @param forceIndex forces re-indexing of block-gzipped genotype matrix if true, otherwise uses existing file
#' @return path to indexed bgzipped matrix
#'
#' @importFrom stringr str_replace
#' @importFrom Rsamtools bgzip indexTabix
#' @importFrom data.table fread
#' @importFrom readr write_tsv
#' @export
indexGenotypeMatrix<-function(genotypeMatrixFile, forceBgzip=F, forceIndex=F) {
    bgzMatrixFile=str_replace(genotypeMatrixFile, ".gz$", ".bgz")
    if (forceBgzip || !file.exists(bgzMatrixFile)) {
        tmp=tempfile()
        gtm=fread(genotypeMatrixFile)
        splitSnp=strsplit(gtm$id, ":")
        gtm=cbind(CHR=sapply(splitSnp, function(x) {x[1]}),
                  POS=sapply(splitSnp, function(x) {x[2]}),
                  #REF=splitSnp[3,], ALT=splitSnp[4,], 
                  SNP=gtm$id, gtm[,-1])

        # make header into a comment
        names(gtm)[1]<-paste0("#",names(gtm)[1])
        write_tsv(gtm, tmp, col_names=T, quote_escape='none')

        rm(gtm)
        bgzip(tmp, bgzMatrixFile, overwrite=TRUE)
        file.remove(tmp)
    }
    if (forceIndex || !file.exists(paste0(bgzMatrixFile, '.tbi')))
        indexTabix(file=bgzMatrixFile,seq=1L,start=2L,end=2L,skip=1L)
    return (bgzMatrixFile)
}

#' Fetch  a named list of eQTL files
#'
#' @param rootDir root directory path to search for eQTL files
#' @return named list of eQTL files with following names:
#'         eqtlResultsFile, indexSnpsFile, eigenMT.indexSnpsFile,
#'         eigenMT.resultsFile, expressionFile, peerExpressionFile
#'         genotypeMatrixFile, genotypeMatrixFileBgz, variantLocationsFile
#'
#' @importFrom stringr str_remove
#' @export
fetchEqtlFiles<-function(rootDir) {
    res=list()
    res$eqtlResultsFile=list.files(rootDir, ".eQTL_results.txt.gz", recursive=T, full.names=T)
    if (length(res$eqtlResultsFile) == 0) {
        stop("Could not find eQTL results for provided directory")
    } else if (length(res$eqtlResultsFile) > 1) {
        stop("More than 1 eqtl results files found for given directory")
    }

    prefix=str_remove(res$eqtlResultsFile, ".eQTL_results.txt.gz$")
    res$indexSnpsFile=paste0(prefix,".index_snps.txt")
    res$eigenMT.indexSnpsFile=paste0(prefix,".eigenMT.index_snps.txt")
    res$eigenMT.resultsFile  =paste0(prefix,".eQTL_eigenMT_results.txt")
    res$peerExpressionFile   =paste0(prefix,".gene_expression_peer.txt")
    res$expressionFile       =paste0(prefix,".gene_expression.txt")
    res$genotypeMatrixFile   =paste0(prefix,".genotype_matrix.txt.gz")
    res$genotypeMatrixFileBgz=paste0(prefix,".genotype_matrix.txt.bgz")
    res$variantLocationsFile =paste0(prefix,".variant_locations.txt")
    res$propertiesFile       =paste0(prefix,".properties")
    return(res)
}

#' Quickly query a tabix-indexed bgzipped tsv file
#'
#' @param tabixFile path to tabix-indexed bgzipped file
#' @param chr chromosome to query
#' @param irange object of class IRanges from Rsamtools defining positions to query
#' @param delim field delimiter of tabixFile defaults to tab
#' @param comment comment character used in tabixFile
#' @return data frame with results of tabix query
#'
#' @importFrom stringr str_remove str_split
#' @importFrom Rsamtools headerTabix scanTabix
#' @importFrom GenomicRanges GRanges
#' @export
queryTabixFile<-function(tabixFile, chr, irange, delim="\t",comment="#") {
    header=headerTabix(tabixFile)$header # pull out header
    header=str_remove(header,paste0("^",comment)) # remove comment character
    header=unlist(str_split(header, delim, simplify=T)) # convert to a vector
    res=scanTabix(tabixFile, param=GRanges(chr, irange))
    res=sapply(unlist(res), function(x) { unlist(str_split(x,delim,simplify=T)) })
    res=t(res)
    res=as.data.frame(res, stringsAsFactors=F)
    colnames(res)=header
    return(res)
}


#' Plot gene expression against a particular SNP using ggplot and ggbeeswarm
#'
#' @param gene name of gene being tested
#' @param snp string of the form {chrom}:{pos}:{ref}:{alt} specifying a SNP
#' @param rsid string giving RSID of SNP being tested (for showing in title/file name outputs)
#' @param genotype named integer vector representing the number of alt alleles at the given SNP for each donor
#' @param files named list of eQTl files returned from the function `fetchEqtlFiles`
#' @param outDir directory to write out any output files generated using naming scheme {paste0(outDir,'/',paste(gene,rsid,experimentName, sep="_")}
#' @param as.obj boolean flag that sets function to return ggplot object being plotted if set to TRUE. Will still write out if outDir or fullPath specified
#' @param fullPath path to output file used to overwrite default location relative to outDir (can leave outDir arg missing if this is specified)
#' @param experimentName if not specified, will use EXPERIMENT field of properties file
#' @return only if as.obj == TRUE, ggplot object
#'
#' @import ggplot2
#' @import ggbeeswarm
#' @importFrom stringr str_split str_glue
#' @export
ggplotGenotypesVsExpression<-function(gene, geneExpression, snp, rsid, genotype, files, outDir, as.obj=F, fullPath=NULL, experimentName=NULL) {
    data=data.frame(GT=as.integer(genotype), Donor=names(genotype), Expression=geneExpression, stringsAsFactors = F)
    splitSnp=as.character(unlist(str_split(snp, ":", simplify=T)))
    ref=splitSnp[3]
    alt=splitSnp[4]
    alleles=c(paste0(ref,'/',ref),paste0(ref,'/',alt),paste0(alt,'/',alt))
    data$alleles=factor(alleles[data$GT+1], alleles)
    fit=lm(Expression ~ GT, data)
    fit.slope=coef(fit)[[2]]
    fit.intercept=coef(fit)[[1]]
    fit.effect_verb=ifelse(abs(fit.slope) >= 0, "increase", "decrease")
    fit.effect_pct = fit.slope / median(data$Expression[data$GT == 1]) * 100
    pval=summary(fit)$coefficients[2,4]

    props=readProperties(propFile=files$propertiesFile)

    if (is.null(experimentName))
        experimentName=props$value[props$key=="EXPERIMENT"]

    plottitle=paste("{experimentName}",
                    "SNP={rsid} Gene={gene}",
                    "R2 = {round(summary(fit)$r.squared,3)}",
                    "p = {format.pval(pval,3)}",sep="\n")
    if (pval < 0.05) {
        plottitle=paste(
                        plottitle,
                        "beta = {round(fit.slope,2)}",
                        "{round(fit.effect_pct, 0)}% {fit.effect_verb} per allele",
                        sep="\n")
    }

    plottitle=str_glue(plottitle)
    ylims=c(floor(min(data$Expression)),
            ceiling(max(data$Expression)))

    p=ggplot(data, aes(x=GT, y=Expression))+
        geom_quasirandom(shape=1,size=3)+
        theme_classic()+
        labs(x='',title = plottitle)+
        scale_x_continuous(breaks = c(0,1,2), labels=alleles)+
        lims(y=ylims)+
        theme(title = element_text(size = 14),
              axis.text = element_text(size = 16),
              axis.title = element_text(size = 16))
    if (pval < 0.05) {
        p=p+geom_abline(slope = fit.slope, intercept=fit.intercept,lty='dashed', col='red')
    }

    if (!missing(outDir) && is.null(fullPath)) {
        ggsave(paste0(outDir,'/',paste(gene,rsid,experimentName, sep="_"), '.png'), p, height = 6, width = 6, units = "in")
    } else if (!is.null(fullPath)) {
        ggsave(fullPath, p, height = 6, width = 6, units = "in")
    }
    if (as.obj) {
        return(p)
    } else {
        invisible()
    }
}

#' Given a vector of root directories and a queryFile, runs set of
#' SNP/Gene queries across eQTL results contained within each root directory
#'
#' @param rootDirs character vector of directories containing one set of eQTL results
#' @param queryFile A tab delimited file with 2 columns.  SNP: the rsID of the SNP of interest.
#' @param outDir output directory to write out files into. Files are written out
#' @param usePeer specifies whether or not to use peer-adjusted expression data for testing queries
#'
#' @importFrom purrr partial
#' @export
plotGenotypesVsExpressionAcrossResults<-function(rootDirs, queryFile, outDir, usePeer=T) {
    if (!dir.exists(outDir)) {
        if (dir.exists(dirname(outDir))) dir.create(outDir)
        else stop("Output directory does not exist and cannot be created")
    }
    rootDirs.exists=sapply(rootDirs, dir.exists)
    if (!all(rootDirs.exists)) {
        stop(paste0("The following rootDirs do not exist: [", paste(rootDirs[!rootDirs.exists], collapse=", "),"]"))
    }
    if (!file.exists(queryFile)) stop("queryFile provided does not exist")
    doAcrossResults(rootDirs, queryFile, partial(ggplotGenotypesVsExpression, outDir=outDir), usePeer=usePeer)
}

#library(DropSeq.eqtl, lib.loc="/broad/mccarroll/software/dropseq/dm_branch/R_LIBS")
#rootDirs=c("/broad/mccarroll/dropulation/analysis/eQTL/CHDI/anteriorCaudate/CHDI.anteriorCaudate.raw.dge.SPN.peer10",
#           "/broad/mccarroll/dropulation/analysis/eQTL/CHDI/anteriorCaudate/CHDI.anteriorCaudate.raw.dge.astrocyte.peer10",
#           "/broad/mccarroll/dropulation/analysis/eQTL/CHDI/anteriorCaudate/CHDI.anteriorCaudate.raw.dge.endothelia.peer10",
#           "/broad/mccarroll/dropulation/analysis/eQTL/CHDI/anteriorCaudate/CHDI.anteriorCaudate.raw.dge.interneuron.peer10",
#           "/broad/mccarroll/dropulation/analysis/eQTL/CHDI/anteriorCaudate/CHDI.anteriorCaudate.raw.dge.microglia.peer10",
#           "/broad/mccarroll/dropulation/analysis/eQTL/CHDI/anteriorCaudate/CHDI.anteriorCaudate.raw.dge.oligodendrocyte.peer10",
#           "/broad/mccarroll/dropulation/analysis/eQTL/CHDI/anteriorCaudate/CHDI.anteriorCaudate.raw.dge.polydendrocyte.peer10")
#rootDirs=paste0(rootDirs, "/maf_0.05_cisDist_10kb")
#queryFile="/broad/mccarroll/dmeyer/analysis/SteveOneoffs/snp_gene_pairs.txt"
#plotOutputDir="/broad/mccarroll/dmeyer/analysis/SteveOneoffs/queryPairPlots"
#for (i in 1:length(rootDirs))
#  plotGenotypesVsExpressionAcrossResults(rootDirs[i],queryFile, plotOutputDir)
