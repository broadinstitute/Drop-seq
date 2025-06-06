# library (LDlinkR)
#
# ldTargetPopulation="EUR"
#
# rootDir="/broad/mccarroll/dropulation/analysis/eQTL/CHDI/anteriorCaudate"
# queryFile="/broad/mccarroll/nemesh/eQTL_Fetch/query.txt"
# outFile="/broad/mccarroll/nemesh/eQTL_Fetch/CHDI_eQTL_results.txt"
# outSummaryFile="/broad/mccarroll/nemesh/eQTL_Fetch/CHDI_eQTL_results.summary.txt"

# rootDir="/broad/mccarroll/dropulation/analysis/eQTL/BA46"
# queryFile="/broad/mccarroll/nemesh/eQTL_Fetch/query.txt"
# #outFile="/broad/mccarroll/nemesh/eQTL_Fetch/BA46_eQTL_results.txt"
# outFile="/broad/mccarroll/nemesh/eQTL_Fetch/BA46_eQTL_results_cleaned.txt"
# outSummaryFile="/broad/mccarroll/nemesh/eQTL_Fetch/BA46_eQTL_results.summary.txt"

#z=read.table("/broad/mccarroll/nemesh/eQTL_Fetch/BA46_eQTL_results.txt", header=T, stringsAsFactors=F)
#zz=unique (z[!is.na(z$explicit_pair_pvalue),c("gene", "rsID")])
#options(warn=2)

#' Aggregate results across many eQTL data sets for a list SNPs/Genes pairs
#'
#' Given a root directory, recursively find all eQTL_result.txt files under that directory
#' For each, extract the empiric and permuted results for the SNP/Gene pairs of interest
#' as well as the best results for each SNP or gene independently
#'
#' @param rootDir A root directory under which to find eQTL files
#' @param queryFile A tab delimited file with 2 columns.  SNP: the rsID of the SNP of interest.
#' GENE: the gene of interest.  Each row corresponds to a SNP/Gene pair.  SNPs and genes may be
#' repeated multple times on different rows.
#' @param outFile The location to write the output file.  This file is fairly verbose,
#' with the explicit SNP/Gene pair and independent SNP and gene results on a single line.
#' @param outSummaryFile The location to write the summary file to.  This file contains the best (most significant
#' q-value or emperic p-value) for each SNP gene pair across all experiments, as well as the best
#' result for each SNP and gene independently.
#' @param ldTargetPopulation if supplied, use LDLinkR to add LD scores to the summary file for the "best gene" result.
#' @export
#' @import LDlinkR
aggregateEqtlResults<-function (rootDir, queryFile, outFile, outSummaryFile, ldTargetPopulation=NULL) {
    #gather up the file names of all the eQTL result files under the root directory
    eQTLFiles=list.files(path=rootDir, pattern="eQTL_results.txt.gz", full.names = T, recursive = T)
    query=read.table(queryFile, header=T, stringsAsFactors = F, sep="\t")
    result=do.call(rbind, lapply(1:length(eQTLFiles), fetchOneResult, eQTLFiles, query))
    write.table(result, outFile, row.names=F, col.names = T, quote=F, sep="\t")
    if (!is.null(outSummaryFile))
      summarizeBestHits(outFile, outSummaryFile, ldTargetPopulation)

}

summarizeBestHits<-function (outFile, outSummaryFile, ldTargetPopulation=NULL) {
    getBestExact<-function (index, query, a) {
        z=a[a$explicit_pair_gene==q[index,]$gene & a$explicit_pair_rsID==q[index,]$rsID,]
        idxQ=which.min(z$explicit_pair_q_value)
        if (length(idxQ)>0) return (z[idxQ,])
        idxE=which.min(z$explicit_pair_pvalue)
        return (z[idxE,])
    }

    getBestGene<-function (gene,a ) {
        z=a[a$best_snp_gene==gene ,]
        idxQ=which.min(z$best_snp_q_value)
        if (length(idxQ)>0) return (z[idxQ,])
        idxE=which.min(z$best_snp_pvalue)
        return (z[idxE,])
    }

    getBestSNP<-function (snp, a) {
        z=a[a$best_gene_rsID==snp ,]
        idxQ=which.min(z$best_gene_q_value)
        if (length(idxQ)>0) return (z[idxQ,])
        idxE=which.min(z$best_gene_pvalue)
        return (z[idxE,])
    }

    a=read.table(outFile, header = T, stringsAsFactors = F, sep="\t")
    q=unique (a[,c("gene", "rsID")])

    #column name manipulation to make all the results have the same fields
    colnames(a)=sub("rsid", "rsID", colnames(a))
    #change the gene and rsID to explicitly say..explicit pair.
    idx=match(c("gene", "rsID"), colnames(a))
    if (length(idx)>0) colnames (a)[idx]=paste("explicit_pair", colnames(a[idx]), sep="_")

    colsAll=c("experiment", "cis_dist", "peer_factors")
    fields=c("gene", "rsID", "pvalue", "beta", "median_expression", "effect_size", "permuted_p", "q_value")

    colsExact=c(colsAll, paste("explicit_pair", fields, sep="_"))
    geneCols=c(colsAll, paste("best_snp", fields, sep="_"))
    snpCols=c(colsAll, paste("best_gene", fields, sep="_"))

    #Exact query results
    exact=do.call(rbind,lapply(1:dim (q)[1], getBestExact, query, a))
    exact=exact[,colsExact]
    colnames (exact)=c(colsAll, fields)
    exact=cbind(query="exact", exact, LD=NA)

    #the best SNP for each gene.
    geneResults=do.call(rbind, lapply(unique (q$gene), getBestGene, a))
    geneResults=geneResults[,geneCols]
    colnames (geneResults)=c(colsAll, fields)
    geneResults=cbind(query="best_snp_for_gene", geneResults)

    #the gene results have an additional column - the LD of the selected SNP to the best SNP.
    #there can be more than 1 SNP that is queried for a gene.
    geneResults=fetchLDInfo(geneResults, q, ldTargetPopulation)

    #the best gene for each SNP
    snpResults=do.call(rbind, lapply(unique (q$rsID), getBestSNP, a))
    snpResults=snpResults[,snpCols]
    colnames (snpResults)=c(colsAll, fields)
    snpResults=cbind(query="best_gene_for_snp", snpResults, LD=NA)

    all=rbind(exact, geneResults, snpResults)
    write.table(all, outSummaryFile, row.names=F, col.names = T, quote=F, sep="\t")
}


fetchLDInfo<-function (geneResults, q, ldTargetPopulation=NULL) {
    #this only needs to annotate results where we're allowing the SNP to change for the given query gene.

    getOneGeneR2<-function (index=1, geneResults, q, ldTargetPopulation) {
        if (is.null(ldTargetPopulation)) return (NA)
        cat (index)
        indexSNP=geneResults[index,]$rsID
        querySNPs=q[which(q$gene %in% geneResults[index,]$gene),]$rsID
        allSNPs=unique (c(indexSNP, querySNPs))
        notIndexSNPs=setdiff(querySNPs, indexSNP)
        r=LDlinkR::LDmatrix(allSNPs, pop = ldTargetPopulation, r2d = "r2", token = Sys.getenv("LDLINK_TOKEN"), file = FALSE)
        r2=r[which(r$RS_number==indexSNP),notIndexSNPs, drop=F]
        paste(colnames(r2), as.numeric(r2), sep=":", collapse=",")
    }

    r=sapply(1:dim(geneResults)[1], getOneGeneR2, geneResults, q, ldTargetPopulation)
    geneResults$LD=r
    return (geneResults)

}

#eQTLFile="/broad/mccarroll/dropulation/analysis/eQTL/BA46/BA46.raw.dge.astrocytes.All.noOutliers.peer0/maf_0.05_cisDist_10kb/BA46.raw.dge.astrocytes.All.noOutliers.peer0.maf_0.05_cisDist_10kb.eQTL_results.txt.gz"
fetchOneResult<-function (index, eQTLFiles, query) {
    eQTLFile=eQTLFiles[index]
    print (paste("Starting analysis [", index, "/", length(eQTLFiles),"] file [", eQTLFile, "]", sep=""))

    expName=sub(".eQTL_results.txt.gz", "", basename(eQTLFile))
    directory=dirname(eQTLFile)

    props=readProperties(expName, directory)

    #read the empiric results
    er=data.table::fread(eQTLFile, sep="\t", data.table=FALSE)
    if (length(grep ("rs", er$id))==0) {
        warning (paste("Data set does not contain rsIDs in the id field, skipping analysis [", eQTLFile,"]", sep=""))
        return ()
    }
    #gather all results.
    result=do.call(rbind, lapply(1:dim (query)[1], getEmpiricResult, query, er))

    #read in the permuted `p-value`s
    eigenMT=readEigenMTResults(expName, directory)

    #add the permuted results
    result2=do.call(rbind, lapply(1:dim(result)[1], addEigenMTResult, result, eigenMT))

    #add the property data.
    p=data.frame(experiment=props[props$key=="EXPERIMENT",]$value, cis_dist=props[props$key=="CIS_DIST",]$value, peer_factors=props[props$key=="EXTRACT_PEER_FACTORS",]$value)

    #final result
    fr=cbind (p, result2)
    return (fr)
}



getEmpiricResult<-function (index, query, er) {
    #query the explicit pair, then the best hit for the gene, then the best hit for the SNP

    #TODO: check what happens when the query returns no results, and act appropriately.
    #does getNAIfEmpty cover the behavior, or does there need to be more protection?

    #explicit test
    q=query[index,]
    z=er[er$id==q$SNP & er$gene==q$GENE,]

    g=er[er$gene==q$GENE,]
    g=g[order(g$`p-value`, decreasing = F),]

    s=er[er$id==q$SNP,]
    s=s[order(s$`p-value`, decreasing = F),]

    result=data.frame(gene=q$GENE, rsID=q$SNP, explicit_pair_pvalue=getNAIfEmpty("p-value", z),
                      explicit_pair_beta=getNAIfEmpty("beta", z),
                      explicit_pair_median_expression=getNAIfEmpty("median_expression", z),
                      explicit_pair_effect_size=getNAIfEmpty("effect_size", z),
                      best_snp_rsid=g[1,]$id, best_snp_gene=g[1,]$gene, best_snp_pvalue=getNAIfEmpty("p-value", g[1,]),
                      best_snp_beta=getNAIfEmpty("beta", g[1,]),
                      best_snp_median_expression=getNAIfEmpty("median_expression", g[1,]),
                      best_snp_effect_size=getNAIfEmpty("effect_size", g[1,]),
                      best_gene_rsid=s[1,]$id, best_gene_gene=s[1,]$gene, best_gene_pvalue=getNAIfEmpty("p-value", s[1,]),
                      best_gene_beta=getNAIfEmpty("beta", s[1,]),
                      best_gene_median_expression=getNAIfEmpty("median_expression", s[1,]),
                      best_gene_effect_size=getNAIfEmpty("effect_size", s[1,]),
                      stringsAsFactors=F
                      )

    return (result)

}


addEigenMTResult<-function (index, result, eigenMT) {
    q=result[index,]

    #gene=q$gene;snp=q$rsID
    getResult<-function (gene, snp, eigenMT, colName="gene_permuted_pvalue") {
        #if either parameter is empty don't bother to search
        if (is.na(gene) | is.na(snp)) return (NA)
        z=eigenMT[eigenMT$id==snp & eigenMT$gene==gene,]
        getNAIfEmpty(colName, z)
    }


    df=data.frame(explicit_pair_permuted_p=getResult(q$gene, q$rsID, eigenMT, "gene_permuted_pvalue"),
                  explicit_pair_q_value=getResult(q$gene, q$rsID, eigenMT, "qvalue"),
                  best_snp_permuted_p=getResult(q$best_snp_gene, q$best_snp_rsid, eigenMT, "gene_permuted_pvalue"),
                  best_snp_q_value=getResult(q$best_snp_gene, q$best_snp_rsid, eigenMT, "qvalue"),
                  best_gene_permuted_p=getResult(q$best_gene_gene, q$best_gene_rsid, eigenMT, "gene_permuted_pvalue"),
                  best_gene_q_value=getResult(q$best_gene_gene, q$best_gene_rsid, eigenMT, "qvalue"),
                  stringsAsFactors=F)

    result=cbind (q, df)
    return(result)

}


getNAIfEmpty<-function (columnName, dataFrameRow) {
    r=dataFrameRow[[columnName]]
    if (length(r)==0) return (NA)
    return (r)
}

readEigenMTResults<-function (expName, directory) {
    eigenFile=paste(directory, "/", expName, ".eQTL_eigenMT_results.txt",sep="")

    if (!file.exists(eigenFile)) {
        stop (paste("Can't find eigenMT file [", eigenFile, "]"))
    }

    a=read.table(eigenFile, header=T, stringsAsFactors = F, sep="\t")
    return (a)

}

readProperties<-function (expName, directory, propFile=NULL) {
    if (!missing(expName) && !missing(directory) && is.null(propFile)) {
      propFile=paste(directory, "/", expName, ".properties",sep="")
    }
    if (!file.exists(propFile)) {
        stop (paste("Can't find properties file [", propFile, "]"))
    }

    #parse the prop file
    getParam<-function (x) {
        xx=strsplit (x, "=", fixed=T)[[1]]
        data.frame(key=xx[1], value=xx[2], stringsAsFactors = F)
    }

    a=read.table(propFile, header=F, stringsAsFactors = F)$V1
    r=lapply(a, getParam)
    r=do.call(rbind, r)
    return (r)
}
