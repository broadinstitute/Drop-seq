#help with standardized file names

getSNPLocationFileName<-function (base.dir, base.name) paste(base.dir, "/", base.name, ".variant_locations.txt", sep="");
getSNPFileName<-function (base.dir, base.name) paste(base.dir, "/", base.name , ".genotype_matrix.txt.gz", sep="");
getExpressionFileName<-function (base.dir, base.name) paste(base.dir, "/", base.name, ".gene_expression.txt", sep="");
getGeneLocationFileName<-function (base.dir, base.name) paste(base.dir, "/", base.name, ".gene_locations.txt", sep="");
getCovariateFileName<-function (base.dir, base.name) paste(base.dir,  "/", base.name, ".covariates.txt", sep="");
getIndexSNPsFileName<-function (base.dir, base.name) paste(base.dir, "/", base.name, ".eigenMT.index_snps.txt", sep="");
getEQTLResults<-function (base.dir, base.name) paste(base.dir, "/", base.name, ".eQTL_results.txt.gz", sep="");

getOutCiseQTLFileName<-function (base.dir, base.name, useQuantileNormalization, cisDist) {
    if (useQuantileNormalization) {
        output_file_name_cis = paste(base.dir, "/", base.name, ".eQTL_results_", (cisDist/1000), "kb.quantile_normalized.txt", sep="");
    } else {
        output_file_name_cis = paste(base.dir, "/", base.name, ".eQTL_results_", (cisDist/1000), "kb.txt", sep="");
    }
    return (output_file_name_cis)
}

#same as getOutCiseQTLFileName, but removes the .txt from the end.
getOutputBaseName<-function (base.dir, base.name, output.dir=base.dir, chromosome=NULL, useQuantileNormalization, cisDist) {
    eQTLInputFile=getOutCiseQTLFileName(output.dir, base.name, useQuantileNormalization, cisDist)
    eQTLInputFile=sub (".txt", "", eQTLInputFile)
    return (eQTLInputFile)
}

getBestSNPPermutationFileName<-function (base.dir, base.name, output.dir=base.dir, chromosome=NULL, useQuantileNormalization, cisDist) {
    eQTLInputFile=getOutputBaseName(base.dir, base.name, output.dir, chromosome, useQuantileNormalization, cisDist)
    r=paste(eQTLInputFile, ".permuted", sep="")
    if (!is.null(chromosome)) r=paste(r, ".chr_", chromosome, sep="")
    r=paste(r, ".txt", sep="")
    return (r)
}

getFullPermutationFileName<-function (base.dir, base.name, output.dir=base.dir, chromosome=NULL, useQuantileNormalization, cisDist) {
    eQTLInputFile=getOutputBaseName(base.dir, base.name, output.dir, chromosome, useQuantileNormalization, cisDist)
    r=paste(eQTLInputFile, ".full_permuted", sep="")
    if (!is.null(chromosome)) r=paste(r, ".chr_", chromosome, sep="")
    r=paste(r, ".txt", sep="")
    return (r)
}

getPermutationOrderFileName<-function (base.dir, base.name, output.dir=base.dir, chromosome=NULL, useQuantileNormalization, cisDist) {
    eQTLInputFile=getOutputBaseName(base.dir, base.name, output.dir, chromosome, useQuantileNormalization, cisDist)
    r=paste(eQTLInputFile, ".permutation_order", sep="");
    if (!is.null(chromosome)) r=paste(r, ".chr_", chromosome, sep="")
    r=paste(r, ".txt", sep="")
    return (r)
}

getSNPLevelReportFile<-function (base.dir, base.name, output.dir=base.dir, chromosome=NULL, useQuantileNormalization, cisDist) {
    eQTLInputFile=getOutputBaseName(base.dir, base.name, output.dir, chromosome, useQuantileNormalization, cisDist)
    r=paste(eQTLInputFile, ".snp_report", sep="");
    if (!is.null(chromosome)) r=paste(r, ".chr_", chromosome, sep="")
    r=paste(r, ".txt", sep="")
    return (r)
}


generatePermutationOrderFile<-function (numPermuations=10000, expressionFileName, permutationOrderFileName, random.seed=1) {
    set.seed(random.seed)
    a=fread(expressionFileName, nrows = 1)
    idx=1:(dim(a)[2]-1)

    getOne<-function (idx) {
        x=sample(idx, length(idx), replace = F)
        return (x)
    }
    r=t(replicate(numPermuations, getOne(idx)))
    #rr=apply(r, 1, function (x) paste(x, collapse=":"))
    write.table(r, permutationOrderFileName, row.names=F, col.names=F, quote=F, sep="\t")
}

#' Convert a data table into matrix eQTL's slice data format.
#'
#' @param dt a data table in the expression data or gentotype data format.
#' @param numRowsPerSlice the number of rows per slice.
#' @return A SlicedData object containing the data
#' @export
convertToSlicedData<-function (dt, numRowsPerSlice=1000) {
    #convert first to matrix.
    rowNames=dt$id
    m=as.matrix(dt[,-1])
    rownames(m)=rowNames
    sd = SlicedData$new();
    sd$CreateFromMatrix(m)
    sd$ResliceCombined(sliceSize = numRowsPerSlice)
    return (sd)
}

#' Compute the minor allele frequency for a set of SNPs
#'
#' Calculates the frequency of the less common allele per SNP.
#'
#' @param snps A set of SNP genotypes. SNPs can either be in a SlicedData object or a data.frame.
#'
#' @return A vector of allele frequencies, names set to the SNP IDs.
#' @export
getAlleleFrequency<-function (snps) {

    if ("SlicedData" %in% class(snps)) return (getAlleleFrequencyFromSlicedData(snps))
    if ("data.frame" %in% class(snps)) return (getAlleleFrequencyFromDF(snps))
}
getAlleleFrequencyFromSlicedData<-function (snps) {
    maf.list = vector('list', length(snps))
    for(sl in 1:length(snps)) {
        slice = snps[[sl]];
        maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
        maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
    }
    maf = unlist(maf.list)
    return (maf)
}

getAlleleFrequencyFromDF<-function (snps) {
    getAF<-function (x) {
        r=sum(x, na.rm=T)/(length(which(!is.na(x)))*2)
        min(r, 1-r)
    }
    mafs=apply(snps[,-1], 1, getAF)
    names (mafs)=snps$id
    return (mafs)
}

getMedianExpression<-function (gene) {
    if ("SlicedData" %in% class(gene)) return (getMedianExpressionFromSlicedData(gene))
    if ("data.frame" %in% class(gene)) return (getMedianExpressionFromDF(gene))
}

getMedianExpressionFromSlicedData<-function (gene) {
    mExp.list = vector('list', length(gene))
    for(sl in 1:length(gene)) {
        slice = gene[[sl]];
        mExp.list[[sl]]=apply (slice, 1, median)
    }
    r = unlist(mExp.list)
    df=data.frame(gene=names(r), medianExpression=as.numeric(r),stringsAsFactors = F)
    return (df)
}

getMedianExpressionFromDF<-function (gene) {
    df=data.frame(gene=gene$id, medianExpression=as.numeric (apply(gene[,-1], 1, median)),stringsAsFactors = F)
    return (df)

}


#' Add median expression per gene and effect size normalized to expression
#'
#' The effect size is the beta of the linear regression / median expression for that gene.
#' @param result The eQTL results
#' @param gene The gene expression data
#'
#' @return The eQTL results with added columns.
#' @export
addMedianExpression<-function (result, gene) {
    mExp=getMedianExpression(gene)
    idx=match(result$gene, mExp$gene)
    result$median_expression=mExp[idx,]$medianExpression
    if ("beta" %in% colnames(result))
        result$effect_size=result$beta/result$median_expression
    return (result)
}

#' Convert SNP names so they are gene-specific
#' @param snps vector of SNP names
#' @param genes vector of gene names
#' @return vector of SNP-name:gene-name
makePerGeneSnpName<-function(snps, genes) {
    return(paste(snps, genes, sep=":"))
}

# Expand SNP locations so that there is a distinct row for every gene for which it will be considered,
# with contig name == gene to be considered, and SNP name augmented with gene name.
snpLocationsToPerGeneCoordinates<-function(snpspos, snp_gene_map) {
    # Create a snpspos row for each {snp,gene} pair in snp_gene_map
    snpspos = merge(snpspos, snp_gene_map, by="snp")
    # Make the contig for each snpspos row be the same as the gene
    snpspos$chr = snpspos$gene
    # Rename each copy of the SNP based on the gene to be co-analyzed
    snpspos$snp = makePerGeneSnpName(snpspos$snp, snpspos$chr)
    return(snpspos)
}

# Convert gene locations so that each gene is on a separate contig, which has the same name as the gene.
geneLocationsToPerGeneCoordinates<-function(genepos, snp_gene_map) {
    # make the contig be the same as the gene
    genepos$chr = genepos$geneid
    return(genepos)
}

# Expand genotypes so that there is a distinct for for every SNP, gene for which the SNP is to be considered.
# with SNP name augmented with gene name.
snpsToPerGeneCoordinates<-function(snps, snp_gene_map) {
    cols = colnames(snps)
    # Create a separate row in genotypes for each SNP according to the
    # gene to be co-analyzed.
    snps = merge(snps, snp_gene_map, by.x="id", by.y="snp")
    # Rename each copy of the SNP based on the gene to be co-analyzed
    snps$id = makePerGeneSnpName(snps$id, snps$gene)
    # drop the gene row
    snps = snps[, match(cols, colnames(snps)), with=FALSE]
    return(snps)
}

# Reverse the transformation above in which :gene name is appended to SNP name
undoPerGeneSnpName1 <- function(snp) {
    fields = strsplit(snp, ':', fixed = TRUE)[[1]]
    # drop last element (the gene name)
    length(fields) = length(fields) - 1
    return(paste(fields, collapse=":"))
}

undoPerGeneSnpName <- function(snpVec) {
    return(vapply(snpVec, undoPerGeneSnpName1, FUN.VALUE = "X"))

}

#' Load genotype data for eQTL analysis
#'
#' @param SNP_file_name genotypes => row: SNP; column: sample; value: {0,1,2,NA}
#' @param snps_location_file_name columns: snp, chr, pos, id
#' @param snp_gene_map columns: snp, gene.  If present, this is joined to genotype data and
#'                               SNP location data.  Both genotype data and location data for a SNP are
#'                               duplicated for each gene that is to be considered with the SNP.  The SNP
#'                               name has the gene name appended to it to make MatrixEQTL work properly.
#' @param minorAlleleFreqThreshold Filter variants by minor allele frequency
#' @return list(snps, snpspos)
loadGenotypeData<-function (SNP_file_name, snps_location_file_name,
                            minorAlleleFreqThreshold=NULL, snp_gene_map=NULL) {
    snps=DropSeq.utilities::fastRead(SNP_file_name)
    snpspos = fread(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);

    if (!is.null(snp_gene_map)) {
        snpspos = snpLocationsToPerGeneCoordinates(snpspos, snp_gene_map)
        snps = snpsToPerGeneCoordinates(snps, snp_gene_map)
    }

    #covert to data slices.
    snps=convertToSlicedData(snps)

    #convert snpspos to data frame.
    setDF(snpspos)

    ## Load genotype data
    # this is the slow-ass way to do it.
    # cat ("Loading genotype Data\n")
    # snps = SlicedData$new();
    # snps$fileDelimiter = "\t";      # the TAB character
    # snps$fileOmitCharacters = "NA"; # denote missing values;
    # snps$fileSkipRows = 1;          # one row of column labels
    # snps$fileSkipColumns = 1;       # one column of row labels
    # snps$fileSliceSize = 10000;      # read file in slices of 2,000 rows
    # snps$LoadFile(SNP_file_name);

    #filter SNPs by minor allele frequency
    r=list(snps=snps, snpspos=snpspos)
    if (is.null(minorAlleleFreqThreshold)) return (r)

    maf.list = vector('list', length(snps))
    for(sl in 1:length(snps)) {
        slice = snps[[sl]];
        maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
        maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
    }
    maf = unlist(maf.list)

    ## Look at the distribution of MAF
    #hist (maf, breaks=50)
    #hist(maf,seq(min(maf),0.1,0.01))
    #hist(maf,seq(min(maf),max(maf),0.01))
    cat('SNPs before filtering:',nrow(snps))
    snps$RowReorder(maf>minorAlleleFreqThreshold);
    cat('SNPs after filtering:',nrow(snps))
    r=list(snps=snps, snpspos=snpspos)
    return (r)
}

#' Load gene expression data
#'
#' Load gene expression data
#'
#' @param expression_file_name The expression file
#' @param useQuantileNormalization Should the data be quantile normalized
#'
#' @return A SlicedData object holding the expression data.
#' @import MatrixEQTL stats
#' @export
loadExpressionData<-function (expression_file_name, useQuantileNormalization=F) {
    gene=fread(expression_file_name)
    #covert to data slices.
    gene=convertToSlicedData(gene)

    #perform quantile normalization if requested
    if (useQuantileNormalization) {
        for( sl in 1:length(gene) ) {
            mat = gene[[sl]];
            mat = t(apply(mat, 1, rank, ties.method = "average"));
            mat = qnorm(mat / (ncol(gene)+1));
            gene[[sl]] = mat;
        }
        rm(sl, mat);
    }
    return (gene)
}



## Load covariates
loadCovariates<-function (covariates_file_name=NULL) {

    cvrt = SlicedData$new();
    cvrt$fileDelimiter = "\t";      # the TAB character
    cvrt$fileOmitCharacters = "NA"; # denote missing values;
    cvrt$fileSkipRows = 1;          # one row of column labels
    cvrt$fileSkipColumns = 1;       # one column of row labels

    if (is.null(covariates_file_name)) return (cvrt)

    a=read.table(covariates_file_name, header=T, stringsAsFactors=F) #read in the file in case its empty.
    if (dim (a)[1]==0) return (cvrt) #this is how matrixEQTL likes it.

    if(length(covariates_file_name)>0) {
        cvrt$LoadFile(covariates_file_name);
    }
    return (cvrt)
}


validateFilesExist<-function (snps_location_file_name, SNP_file_name, expression_file_name, gene_location_file_name, covariates_file_name) {

    validateSingleFileExists(snps_location_file_name)
    validateSingleFileExists(SNP_file_name)
    validateSingleFileExists(expression_file_name)
    validateSingleFileExists(gene_location_file_name)
    validateSingleFileExists(covariates_file_name)
}


#' Check if a file exists, and if not print a nice error message and stop.
#'
#' @param x The file to check.  If NULL this is a no-op.
#' @export
validateSingleFileExists<-function (x) {
    vName=deparse(substitute(x))
    if (!is.null(x)) {
        if (!file.exists (x)) {
            print (paste ("Expected file ",x , sep=""))
            stop((paste(vName, "is not null and is missing")))
        }
    }
}
