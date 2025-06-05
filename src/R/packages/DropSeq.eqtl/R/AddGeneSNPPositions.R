


#' Add distance from SNP to gene in eQTL data.
#'
#' In an eQTL data frame, calculate the distance from the SNP to the gene.  SNPs found 5' of the gene have a negative distance, SNPs found 3' have a positive distance.  SNPs inside the gene are distance 0.
#'
#' @param eQTLData data frame that has a SNP column for the SNP name and a gene column for the gene name
#' @param snps_location_file_name file in the standard eQTL format containing the columns snp,chr,pos that encode the snp name, chromosome, and position.
#' @param gene_location_file_name file in the standard eQTL format containing the columns geneid,chr,s1,s2 that encode the snp name, chromosome, and start position, end position.
#' @return the eQTLData data frame annotated with the gene start/end positon, snp position, and distance to gene.
#' @import data.table
#' @export

getDistanceToGene<-function (eQTLData, snps_location_file_name, gene_location_file_name) {
    validateSingleFileExists(gene_location_file_name)
    validateSingleFileExists(snps_location_file_name)

    geneLoc=fread(gene_location_file_name)
    snpLoc=fread(snps_location_file_name)

    d=addLocationAnnotationToSNPs(eQTLData, snpLoc)
    d=addLocationAnnotationToGenes(d, geneLoc)

    df=data.frame(distStart=abs(d$snp_pos-d$gene_start_pos), distEnd=abs(d$snp_pos-d$gene_end_pos))
    minDist=apply (df, 1, min)
    d$dist_to_gene=minDist
    idxInGene=which(d$snp_pos>=d$gene_start_pos & d$snp_pos<=d$gene_end_pos)
    if (length(idxInGene)>0) d[idxInGene,]$dist_to_gene=0
    idxNegative=which(d$snp_pos<d$gene_start_pos)
    if (length(idxNegative)>0) d[idxNegative,]$dist_to_gene=d[idxNegative,]$dist_to_gene*as.integer(-1)
    return (d)
}

#' Adds the SNP location to an eQTL result data frame.
#'
#' @param eQTLData a data frame that has a column "SNP" containing the SNP ID.
#' @param snpLoc Data frame of SNP locations in the standard eQTL format containing the columns snp,chr,pos that encode the snp name, chromosome, and position.
#'
#' @return The data frame with an additional column snp_pos
#' @export
addLocationAnnotationToSNPs<-function (eQTLData, snpLoc) {
    idx=match(eQTLData$SNP, snpLoc$snp)
    if (length(which(is.na(idx)))>0) cat ("SNPs in report that I don't have locations of.")
    eQTLData$snp_pos=snpLoc[idx,]$pos
    eQTLData$snp_end=snpLoc[idx,]$end
    eQTLData$id=snpLoc[idx,]$id
    return (eQTLData)
}

#' Adds the gene location to an eQTL result data frame.
#'
#' @param eQTLData a data frame that has a column "gene" containing the gene ID.
#' @param geneLoc Data frame of gene locations in the standard eQTL format containing the columns geneid,chr,s1,s2 that encode the snp name, chromosome, and start position, end position.
#'
#' @return The data frame with an additonal columns gene_start_pos and gene_end_pos
#' @export
addLocationAnnotationToGenes<-function (eQTLData, geneLoc) {
    idx=match(eQTLData$gene, geneLoc$geneid)
    if (length(which(is.na(idx)))>0) cat ("Genes in report that I don't have locations of.")
    eQTLData$gene_start_pos=geneLoc[idx,]$s1
    eQTLData$gene_end_pos=geneLoc[idx,]$s2
    return (eQTLData)
}
