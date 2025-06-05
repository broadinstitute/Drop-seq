

# Standard analysis of single phenotype results
# Manhatten plot [qqman::manhattan()]
# QQ Plot [qqman:qq]
# Region plots for significant SNPs [RACER package]
#

#library (data.table)
#library (qqman)
#install_github("oliviasabik/RACER")
#library(RACER)
#library (logger)

#gwasResult="/broad/mccarroll/nemesh/Levy_Perturbation/single_phenotype_GWAS/cholesterol/out/Neuron_CLOZ_phenotype_sum.gwas_results.txt"
#expName="Neuron CLOZ sum"

# gwasResult="/broad/mccarroll/nemesh/Levy_Perturbation/single_phenotype_GWAS/cholesterol/out/Neuron_CLOZ_minus_DMSO_24hr_phenotype_sum.gwas_results.txt"
# expName="Neuron [CLOZ - DMSO] sum"
# variantLocationFile="/broad/mccarroll/dropulation/analysis/eQTL/Levy_Perturbation/Levy_village_Neuron.permissive_44_donor_WGS/maf_0.20_cisDist_10kb/Levy_village_Neuron.permissive_44_donor_WGS.maf_0.20_cisDist_10kb.variant_locations.txt"
# genotypeFile="/broad/mccarroll/dropulation/analysis/eQTL/Levy_Perturbation/Levy_village_Neuron.permissive_44_donor_WGS/maf_0.20_cisDist_10kb/Levy_village_Neuron.permissive_44_donor_WGS.maf_0.20_cisDist_10kb.genotype_matrix.txt.gz"
# phenotypeFileName="/broad/mccarroll/nemesh/Levy_Perturbation/single_phenotype_GWAS/cholesterol/Neuron_CLOZ_minus_DMSO_24hr_phenotype_sum.txt"
# outPDF="/broad/mccarroll/nemesh/Levy_Perturbation/single_phenotype_GWAS/cholesterol/Neuron_CLOZ_minus_DMSO_24hr_phenotype_sum.gwas.pdf"
# pValueThreshold=1e-5; population="EUR"
#
#
# eQTLResult="/broad/mccarroll/dropulation/analysis/eQTL/Levy_Perturbation/Levy_village_Neuron.permissive_44_donor_WGS/maf_0.20_cisDist_10kb/Levy_village_Neuron.permissive_44_donor_WGS.maf_0.20_cisDist_10kb.eQTL_results.txt.gz"
# variantLocationFile="/broad/mccarroll/dropulation/analysis/eQTL/Levy_Perturbation/Levy_village_Neuron.permissive_44_donor_WGS/maf_0.20_cisDist_10kb/Levy_village_Neuron.permissive_44_donor_WGS.maf_0.20_cisDist_10kb.variant_locations.txt"
#
# eQTLResult="/broad/mccarroll/dropulation/analysis/eQTL/Levy_Perturbation/Levy_village_Neuron_CLOZ.permissive_44_donor_WGS/maf_0.20_cisDist_10kb/Levy_village_Neuron_CLOZ.permissive_44_donor_WGS.maf_0.20_cisDist_10kb.eQTL_results.txt.gz"
# variantLocationFile="/broad/mccarroll/dropulation/analysis/eQTL/Levy_Perturbation/Levy_village_Neuron_CLOZ.permissive_44_donor_WGS/maf_0.20_cisDist_10kb/Levy_village_Neuron_CLOZ.permissive_44_donor_WGS.maf_0.20_cisDist_10kb.variant_locations.txt"

#DropSeq.eqtl::plotGenotypesVsExpression

#' Standard Analysis for a single phenotype GWAS scan
#'
#' Generates manhatten plot for empiric and permuted p-values
#' Generates a region plot for each locus that contains a SNP above the pvalue threshold
#' Generates a genotype vs phenotype plot for the best SNP in the locus
#'
#' @param gwasResult A file containing the permuted output of the testSinglePhenotype function
#' @param variantLocationFile The location of each SNP.  Encodes the snp name, chromosome, position, rsID.
#' @param population If supplied, fetches LD for the SNP from this thousand genomes population. See the LD link website: https://ldlink.nci.nih.gov/
#' @param pValueThreshold Pvalue threshold to define which loci are plotted.
#'
#' @export
#'
#' @seealso DropSeq.eqtl::testSinglePhenotype
#' @import data.table qqman RACER logger
singlePhenotypeStandardAnalysis<-function (gwasResult, variantLocationFile, population="EUR", pValueThreshold=1e-5) {
    d=readGWASInput(gwasResult, variantLocationFile)
    #qqman::qq(pvector = d$P)
    pdf(outPDF)
    par(mfrow=c(2,1), mar=c(4,4,3,2))
    manhattanPlot(d, pvalueColumn="empiric_p", strTitle=paste(expName, "\nempiric p-values"))
    manhattanPlot(d, pvalueColumn="P", strTitle=paste(expName, "\npermuted p-values"))
    par(mfrow=c(1,1), mar=c(4.5,4.5,4.5,4.5))

    allLeadSNPs = plotRegions(d, population, pValueThreshold)
    plotManyGenotypeVsPhenotype (allLeadSNPs, d, phenotypeFileName, genotypeFile, expName)

    dev.off()

}


#for an eQTL result, was the GWAS SNP or the GWAS SNP haplotype tested?
geteQTLResultForGwasSnp<-function (eQTLResult, variantLocationFile, population="EUR", leadSnpRSID="rs10182220", chr=2, pos=152145280, windowSize=250000) {
    d=readEQTLInput(eQTLResult, variantLocationFile)
    #restrict data to the window around the SNP of interest.
    dd=d[d$CHR==chr & d$POS>(pos-windowSize) & d$POS<(pos+windowSize),]
    genes=unique (dd$gene)

    plotOneGene<-function (geneName, dd, population="EUR", leadSnpRSID="rs10182220", chr=2, pos=152145280, windowSize=250000) {
        z=dd[dd$gene==geneName,]
        #extend the window size to encompass the gene being tested.
        #this plotting code is more and more lame.
        startDistance=abs(min (z$gene_start_pos-pos))
        endDistance=abs(max (z$gene_end_pos-pos))
        windowSizeFinal=windowSize
        if (startDistance>windowSizeFinal) windowSizeFinal=startDistance
        if (endDistance>windowSizeFinal) windowSizeFinal=endDistance
        plotRegionArbitrarySNP(z, population, leadSnpRSID, chr, pos, windowSizeFinal)
    }

    sapply(genes, plotOneGene, dd, population, leadSnpRSID, chr, pos, windowSize)



}



manhattanPlot<-function (d, pvalueColumn="P", strTitle="") {
    maxY=max (max (-log10(d[[pvalueColumn]])), 9)
    qqman::manhattan(d, chr="CHR", bp="POS", p=pvalueColumn, snp="RS_ID", ylim=c(0, maxY),
                     main=strTitle, col=c("darkslateblue","cadetblue"))
}


plotRegions<-function (d, population="EUR", pValueThreshold=1e-5) {
    #select candidate SNPs of interest.
    candidateSNPs=d[d$P<pValueThreshold,]
    candidateSNPs=candidateSNPs[order(candidateSNPs$P),]

    #the plots are 50kb around the SNP of interest, so as plots are generated remove any
    #overlapping candidates.
    allLeadSNPs=c()
    while (dim (candidateSNPs)[1]>0) {
        leadSnpRSID=candidateSNPs[1,]$RS_ID
        allLeadSNPs=c(allLeadSNPs, leadSnpRSID)

        r=plotRegion(d, population, leadSnpRSID=leadSnpRSID)
        plottedSNPs=candidateSNPs[candidateSNPs$CHR==r$chr & candidateSNPs$POS>=r$start & candidateSNPs$POS<=r$end,]$SNP
        idx=match(plottedSNPs, candidateSNPs$SNP)
        candidateSNPs=candidateSNPs[-idx,]
    }
    return (allLeadSNPs)

}

#is there any effect of the synaptic genes on the neurons for the CLOZ treatment?

#wrapper around RACER::ldRACER and RACER::singlePlotRACER.
#The start/end positions of the plots are hard coded,
plotRegion<-function (d, population="EUR", leadSnpRSID="rs10182220") {
    leadSnpRSIDFinal=strsplit (leadSnpRSID, ";")[[1]][1]
    #add LD info.  Grabs roughly 1 MB window of LD, with the lead SNP centered.
    dd = tryCatch({
        dd=RACER::ldRACER(d, rs_col=getColumnIndex("RS_ID", d), lead_snp=leadSnpRSIDFinal, pops=population)
    }, error = function(e) {
        log_warn("Unable to fetch LD for [",leadSnpRSID, "]")
        dd=d
    })

    chr=unique(dd[dd$RS_ID==leadSnpRSID,]$CHR)
    z=RACER::singlePlotRACER(dd, chr=chr, build="hg38", set="protein_coding", plotby="snp", snp_plot=leadSnpRSID, label_lead=T)
    print (z)
    #this is hard coded in singlePlotRacer:
    start=dd[dd$RS_ID==leadSnpRSID,]$POS-5e+05
    end=dd[dd$RS_ID==leadSnpRSID,]$POS+5e+05
    return (data.frame(chr=chr, start=start, end=end))

}

# plotRegionArbitrarySNP(d, population="EUR", leadSnpRSID="rs61896127", chr=11, pos=47078316, windowSize=250000)
# plotRegionArbitrarySNP(d, population="EUR", leadSnpRSID="rs2387414", chr=19, pos=50530986, windowSize=250000)
# plotRegionArbitrarySNP(d, population="EUR", leadSnpRSID="rs11656775", chr=17, pos=17751005, windowSize=250000)
# plotRegionArbitrarySNP(d, population="EUR", leadSnpRSID="rs5751191 ", chr=22, pos=41974987, windowSize=250000)
plotRegionArbitrarySNP<-function (d, population="EUR", leadSnpRSID="rs2387414", chr=19, pos=50530986, windowSize=250000) {

    #head (d[d$CHR==chr & d$POS>(pos-windowSize) & d$POS<(pos+windowSize),])

    #if the SNP isn't present in the tested data add it.
    #maybe you don't need this?
    dd=d
    if (!leadSnpRSID %in% d$RS_ID) {
        pad=d[1,]
        pad[!is.na(pad)]<-NA
        pad[c("CHR", "POS", "SNP", "P", "RS_ID", "LOG10P")]<-c(chr, pos, leadSnpRSID, 1, leadSnpRSID, 0)
        #pad=data.frame(CHR=chr, POS=pos, SNP=leadSnpRSID, empiric_p=1, excceds_empiric_p=NA, total_permutations=NA, P=1, RS_ID=leadSnpRSID, LOG10P=0)
        dd=rbind(d, pad)
    }

    leadSnpRSIDFinal=strsplit (leadSnpRSID, ";")[[1]][1]
    #plotRegion(z, population, leadSnpRSID)

    #add LD info.  Grabs roughly 1 MB window of LD, with the lead SNP centered.
    dd = tryCatch({
        dd=RACER::ldRACER(dd, rs_col=getColumnIndex("RS_ID", z), lead_snp=leadSnpRSIDFinal, pops=population)
    }, error = function(e) {
        log_warn("Unable to fetch LD for [",leadSnpRSID, "]")
        dd=dd
    })

    outPlot=RACER::singlePlotRACER(dd, chr=chr, build="hg38", set="protein_coding", plotby="coord", start_plot=(pos-windowSize), end_plot=(pos+windowSize), label_lead=T)
    print (outPlot)


}



plotManyGenotypeVsPhenotype<-function (allLeadSNPs, d, phenotypeFileName, genotypeFile, expName) {
    p=read.table(phenotypeFileName, header=T, stringsAsFactors=F, sep="\t", check.names = F)
    genotypes=DropSeq.utilities::fastRead(genotypeFile, check.names=F)
    p=DropSeq.eqtl::validatePhenotypeGenotypeInputs(genotypes, p)
    z=sapply(allLeadSNPs, plotGenotypesVsPhenotype, d, p, genotypes, expName)
}

# leadSnpRSID="rs10182220"; phenotypeValues=p

plotGenotypesVsPhenotype<-function (leadSnpRSID, d, phenotypeValues, genotypes, expName="") {

    indexSNP=d[d$RS_ID==leadSnpRSID,]$SNP
    permuted_pvalue=d[d$RS_ID==leadSnpRSID,]$P
    g=as.numeric(genotypes[genotypes$id==indexSNP,-1])
    e=as.numeric(phenotypeValues[-1])
    refAllele=strsplit (indexSNP, ":")[[1]][3]
    altAllele=strsplit (indexSNP, ":")[[1]][4]
    #if the SNP is missing from the data set, don't run the regression, and generate an empty plot.
    r=lm(e ~ g)
    r2=round (summary(r)$adj.r.squared,3)

    #set up the complicated title string.
    strTitle=paste(expName, leadSnpRSID, sep=" ")

    #Add the permuted pvalue
    strTitle=paste(strTitle, "\n perm p-value", format(permuted_pvalue ,digits=3, scientific=T))

    par(mar=c(5,5,6,1))
    plot (e ~ g, col="black", main=strTitle, ylab="Phenotype", xlab="GENOTYPE", xlim=c(-0.2,2.2), axes=F, cex.lab=1.5, cex.main=1.5)
    #axis(1, at=c(0,1,2), labels=c("RR", "RA", "AA"), cex.lab=2, cex.axis=1.5)
    axis(1, at=c(0,1,2), labels=c(paste(refAllele,refAllele,sep=""), paste(refAllele,altAllele,sep=""), paste(altAllele,altAllele,sep="")), cex.lab=2, cex.axis=1.5)
    axis(2, cex.lab=2, cex.axis=1.5)
    if (!is.null(r)) {
        yStart=r$coefficients[[1]]
        yEnd=2*r$coefficients[[2]]+yStart
        lines (x=c(0,2), y=c(yStart, yEnd), lty=2, col="black")
    }
    par(mar=c(5, 4, 4, 2) + 0.1)
}

#' Read the GWAS result from permutation testing and annotate / clean up the input for analysis
#'
readGWASInput<-function (gwasResult, variantLocationFile) {
    a=data.table::fread(gwasResult)
    v=data.table::fread(variantLocationFile)

    idx=match(a$SNP, v$snp)
    a$rs_id=v[idx,]$id

    #if there are multiple SNPs at the same position, pick the most significant one.
    #if there are multiple SNPs that share the same RS_ID, this screws up the 3rd party code
    #that relies on the identifier being unique.
    getDuplicatedSNPID<-function (rsID, dupesDF) {
        x=dupesDF[dupesDF$rs_id==rsID,]
        idx=which.min(x$snp_permuted_pvalue)
        x[-idx,]$SNP
    }

    idx=which(duplicated(a$rs_id))
    dupes=a[idx,]$rs_id
    dupesDF=a[a$rs_id %in% dupes,]
    dupesSNPIDs=unique(as.character(unlist(lapply(dupes, getDuplicatedSNPID, dupesDF))))

    #remove the dupes
    idxRemove=match(dupesSNPIDs, a$SNP)
    log_info(paste("Removing [", length(idxRemove), "] SNPs that have duplicate rsIDs of total [", dim (a)[1], "]"))
    if (length(idxRemove)>0)
        a=a[-idxRemove,]

    #rename columns
    setnames(a, "P", "empiric_p")


    #format with RACER
    aa=formatRACER(a, chr_col=getColumnIndex("CHR" ,a),
                  pos_col=getColumnIndex("BP", a),
                  p_col=getColumnIndex("snp_permuted_pvalue", a),
                  rs_col=getColumnIndex("rs_id", a)
                  )
    return (aa)

}

readEQTLInput<-function (eQTLResult, variantLocationFile) {
    a=data.table::fread(eQTLResult)
    v=data.table::fread(variantLocationFile)
    idx=match(a$SNP, v$snp)
    a$rs_id=v[idx,]$id

    #extract the chr/pos

    a$CHR=sapply(strsplit (a$SNP, ":", fixed=T), function (x) x[1])
    a$CHR=sub("chr", "", a$CHR)
    a$BP=sapply(strsplit (a$SNP, ":", fixed=T), function (x) as.numeric(x[2]))
    a[a$CHR=="X",]$CHR<-23

    aa=formatRACER(a, chr_col=getColumnIndex("CHR" ,a),
                   pos_col=getColumnIndex("BP", a),
                   p_col=getColumnIndex("p-value", a),
                   rs_col=getColumnIndex("rs_id", a)
    )
    return (aa)
}

#RACER requires the index to be numeric, not integer.  Lame.
getColumnIndex<-function (colName, a) as.numeric (which(colnames (a)==colName))

