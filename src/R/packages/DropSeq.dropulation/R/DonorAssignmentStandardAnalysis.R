# MIT License
#
# Copyright 2018 Broad Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


#This is for R CMD CHECK
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))


plotCellsByProbabilityNew<-function(df, fdrThreshold=0.05) {
	#TO make R CMD CHECK happy.  Blame ggplot2.
	FDR_pvalue=NULL
	
	plotOne<-function (a, strTitle) {
		p=ggplot(a, aes(x = log10(FDR_pvalue))) +
			geom_histogram(bins = 100, fill = "skyblue", color = "black") +
			labs(title = strTitle, x = "FDR p-value [log10]", y = "Frequency") +
			theme(plot.title = element_text(size = 14), 
				  axis.text = element_text(size = 10), 
				  axis.title = element_text(size = 12))
		print (p)
	}
	
	
    totalNumCells=dim(df)[1]
    strTitle=paste("All FDR corrected pvalues\n", round(length(which(df$FDR_pvalue< fdrThreshold))/dim(df)[1]*100,1), "% of cells with FDR < ", fdrThreshold)
    plotOne(df, strTitle)
    
    #remove doublets.
    cellsToKeep=df[df$label_simple=="singlet",]$cell_barcode
    
    if (!is.null(cellsToKeep) & length(cellsToKeep)>0 ) {
        idx=sort(match(cellsToKeep, df$cell_barcode))
        if (length(idx)>0) {
        	a=df[idx,]
        	cellsPassFDR=length(which(a$label=="single_pass_FDR"))
        	pct=paste(round(cellsPassFDR/totalNumCells*100,1), "% of cells", sep="")
        	strTitle=paste("All FDR corrected pvalues after filtering doublets\n", pct, "[", cellsPassFDR, "/", dim(a)[1], "] with FDR < ", fdrThreshold)
        	plotOne(a, strTitle)
        	
        }
    }
}

validateFilesExist<-function (likelihoodSummaryFile, doubletLikelihoodFile, dgeRawSummaryFile, dgeSummaryFile, censusFile, expectedSamplesFile) {
    f<-function (x) {
        vName=deparse(substitute(x))
        if (!is.null(x)) {
            if (!file.exists (x)) {
                stop((paste(vName, "is not null and is missing")))
            }
        }
    }
    f(likelihoodSummaryFile)
    f(doubletLikelihoodFile)
    f(dgeRawSummaryFile)
    f(dgeSummaryFile)
    f(censusFile)
    f(expectedSamplesFile)
}

#calculate the probability of being a singlet.
#returns cells unlikely to be a singlet (AKA some doublet.)
plotDoubletProbability<-function (df, expName="", summaryStats, doubletPvalueThreshold=0.9) {
	#TO make R CMD CHECK happy.  Blame ggplot2.
	doublet_pval=NULL
	
	idx=which(df$label_simple=="singlet")
	
    numSinglets=length(idx)
    pctSinglets=round (numSinglets/dim (df)[1]*100,1)
    
    # Create the histogram
    
    # There's a error with the new ggplot2 code (V 3.5) where this produces a warning if you use geom_segment directly.
    # This is similar to: https://github.com/tidyverse/ggplot2/issues/5762
    # Using hacky fix for now.
    p=ggplot(df, aes(x = doublet_pval)) +
    	geom_histogram(bins = 100, fill = "skyblue", color = "black", show.legend = FALSE) +
    	labs(x = "Doublet Probability", y = "Frequency") +
    	ggtitle(paste(expName, "\nCells remaining [", numSinglets, "] [", pctSinglets, "%]")) +
    	annotate("segment", x =doubletPvalueThreshold , y = 0, xend = doubletPvalueThreshold, yend = Inf, colour = "red") + 
    	annotate("text", x = 0.5, y = Inf, label = paste("Confident doublet rate", round(summaryStats$pct_confident_doublets, 2), "%"),hjust = 0.5, vjust = 2, size = 7) +
    	theme(plot.title = element_text(size = 14, hjust = 0.5, face="bold"), axis.text = element_text(size = 12), axis.title = element_text(size = 12))
    
    print (p)
}

plotCommonDonors<-function (df, minimumFractionDonor=0.005) {
	#TO make R CMD CHECK happy.  Blame ggplot2.
	index=value.Freq=NULL
	
    plotOne <- function(r, cex.axis = 0.5, main) {
    	# Exit early if there's no data.
    	if (length(r) == 0) return(NULL)
    	
    	# Create a data frame from the input vector
    	df <- data.frame(index = seq_along(r), value = r)
    	
    	maxY=max (df$value.Freq)*1.1
    	
    	# Define plot aesthetics
    	p <- ggplot(df, aes(x = index, y = value.Freq)) +
    		geom_line(linetype="dashed") +   # Add lines
    		geom_point(size=4) + 
    		labs(x = "", y = "Number of cells assigned to donor") +
    		scale_x_continuous(breaks = seq_along(r), labels = names(r), minor_breaks = NULL) +
    		scale_y_continuous(expand = c(0, 0), limits = c(0, maxY)) +
    		ggtitle(main) +
    		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10), 
    			  axis.text.y = element_text(size = 10),
    			  axis.title = element_text(size = 12),
    			  plot.title = element_text(size = 14)) +
    		theme(panel.background = element_rect(colour = "black"))
    		
		print (p)
    }

    result=df[df$label=="single_pass_FDR",]
    rUnfiltered=sort(table(result$donor), decreasing=T)

    minNum=sum(rUnfiltered)* minimumFractionDonor
    rFiltered=rUnfiltered[rUnfiltered>minNum]

    plotOne(rUnfiltered, cex.axis=0.75, main=paste("All Donor Assignments seen at least once\n[FDR passing cells] [", sum(rUnfiltered), "]"))
    plotOne(rFiltered, cex.axis=0.75, main =paste("Common Single Donor Assignment\n[FDR passing cells] [", sum(rFiltered), "]"))
	
    df=data.frame(donor=names(rFiltered), count=as.numeric(rFiltered), stringsAsFactors=F)
    df2=result[,c("cell_barcode", "donor")]

    #filter df2 to be consistent with rFiltered
    df2=df2[df2$donor %in% names(rFiltered),]
    #change column names to expected
    colnames(df2)=c("cell", "bestSample")
    result=list(summary=df, cellDonorMap=df2)
    return (result)

}

plotCommonDonorAssignmentsWithUnexpected<-function (donors, es, expName) {
	#TO make R CMD CHECK happy.  Blame ggplot2.
	index=count=expected=NULL
	
	#Skip this plot if there's no data.
	if (dim (donors)[1]==0) return (NULL)
	
	idx=match(donors$donor, es$V1)
	donors$expected=T
	if (any(is.na(idx))) donors[is.na(idx),]$expected=F
	missing=setdiff(es$V1, donors$donor)
	if (length(missing)>0) {
		d2=data.frame(donor=missing, count=0, expected=T)
		d2=rbind(donors, d2)
	} else {
		d2=donors
	}
	
	plotOne<-function (d2, expName="") {
		d=vegan::diversity(d2$count)
		eq=d/log(dim(d2)[1])
		strTitle=paste (expName, "\n", "SW Div: ", sprintf("%.2f",round(d,2)), "; SW Eq: ", sprintf("%.2f", round (eq,3)), sep="")
		d2$index=1:dim(d2)[1]
		maxY=max(d2$count)*1.1
		p=ggplot(d2, aes(x = index, y = count, color = expected)) +
			geom_line(linetype = "dashed", color = "black") +   # Add lines
			geom_point(data = d2, aes(color = expected), size = 4) +
			labs(x = "", y = "Number of cells assigned to donor", color = "Expected") +
			scale_x_continuous(breaks = seq(from = 1, to = dim(d2)[1]), labels = d2$donor, minor_breaks = NULL) +
			scale_y_continuous(expand = c(0, 0), limits = c(0, maxY)) +
			ggtitle(strTitle) +
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10), 
				  axis.text.y = element_text(size = 10),
				  axis.title = element_text(size = 12),
				  plot.title = element_text(size = 14),
				  legend.position = c(0.8, 0.8)) +
				  # This is for ggplot 3.5, which we don't want to force yet.
#				  legend.position = "inside", legend.position.inside = c(0.8, 0.8)) +
			theme(panel.background = element_rect(colour = "black")) +
			scale_color_manual(values = c("black", "red"), 
							   breaks = c(TRUE, FALSE),
							   labels = c("EXPECTED", "UNEXPECTED"),
							   name="Donor Status")	
		print (p)	
	}
	
	plotOne(d2, expName)
	
}


#' Plot the single cell FDR distribution of each donor
#'
#' @param df The data frame generated by getSingletDoubletDF
#' @inheritParams donorAssignmentQC
#' @noRd
plotCommonDonorsFdrDistribution<-function (df, minimumFractionDonor=0.005) {
	#TO make R CMD CHECK happy.  Blame ggplot2.
	donor=FDR_pvalue=label=value=NULL
	
	# Subset the data
	rUnfiltered=sort(table(df[df$label_simple=="singlet",]$donor), decreasing=T)
	
	minNum=sum(rUnfiltered)* minimumFractionDonor
	rFiltered=rUnfiltered[rUnfiltered>minNum]
	commonDonors=names (rFiltered)
	
	subset_data=df[df$donor %in% commonDonors,]
	
	# Create the boxplot with -log10 transformation
	ggplot(subset_data, aes(x = donor, y = -log10(FDR_pvalue))) +
		geom_violin() +
		labs(x = "Donor", y = "-log10(FDR p-value)") +
		ggtitle("Singlet cell -log10(FDR p-values) distribution for each common donor")
	
}

plotCensusComparison<-function (censusFile, summary, expName) {

    summary$frac=summary$count/sum(summary$count)
    census=read.table(censusFile, header=T, stringsAsFactors=F)

    both=union (summary$donor, census$DONOR)
    summary=summary[match(both,summary$donor),]
    summary$donor=both
    idx=which(is.na(summary$frac))
    if (length(idx)>0) summary[idx,]$frac<-0

    census=census[match(both, census$DONOR),]
    census$DONOR=both
    idx=which(is.na(census$REPRESENTATION))
    if (length(idx)) census[idx,]$REPRESENTATION<-0

    idx=match(summary$donor, census$DONOR)
    census=census[idx,]

    df=data.frame(donor=summary$donor, dropulation=summary$frac, census=census$REPRESENTATION, stringsAsFactors = F)
    df$foldChange=exp(abs(log(df$census/df$dropulation)))
    foldChangeScore=exp (median (abs(log(df$census/df$dropulation))))
    lim=max(c(df$dropulation,df$census))

    par(mar=c(5, 4, 4, 2) + 0.1)
    plot (df$dropulation, df$census, xlim=c(0, lim), ylim=c(0, lim), xlab="dropulation", ylab="census", cex.axis=1.5, cex.lab=1.5, pch=16, cex=2)
    title(sub=paste("Median fold change: ", round(foldChangeScore, 3)))
    title(paste(expName, "Dropulation and Census Results"))
    abline(0,1, col='red', lty=2)
    return (df)
}

calculateReadsPerUMI<-function (readsPerCellFile=NULL, dgeSummaryFile=NULL, cellDonorMap) {
	#if there's no data, exit and return NA.
	if (dim(cellDonorMap)[1]==0)
		return (NA)
	
	if (is.null(readsPerCellFile) | is.null(dgeSummaryFile))
		return (NA)
	
    a=read.table(readsPerCellFile, header=F, stringsAsFactors = F)
    ds=read.table(dgeSummaryFile, header=T, stringsAsFactors = F)
    idx=match(cellDonorMap$cell, ds$CELL_BARCODE)
    ds=ds[idx,]
    idx=match(ds$CELL_BARCODE, a$V2)
    a=a[idx,]
    #plot (a$V1, ds$NUM_TRANSCRIPTS, xlab="reads", ylab="UMIs")
    readsPerUMI=a$V1/ds$NUM_TRANSCRIPTS
    readsPerUMI=readsPerUMI[!is.na(readsPerUMI)]
    medianReadsPerUMI=median(readsPerUMI)

    df=data.frame(readsPerUMI=readsPerUMI)
    p=ggplot(data=df, aes(y = readsPerUMI, x=0)) +
    	geom_boxplot(width=0.35) +
    	labs(
    		title = paste("Reads per UMI in assigned singleton cells\nMedian=", round(medianReadsPerUMI, 2)),
    		y = "Reads/UMIs"
    	) +
    	ylim(0, max(df$readsPerUMI)) +
    	xlim(-1,1) +
    	theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
    	
    print (p)	
    
    return (medianReadsPerUMI)

}

calculateMeanUMIsPerDonor<-function (expName, cellDF, dgeSummaryFile=NULL, donors, cellsToKeep=NULL) {
	#TO make R CMD CHECK happy.  Blame ggplot2.
	donor=medianUMIs=label=variable=totalUMIs=count=NULL
	
	#if there's no data quit early.
	if ((!is.null(cellsToKeep) & length(cellsToKeep)==0) | is.null(dgeSummaryFile)) {
		r=list(stats=data.frame(diversity=NA, equitability=NA), donors=donors)
		return (r)
	}
		
	df=data.table::copy(cellDF)
    if (!is.null(cellsToKeep)) {
        idx=sort(match(cellsToKeep, df$cell_barcode))
        if (length(idx)>0) df=df[idx,]
    }
	
	#filter to the final donor list (which may be filtered)
    idx=which(!is.na(match(df$donor, donors$donor)))
    df=df[idx,]
    
    #if there's no data, 
    dge=read.table(dgeSummaryFile, header=T, stringsAsFactors=F, check.names=F)
    
    getUMIsPerDonor<-function (donor, df, dge) {
        cellBarcodes=df[df$donor==donor,]$cell_barcode
        dgeDonor=dge[match(cellBarcodes, dge$CELL_BARCODE),]
        return (list(dgeDonor$NUM_TRANSCRIPTS))
    }

    umisPerDonor=sapply(donors$donor, getUMIsPerDonor, df, dge)
    medianUMisPerDonor=unlist(lapply(umisPerDonor, median, na.rm=T))
    medianAbsoluteDeviationPerDonor=unlist(lapply(umisPerDonor, mad, na.rm=T))
    totalUMIsPerDonor=unlist(lapply(umisPerDonor, sum, na.rm=T))

    sortedIndex=sort(medianUMisPerDonor, index.return=T)$ix

    umisPerDonor= umisPerDonor[sortedIndex]
    medianUMisPerDonor= medianUMisPerDonor[sortedIndex]
    totalUMIsPerDonor=totalUMIsPerDonor[sortedIndex]

    #add to dataframe.
    idx=match(donors$donor, names(medianUMisPerDonor))
    donors$medianUMIs=unlist(medianUMisPerDonor)[idx]
    idx=match(donors$donor, names(medianAbsoluteDeviationPerDonor))
    donors$medianAbsoluteDeviation=unlist(medianAbsoluteDeviationPerDonor)[idx]
    idx=match(donors$donor, names(totalUMIsPerDonor))
    donors$totalUMIs=unlist(totalUMIsPerDonor)[idx]
	
    #order the data by median UMI from smallest to largest.
    donors=donors[order(donors$medianUMIs),]    
    donors$donor=factor(donors$donor, levels=donors$donor)
    maxY=max(donors$medianUMIs)*1.1
    
    p=ggplot(donors, aes(x = donor, y = medianUMIs)) +
    	geom_point(size = 4) +
    	scale_y_continuous(expand = c(0, 0), limits = c(0, maxY)) +
    	labs(x = "", y = "median UMIs per donor") +
    	ggtitle(expName) +
    	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
    		  axis.text.y = element_text(size = 10),
    		  axis.title = element_text(size = 12),
    		  plot.title = element_text(size = 14),
    		  legend.position = c(0.8, 0.8)) +
    	theme(panel.background = element_rect(colour = "black"))
	print (p)
    
	
	#stack up the UMIs per donor into a dataframe suitable for ggplot2
    df=stack(umisPerDonor)
    names(df) <- c("umisPerDonor", "donor")
    
    sortedDonorNames=names (sort (tapply(df$umisPerDonor, FUN=median, INDEX=df$donor)))
    df$donor=factor(df$donor, levels=sortedDonorNames)
    
    plotOne<-function (df, expName, maxY=NULL) {
    	if (is.null(maxY))
    		maxY=max(df$umisPerDonor)
    	p=ggplot(data=df, aes(x=donor, y = umisPerDonor)) +
    		geom_boxplot() +
    		labs(y = "UMIs") +
    		ggtitle(expName) +
    		coord_cartesian(ylim=c(0, maxY)) +
    		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10), 
    			  axis.text.y = element_text(size = 10),
    			  axis.title = element_text(size = 12),
    			  plot.title = element_text(size = 14))
    	print (p)
    }
    
    #plotOne(df, expName)
    #plotOne(df, expName, maxY=2*max(unlist(medianUMisPerDonor)))
    
    #donor UMI distribution
    d=vegan::diversity(donors$totalUMIs)
    eq=d/log(dim(donors)[1])
	
    donors=donors[order(donors$totalUMIs, decreasing = T),]    
    donors$donor=factor(donors$donor, levels=donors$donor)
    
    strTitle=paste("Distribution of UMIs across donors\n",dim (donors)[1], " donors; ", round (sum (donors$totalUMIs)/1e6,1), "M UMIs; SW Div: ", sprintf("%.2f",round(d,2)), "; SW Eq: ", sprintf("%.2f", round (eq,3)), sep="")

    p=ggplot(donors, aes(x = donor, y = totalUMIs / 1e6, fill = 'light blue')) +
    	geom_bar(stat = "identity", color = "black") +
    	labs(y = "Total UMIs [millions]", x = "", fill = "") +
    	ggtitle(strTitle) +
    	scale_fill_identity() +
    	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10), 
    		  axis.text.y = element_text(size = 10),
    		  axis.title = element_text(size = 12),
    		  plot.title = element_text(size = 14))
    
    print (p)
    
    p=ggplot(donors, aes(x = count, y = medianUMIs)) +
    	geom_point(size = 4) +
    	labs(x = "number of cels per donor", y = "median UMIs per donor") +
    	coord_cartesian(ylim=c(0, max (donors$medianUMIs)*1.1)) +
    	ggtitle(expName) +
    	theme(axis.text.x = element_text(size = 10), 
    		  axis.text.y = element_text(size = 10),
    		  axis.title = element_text(size = 12),
    		  plot.title = element_text(size = 14)) +
    	theme(panel.background = element_rect(colour = "black"))
    
    print (p)	
        
    r=list(stats=data.frame(diversity=d, equitability=eq), donors=donors)
    return (r)
}


# Why a pre-cellbender DGE file?
# Donor assignment does not benefit from UMI cleanup that is performed by cellbender,
# so we need the denominator to reflect the number of UMIs in the experiment that donor assignment 
# has access to.
plotRatioUMIsCaptuedToCellSize<-function (dgeRawSummaryFile, cellDF, cellsToKeep=NULL) {
	#TO make R CMD CHECK happy.  Blame ggplot2.
	density=.x=NUM_TRANSCRIPTS=fractionSNPUMIs=donor=NULL
	
	#quit early if no data.
	if ((!is.null(cellsToKeep) & length(cellsToKeep)==0) | is.null(dgeRawSummaryFile))
		return (NULL)
	
    dge=read.table(dgeRawSummaryFile, header=T, stringsAsFactors=F, sep="\t")
    df=data.table::copy(cellDF)
    if (!is.null(cellsToKeep)) {
        idx=sort(match(cellsToKeep, df$cell_barcode))
        if (length(idx)>0) df = df[idx,]
       
    }

    idx=match(dge$CELL_BARCODE, df$cell_barcode)
    
    dge$snpUMIs=df[idx,]$singlet_num_inform_umis
    dge$donor=df[idx,]$donor
    dge$likePerUMI=df[idx,]$best_like/df[idx,]$singlet_num_inform_umis
    dge$fractionSNPUMIs=dge$snpUMIs/dge$NUM_TRANSCRIPTS
    idxNA=which(is.na(dge$snpUMIs))
    if (length(idxNA)>0) dge=dge[-idxNA,]

    # after_stat(!!str2lang("density")) is a replacement for using ..density.. that doesn't make R CMD check mad.
    # https://stackoverflow.com/questions/74756655/how-to-replace-the-dot-dot-notation-in-ggplot2geom-histogramy-density
    
    plotHistWithlog10Axis <- function(df, var, xlab, ylab, main, xLims) {
    	p=ggplot(df, aes(x = !!sym(var))) + 
    		geom_histogram(aes(y = after_stat(!!str2lang("density"))), bins = 100, colour = 1, fill = "white") +
    		geom_density() +
    		scale_x_log10(labels = scales::trans_format('log10', scales::math_format(10^.x)), limits=xLims) +
    		labs(x = xlab, y = ylab, title = main) +
    		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10), 
    			  axis.text.y = element_text(size = 10),
    			  axis.title = element_text(size = 12),
    			  plot.title = element_text(size = 14))
    	
    	return (p)
    }
    
    maxXLog10=ceiling(log10(max(dge$NUM_TRANSCRIPTS, dge$singlet_num_inform_umis)))+1
    xLims=c(1, 10^maxXLog10)
    p1=plotHistWithlog10Axis(dge, var='NUM_TRANSCRIPTS', xlab="UMIs per cell [log10]",  ylab="number of cells", main="Expressed UMIs", xLims=xLims)
    p2=plotHistWithlog10Axis(df, var='singlet_num_inform_umis', xlab="Informative UMIs per cell [log10]",  ylab="number of cells", main="Informative UMIs", xLims=xLims)
    
    maxY=max(1, max (dge$fractionSNPUMIs)*1.1)
    p3=ggplot(dge, aes(x = NUM_TRANSCRIPTS, y = fractionSNPUMIs)) +
    	geom_point(size = 0.25) +
    	labs(x = "num transcripts", y = "fraction of UMIs that are informative") +
    	coord_cartesian(ylim=c(0, maxY)) +
    	ggtitle("Informative UMIs") +
    	theme(axis.text.x = element_text(size = 10), 
    		  axis.text.y = element_text(size = 10),
    		  axis.title = element_text(size = 12),
    		  plot.title = element_text(size = 14)) +
    	theme(panel.background = element_rect(colour = "black")) +
    	scale_x_log10(labels = scales::trans_format('log10', scales::math_format(10^.x)))
    
    
    sortedDonorNames=names (sort (tapply(dge$fractionSNPUMIs, FUN=median, INDEX=dge$donor)))
    dge$donor=factor(dge$donor, levels=sortedDonorNames)
    
    p4=ggplot(dge, aes(x = donor, y = fractionSNPUMIs)) +
    	geom_boxplot() +
    	ggtitle("Informative UMIs per donor") +
    	labs(x = "", y = "fraction of UMIs that are informative") +
    	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
    		  axis.text.y = element_text(size = 10),
    		  axis.title = element_text(size = 12),
    		  plot.title = element_text(size = 14))
    
    # the xlimits on the histograms are slightly cranky but work fine.
    suppressWarnings(gridExtra::grid.arrange(p1, p3, p2, p4, nrow = 2))
    
}

#for each donor, sum all cells expression into a single "meta cell".
generateMetaCells<-function (cellDonorMap, dgeFile, outMetaCellFile, selectedCellsForMetaCellFile=NULL) {

    if (!is.null(selectedCellsForMetaCellFile)) {
        z=read.table(selectedCellsForMetaCellFile, header=F, stringsAsFactors=F)$V1
        z=intersect(z, cellDonorMap$cell)
        cellDonorMap=cellDonorMap[match(z, cellDonorMap$cell),]
    }
    #donor=donors[1]
    makeMetaCell<-function (donor, cellDonorMap, a) {
        cells=cellDonorMap[cellDonorMap$bestSample==donor,]$cell
        idx=match(cells, colnames(a))
        idx=idx[!is.na(idx)]
        x=a[,idx, with=F]
        expression=rowSums(x)
        names(expression)=a$GENE
        return (expression)
    }
    #If there's no data, skip this method.
	if (dim (cellDonorMap)[1]==0)
		return (NULL)
    
    a=DropSeq.utilities::read_dge_txt(dgeFile)
    donors=sort(unique(cellDonorMap$bestSample))
    r=lapply(donors, makeMetaCell, cellDonorMap, a)
    r=do.call(cbind, r)
    colnames(r)=donors
    #put the genes as the first column to look like DGE data.
    result=data.frame(GENE=rownames(r), r, stringsAsFactors = F, check.names = F)
    write.table(result, outMetaCellFile, row.names=F, col.names=T, quote=F, sep="\t")
}

writeSummaryStats<-function (summaryStats, outSummaryStatsFile) {
    #formatting
	summaryStats$cell_equitability=round(summaryStats$cell_equitability, 2)
    summaryStats$diversity=round(summaryStats$diversity,2)
    summaryStats$equitability=round(summaryStats$equitability,2)
    summaryStats$reads_per_umi=round(summaryStats$reads_per_umi,2)
    if (!is.null(outSummaryStatsFile)) write.table(summaryStats, outSummaryStatsFile, row.names=F, col.names = T, quote=F, sep="\t")
}

#' Standard set of analysis for single cell donor assignment.
#' 
#' This generates a QC report with a number of useful plots, along with optional
#' outputs that are useful for downstream analysis.  
#'
#' @param expName The experiment name.  This is used to label plots.
#' @param likelihoodSummaryFile The dropulation single donor likelihood file.  
#' @param doubletLikelihoodFile The dropulation doublet likelihood file.  
#' @param dgeSummaryFile The DGE summary file.  If cleanup steps have been run on the 
#' data, the most processed data should be used at this step to preserve that in the metacell summary.  
#' A DGE summary file is a tab separated file with the follow columns:
#' * CELL_BARCODE: the cell barcode sequence CELL_BARCODE that identifies a single cell uniquely
#' * NUM_TRANSCRIPTS: the total number of UMIs captured by a cell.  This is the sum of expression across all genes. 
#' (optional)
#' @param dgeRawSummaryFile The DGE summary file.  This file should contain the number of UMIs original observed in the 
#' data before any cleanup steps like cellbender. See dgeSummaryFile for file format.  (optional)
#' @param readsPerCellFile The reads per cell file. This file is tab delimited, has no header, and has two columns:
#' 1. The number of reads for a cell barcode
#' 2. The cell barcode identity.  
#' (optional)
#' @param censusFile The census results file. (optional)
#' @param expectedSamplesFile A file containing the expected list of donors in the experiment. The file has no header,
#' and a single donor ID on each line.  
#' @param outFileLikelyDonors Output of which donors are in the experiment based on the Dropulation results.
#' @param outDonorToCellMap A map from each cell barcode to the best donor.  This is restricted to cells
#' that are not doublets and have passed FDR thresholds.  This represents subset of input cells that should
#' be used in downstream analysis that requires donor labels.
#' @param outPDF The output PDF file. (optional)
#' @param outSummaryStatsFile The output summary statistics file. (optional)
#' @param minimumFractionDonor To determine which donors are in the data, what fraction of cells must contain a donor for the 
#' donor to be included in the output.  This controls for a low level of miscalls of donor identity.
#' @param alpha The p-value threshold for doublets.
#' @param rescueDiffuseDoublets [EXPERIMENTAL] Should doublets where the 2nd donor of the pair isn't clear be rescued when possible?  These
#' cells usually have higher amounts of ambient RNA that makes many potential doublet pairs equally likely.  If set to true,
#' attempts to differentiate cells assigned to expected donors from noisy cells to recover a subset of barcodes with marginal scores
#' as singlets.   This method needs a pool where many donors that are evaluated are not expected in the pool (a set of null results).  
#' Without that null result, this option will erroneously rescue cells as singlets - use with caution!
#' @param minNumUMIs Filter cells with fewer than this number of UMIs. (optional)
#' @param dgeFile The input DGE file to generate meta cell output. If cleanup steps have been run on the 
#' data, the most processed data should be used at this step to preserve that in the metacell summary.  This is a tab delimited
#' flat file with each row containing a gene and each column containing a cell barcode.  The first column name is "GENE", followed by 
#' cell barcode identifiers.  (optional)
#' @param outMetaCellFile Sums the expression across cells of each donor to generate a new expression file for donors instead 
#' of cells. (optional)
#' @param selectedCellsForMetaCellFile Limits meta cell expression to a subset of donors.
#' @param outCellBarcodesFile cell barcodes for donors in the experiment (optional)
#' @param anonymizeDonors If set to true, donor IDs are changed to DONOR_1, DONOR_2, etc.  Purely for sharing 
#' visualizations of not-yet-public data.
#' @import grDevices graphics utils stats data.table RColorBrewer  vegan
#' @suggest DropSeq.utilities
#' @export
donorAssignmentQC<-function (expName="", likelihoodSummaryFile, doubletLikelihoodFile, dgeSummaryFile=NULL, dgeRawSummaryFile=NULL,
readsPerCellFile=NULL, censusFile=NULL, expectedSamplesFile, outFileLikelyDonors=NULL, outDonorToCellMap=NULL, outPDF=NULL,
outSummaryStatsFile=NULL, minimumFractionDonor=0.002, alpha=0.05, rescueDiffuseDoublets=F, minNumUMIs=0, dgeFile=NULL,
    outMetaCellFile=NULL, selectedCellsForMetaCellFile=NULL, outCellBarcodesFile=NULL, anonymizeDonors=FALSE) {

    validateFilesExist(likelihoodSummaryFile, doubletLikelihoodFile, dgeRawSummaryFile, dgeSummaryFile, censusFile, expectedSamplesFile)

    summaryStats=data.frame(expName=expName)
    #TODO: maybe not hard code this later?
    doubletPvalueThreshold=0.9

    if (!is.null(outPDF)) pdf(outPDF)

    if (is.null(doubletLikelihoodFile)) {
        # plotCellsByProbabilityNew(likelihoodSummaryFile, cellsToKeep=NULL, fdrThreshold=alpha)
    } else {
		
    	#all filtering/labeling of doublet/singlet status occurs here.    
    	#note: there's a parameter to convert all diffuse doublet to singlets, but it's not used.
    	z=getSingletDoubletDF(likelihoodSummaryFile, doubletLikelihoodFile, expectedSamplesFile=expectedSamplesFile, 
    							   doubletPvalue=doubletPvalueThreshold, bestPairPvalue=doubletPvalueThreshold, 
    							   ignoreDiffuseDoublets=FALSE, fdrThreshold=alpha, anonymizeDonors=anonymizeDonors)
    	
    	#pull out the dataframe and expected samples, which may have been remapped to anon IDs.
    	cellDF=z$cellDF
    	expected_samples=z$expected_samples
    	
    	summaryStats$total_cells=dim(cellDF)[1]
        if (!is.na(minNumUMIs)) {
        	cellDF=cellDF[cellDF$singlet_num_inform_umis>=minNumUMIs,]
        }

    	#this modifies the data frame to set some diffuse doublets to singlet
    	if (rescueDiffuseDoublets) {
    		plotAmbientRescue(cellDF, fdrThreshold=alpha)	
    		plotAmbientRescueDonorLikelihoods(cellDF)
    		z=plotAmbientRescueFDRPlusDonorLike(cellDF, fdrThreshold = alpha)
    		cellDF=plotOptimizedAmbientDonor (cellDF, fdrThreshold=alpha)	
    	}
    	
        doubletStats=evaluateDoubletRateSimple (cellDF, fdrThreshold=alpha)
        pctDS=round(doubletStats[c("fracDoubletAll", "fracDiffuseContaminationDoublets", "fracConfidentDoublets", "fracAllImpossibleDonors", "fracFdrFilteredImpossibleDonors", "fracAllDoubletsFilteredImpossibleDonors")]*100,2)
    	names(pctDS)=c("pct_all_doublets", "pct_diffuse_contam_doublets", "pct_confident_doublets", "pct_impossible_donors", "pct_fdr_impossible_donors", "pct_doublet_filtered_impossible_donors")
        summaryStats=cbind (summaryStats, pctDS)
        
        plotDoubletProbability(df=cellDF, expName = expName, summaryStats=summaryStats, doubletPvalueThreshold=doubletPvalueThreshold)
        
        #add the number of singlets
        singlets=cellDF[cellDF$label_simple=="singlet",]$cell_barcode
        summaryStats$singlets=length(singlets)
		
        plotSummaryStats(summaryStats)
        	
        plotCellsByProbabilityNew(df=cellDF, fdrThreshold=alpha)
		
        #plot some estimates of error rates.
        plotFractionImpossibleAllelesFromDoublets (df=cellDF, expName, fdrThreshold=alpha)
        
        #plot the average likelihood per UMI of singlets, doublets, etc.
        plotAverageLikelihood (df=cellDF)
        
        plotFractionConfidentDoubletsFromSingleLikelihoodFit(cellDF)
        
        #plot the average likelihood partitioned by donor assignment.
        # plotAverageLikelihoodPerDonor(df=cellDF, minimumFractionDonor)
    }

    r=plotCommonDonors(df=cellDF, minimumFractionDonor)
    if (!is.null(censusFile)) {
        plotCensusComparison(censusFile, r$summary, expName)
    }
	
    plotCommonDonorsFdrDistribution(df=cellDF, minimumFractionDonor)
    
    donors=r$summary
    cellDonorMap=r$cellDonorMap
    summaryStats$assignable_singlets=dim(cellDonorMap)[1]

    #get cell diversity and add to summary stats
    cellEquitability = getCellEquitability(cellDonorMap)
    summaryStats$cell_equitability=cellEquitability
    
    plotCommonDonorAssignmentsWithUnexpected(donors, expected_samples, expName)
    
    cellsToKeep=cellDonorMap$cell
    zz=calculateMeanUMIsPerDonor(expName, cellDF, dgeSummaryFile, donors, cellsToKeep=cellsToKeep)
    
    donors=zz$donors
    summaryStats=cbind(summaryStats, zz$stats)
    summaryStats$totalUMIs=sum (donors$totalUMIs)
    
    plotRatioUMIsCaptuedToCellSize(dgeRawSummaryFile, cellDF, cellsToKeep=cellDonorMap$cell)
    summaryStats$reads_per_umi=calculateReadsPerUMI(readsPerCellFile, dgeSummaryFile, cellDonorMap)
    if (!is.null(outSummaryStatsFile)) writeSummaryStats(summaryStats, outSummaryStatsFile)
    if (!is.null(outPDF)) dev.off()
    if (!is.null(outFileLikelyDonors)) write.table(donors, outFileLikelyDonors, row.names=F, col.names=T, quote=F, sep="\t")
    if (!is.null(cellDonorMap)) write.table(cellDonorMap, outDonorToCellMap, row.names=F, col.names=T, quote=F, sep="\t")
    if (!is.null(dgeFile) && !is.null(outMetaCellFile)) {
        generateMetaCells(cellDonorMap, dgeFile, outMetaCellFile, selectedCellsForMetaCellFile=selectedCellsForMetaCellFile)
    }
    if (!is.null(outCellBarcodesFile)) {
        write.table(cellDonorMap$cell, outCellBarcodesFile, row.names=F, col.names=F, quote=F)
    }
}



#' Compute and plot the conditional doublet rate.
#'
#' Does the doublet rate change with the number of UMIs?
#' It might be over-estimated for cells with few UMIs.
#' From internal synthetic mixture experiments, the doublet false positive rate tends 
#' to be higher for cells with fewer than ~ 100 or ~ 150 informative UMIs.

#' @param cellDF dataframe of cell statistics and labels generated by getSingletDoubletDF
#'
#' @return
plotConditionalDoubletRate<-function(cellDF) {
	#break the number of informative UMIs into deciles.
	cellDF$umiDecile=cut(cellDF$doublet_num_inform_umis, breaks=quantile(cellDF$doublet_num_inform_umis, probs=seq(0,1,0.1), na.rm=T), include.lowest=T)
	
	breaks <- quantile(cellDF$doublet_num_inform_umis, probs=seq(0, 1, 1/numBreaks), na.rm=TRUE)
	labels <- scales::label_number()(breaks)
	labels <- sapply(seq_along(labels)[-length(labels)], function(i) paste0("[", labels[i], ",", labels[i+1], "]"))
	
	# Create deciles
	cellDF$umiDecile <- cut(cellDF$doublet_num_inform_umis, 
							breaks=breaks, 
							labels=labels, 
							include.lowest=TRUE)
	
	#a function that gets the doublet rate for a given decile.
	getDoubletRate<-function (decile=NULL, df) {
		if (is.null(decile)) {
			idx=1:dim(df)[1]
			decile=NA
		} else {
			idx=which(df$umiDecile==decile)	
		}
		df2=df[idx,]
		doublets=sum(df2$label=="confident_doublet")
		diffuse=sum(df2$label=="diffuse_contamination")
		denominator=dim(df2)[1]
		df=data.frame(decile=decile, doubletRate=doublets/denominator, diffuseRate=diffuse/denominator, num_cells=dim(df2))
		return (df)
	}
	
	#apply the function to each decile.
	r=lapply(levels(cellDF$umiDecile), getDoubletRate, cellDF)
	r=do.call(rbind, r)
	r$doubletRateSmooth=lowess(r$doubletRate)$y
	r$diffuseRateSmooth=lowess(r$diffuseRate)$y
	r$decile=factor(r$decile, levels=levels(cellDF$umiDecile))
	
	#plot the results.
	ggplot(r, aes(x = decile)) +
		geom_point(aes(y = doubletRate, color = "doublet"), size = 3) +
		geom_point(aes(y = diffuseRate, color = "diffuse"), size = 3) +
		geom_line(aes(y = doubletRateSmooth, color = "doublet", group = 1, linetype = "Smoothed"), size = 1) +
		geom_line(aes(y = diffuseRateSmooth, color = "diffuse", group = 2, linetype = "Smoothed"), size = 1) +
		scale_color_manual(values = c("doublet" = "red", "diffuse" = "purple")) +
		scale_linetype_manual(values = c("Smoothed" = "dotted")) +
		labs(x = "number of informative UMIs", y = "Rate", color = "Type", linetype = "") +
		theme_minimal() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1))
	
}


#' Plot how well separated the singlets and doublets are.
#'
#' Fits the singlet likelihoods to a constrained B-spline model and then plots the residuals of the doublets.
#' The residuals are the difference between the observed likelihood and the expected likelihood based on the singlet model.
#' 
#' @param cellDF dataframe of cell statistics and labels generated by getSingletDoubletDF
#' @import cobs
plotFractionConfidentDoubletsFromSingleLikelihoodFit<-function (cellDF) {
	
	df=data.frame(cellDF)
	df$normalizedBestLikelihood=df$best_like/df$singlet_num_inform_umis
	
	dfSinglet=df[df$label=="single_pass_FDR",]
	
	# Fit a constrained B-spline model ensuring monotonicity
	fit <- cobs::cobs(log10(dfSinglet$singlet_num_inform_umis), dfSinglet$normalizedBestLikelihood,
				constraint = "increase",
				nknots = 50,  # Adjust the number of knots as needed
				degree = 2,   # Degree must be 1 or 2
				print.mesg=FALSE) # quite, please.  
	
	# Extract the fitted values get residuals for singlets
	prediction <- predict(fit, log10(dfSinglet$singlet_num_inform_umis))
	fitted_values <- prediction[, "fit"]
	dfSinglet$fitted_values <- prediction[, "fit"]
	dfSinglet$residuals <- dfSinglet$normalizedBestLikelihood - fitted_values
	
	# predict the doublets, calculate the residuals, and flag outliers
	dfDoublet=df[df$label=="confident_doublet",]
	prediction <- predict(fit, log10(dfDoublet$singlet_num_inform_umis))
	fitted_values <- prediction[, "fit"]
	
	# Calculate the residuals
	dfDoublet$residuals <- dfDoublet$normalizedBestLikelihood - fitted_values
	
	dfSinglet$residual_type <- "Singlet"
	dfDoublet$residual_type <- "Doublet"
	
	# Combine the residuals into a single data frame
	df_combined <- rbind(
		data.frame(residuals = dfSinglet$residuals, type = dfSinglet$residual_type, singlet_num_inform_umis=dfSinglet$singlet_num_inform_umis),
		data.frame(residuals = dfDoublet$residuals, type = dfDoublet$residual_type, singlet_num_inform_umis=dfDoublet$singlet_num_inform_umis)
	)
	
	
	theme_custom <- function() {
		theme_minimal() +
			theme(
				axis.title.x = element_text(size = 8),  # Smaller X-axis title
				axis.title.y = element_text(size = 8),  # Smaller Y-axis title
				plot.title = element_text(size = 10),    # Smaller main title
				legend.position = "top"                  # Legend at the top
			)
	}
	
	fitPlot <- ggplot() +
		# Plot Singlet points first
		geom_point(data = dfSinglet, aes(x = singlet_num_inform_umis, y = normalizedBestLikelihood), 
				   color = "green", size = 0.25) +
		# Plot Doublet points second
		geom_point(data = dfDoublet, aes(x = singlet_num_inform_umis, y = normalizedBestLikelihood), 
				   color = "red", size = 0.1) +
		# Add the fitted line
		geom_line(data = dfSinglet, aes(x = singlet_num_inform_umis, y = fitted_values), color = "black") +
		# Apply log transformation to the x-axis
		scale_x_continuous(trans = 'log10') +
		# Add labels and title
		labs(title = "Singlet likelihood fit",
			 x = "Number of informative UMIs",
			 y = "Normalized likelihood (likelihood/num inform UMIs)") +
		theme_custom()
	
	residualLimits=range(df_combined$residuals)*1.2
	
	residualDensityUnfiltered <- ggplot(df_combined, aes(x = residuals, fill = type)) +
		geom_density(alpha = 0.5) +
		scale_fill_manual(values = c("Singlet" = "green", "Doublet" = "red")) +
		labs(title = "Residuals of Singlets and Doublets",
			 x = "Residuals",
			 y = "Density") +
		theme_custom() +
		xlim(residualLimits) + 
		theme(legend.position = "top")
	
	##############################################################
	#the residuals fit poorly at the low end of the model.
	##############################################################
	
	#capture where all points have less than 3SD from the mean.
	outliers <- which(abs(dfSinglet$residuals) > 3 * sd(dfSinglet$residuals))
	log10Threshold=log10(quantile (dfSinglet[outliers,]$singlet_num_inform_umis, seq(0,1,0.1))["90%"])
	
	#find a threshold that removes 90% of the outliers.
	# plot (log10(dfSinglet$singlet_num_inform_umis), dfSinglet$residuals)
	# points(log10(dfSinglet$singlet_num_inform_umis)[outliers], dfSinglet$residuals[outliers], col = "red")
	# title(main=paste("Num outliers", length(outliers)))
	# abline (v=log10Threshold, col="red")
	
	dfSinglet$outliers <- ifelse(1:nrow(dfSinglet) %in% outliers, TRUE, FALSE)
	
	outlierPlot <- ggplot(dfSinglet, aes(x = singlet_num_inform_umis, y = residuals)) +
		geom_point(size=0.5) +
		geom_point(data = subset(dfSinglet, outliers), aes(x = singlet_num_inform_umis, y = residuals), color = "red", size=0.5) +
		geom_vline(xintercept = 10^log10Threshold, color = "red") +
		labs(title = paste("Num outliers:", length(outliers)),
			 x = "Number of informative UMIs",
			 y = "Residuals") +
		scale_x_continuous(trans = 'log10') +
		theme_custom()
	
	#print(outlierPlot)
	
	#replot the residuals with the threshold.
	df_combined_filtered=df_combined[log10(df_combined$singlet_num_inform_umis)>log10Threshold,]
	
	residualDensityFiltered <- ggplot(df_combined_filtered, aes(x = residuals, fill = type)) +
		geom_density(alpha = 0.5) +
		scale_fill_manual(values = c("Singlet" = "green", "Doublet" = "red")) +
		labs(title = "Filtered Residuals of Singlets and Doublets",
			 x = "Residuals",
			 y = "Density") +
		theme_custom() + 
		xlim(residualLimits) +
		theme(legend.position = "top")
	
	
	gridExtra::grid.arrange(fitPlot, residualDensityUnfiltered, outlierPlot, residualDensityFiltered, nrow=2)
	
}


#maps singlet/doublet labels to colors for consistently colored plots
getLabelColors<-function () {
	cols=c('single pass FDR'='green', 'single fail FDR'='blue', 'confident doublet'='red', 'diffuse contamination'='purple')
	return (cols)
}

#' Plot the average penality per UMI for each cell
#' 
#' Color cells by their class.
#'
#' @param df The data frame generated by getSingletDoubletDF
#' @import ggplot2
#' @noRd
plotAverageLikelihood<-function (df) {
	#TO make R CMD CHECK happy.  Blame ggplot2.
	singlet_num_inform_umis=normalizedBestLikelihood=label=value=NULL
	
	df2=data.table::copy(df)
	df2$label=gsub("_", " ", df2$label)
	df2$normalizedBestLikelihood=df2$best_like/df2$singlet_num_inform_umis
	factorOrder=c("single pass FDR", "single fail FDR", "diffuse contamination","confident doublet")
	factorOrder=intersect(factorOrder, df2$label)
	df2$label=factor(df2$label, levels=factorOrder)
	
	p1=ggplot(df2, aes(x = singlet_num_inform_umis, y = normalizedBestLikelihood, color=label)) +
		geom_point(alpha=0.5, size=0.5) +
		scale_x_continuous(trans='log10') +
		xlab("Number of informative UMIs") +
		ylab("Normalized likelihood (likelihood/num inform UMIs)") +
		ggtitle("Donor assignment normalized likelihood") +
		theme_grey(base_size = 10) +
		guides(colour = guide_legend(override.aes = list(size=3))) +
		scale_colour_manual(values=getLabelColors(), name="")
	
	#calculate the average penalty by class.
	avg=tapply(X=df2$normalizedBestLikelihood, INDEX = df2$label, FUN = median)
	
	avg=data.frame(label=names(avg), value=round(as.numeric(avg),4))
	avg$label=factor(avg$label, levels=factorOrder)
	
	p2=ggplot(data=avg, aes(x=label, y=value, fill=label)) + 
		geom_bar(position = 'dodge', stat='identity', alpha=0.5) +
		geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25, size=4) +
		theme_grey(base_size = 10) +
		theme(legend.position="none") +
		ylab("Median Score") +
		xlab("") +
		scale_fill_manual(values=getLabelColors(), name="")
	
	gridExtra::grid.arrange(p1, p2, nrow = 2, heights=c(0.7, 0.3))
	
}

#' Plot the single cell FDR distribution of each donor
#'
#' @param df The data frame generated by getSingletDoubletDF
#' @inheritParams donorAssignmentQC
#' @noRd
plotAverageLikelihoodPerDonor<-function (df, minimumFractionDonor=0.005) {
	#TO make R CMD CHECK happy.  Blame ggplot2.
	donor=normalizedBestLikelihood=NULL
	
	
	rUnfiltered=sort(table(df$donor), decreasing=T)
	minNum=sum(rUnfiltered)* minimumFractionDonor
	rFiltered=rUnfiltered[rUnfiltered>minNum]
	commonDonors=names (rFiltered)
	
	df2=data.table::copy(df)
	df2$normalizedBestLikelihood=df2$best_like/df2$singlet_num_inform_umis
	subset_data=df2[df2$donor %in% commonDonors,]
	
	# Create the boxplot with -log10 transformation
	p=ggplot(subset_data, aes(x = donor, y = normalizedBestLikelihood)) +
		geom_violin() +
		labs(x = "", y = "Normalized likelihood (likelihood/num inform UMIs)") +
		ggtitle("Normalized likelihood distribution for each common donor") +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10), 
			  axis.text.y = element_text(size = 10),
			  axis.title = element_text(size = 12),
			  plot.title = element_text(size = 14)) +
		theme(panel.background = element_rect(colour = "black"))
	print (p)
}


plotSummaryStats<-function (summaryStats) {
	# quick plot of impossible donors
	# Now combine with the doublet rates.
	par(mfrow=c(2,1))
	stats=c(summaryStats$pct_all_doublets, summaryStats$pct_diffuse_contam_doublets, summaryStats$pct_confident_doublets)
	z=barplot (stats, main="Doublet Rates", names.arg=c("all", "diffuse contamination", "confident"), col="light blue")
	text (y=max(stats)/2, x=z[,1], labels=paste(round(stats,1), "%", sep=""), cex=1.5)        
	stats=c(summaryStats$pct_impossible_donors, summaryStats$pct_fdr_impossible_donors, summaryStats$pct_doublet_filtered_impossible_donors)
	z=barplot (stats, names.arg=c("all", "+fdr filtered", "+doublet filtered"), main="% impossible donors")
	text (y=max(stats)/2, x=z[,1], labels=paste(round(stats,1), "%", sep=""), cex=1.5)        
	par(mfrow=c(1,1))
	
}

getCellEquitability<-function (cellDonorMap) {
	donorCellCounts=table (cellDonorMap$bestSample)
	d=vegan::diversity(donorCellCounts)
	eq=d/log(length(donorCellCounts))
	return (eq)
} 

labelSingletsInDonorList<-function (b, donorListFile) {
	d=read.table(donorListFile, header = F, stringsAsFactors = F)$V1
	lastIdx=which(colnames(b)=="population_average_likelihood")
	cols=c(colnames (b)[1:lastIdx], "in_donorlist", colnames(b) [(lastIdx+1):dim (b)[2]])
	b$in_donorlist=F
	b=b[,match(cols, colnames(b))]
	idx=which(!is.na(match(b$bestSample, d)))
	b[idx,]$in_donorlist=T
	return (b)
}

############################################################################
# Fraction of observed alleles in doublet/singlet call that are impossible.
############################################################################

#' Plots that describe the errors in the doublet donor assignment
#' 
#' Given the doublet analysis, partition cell barcodes into confident singlets,
#' confident doublets, and ambient doublets.
#' For each cell, calculate the fraction of "impossible alleles" - that is, 
#' alleles that are observed that could not have been generated by the donor.
#' For single cells, this is an approximation of the error rate, though perhaps 
#' not the most sensitive one.
#' @param df A data frame of cell level metrics generated by getSingletDoubletDF
#' @param fdrThreshold The desired false discovery rate
#' @param expName An experiment name to use as a title for the plot
#' @import ggplot2
#' @noRd
plotFractionImpossibleAllelesFromDoublets<-function (df, expName, fdrThreshold=0.05) {
	#TO make R CMD CHECK happy.  Blame ggplot2.
	frac_s1_wrong=label=frac_s2_wrong=NULL
	
    x=data.table::copy(df)
    
    #small adjustment for doublet files with faked doublet calls
    x[is.na(x)] = 0
    
    #restrict to cells with a moderate numbers or large of observations (exclude first quantile)
    x2=x[which(x$doublet_num_inform_umis>=quantile (x$doublet_num_inform_umis) [2]) ,]
	
    #restrict singlets to those that pass FDR
    idx=which( (x$label=="single_pass_FDR") | x$label_simple!="singlet")
    x=x[idx,]
    
    df2=data.table::copy(x)
    df2$label=gsub("_", " ", df2$label)
    factorOrder=c("single pass FDR", "single fail FDR", "diffuse contamination","confident doublet")
    factorOrder=intersect(factorOrder, df2$label)
    df2$label=factor(df2$label, levels=factorOrder)
    
    singletErrorRateString=paste(round (median (x[x$label=="single_pass_FDR",]$frac_s1_wrong)*100,2), "%", sep="")
    
    strTitle=paste("Singlets (FDR<=", fdrThreshold, ") + doublets", "\nSinglet Error Rate ", singletErrorRateString, sep="")
    p1 <- ggplot2::ggplot(data=df2, aes(x=frac_s1_wrong, group=label, fill=label)) +
        ggplot2::geom_density(adjust=1.5, alpha=.2) +
        ggplot2::ggtitle (paste(expName, strTitle, sep="\n")) +
        xlab("fraction of alleles that could not be generated by donor one") +
    	scale_fill_manual(values = getLabelColors())

    p2 <- ggplot2::ggplot(data=x, aes(x=frac_s1_wrong, y=frac_s2_wrong, group=label, col=label)) +
        ggplot2::geom_point(alpha=0.2, size=1) +
        ggplot2::ggtitle (paste("Singlet Error Rate", singletErrorRateString)) +
        xlab("fraction donor one allele errors") +
        ylab("fraction donor two allele errors") +
    	scale_fill_manual(values = getLabelColors())

    #gridExtra::grid.arrange(p1, p2, nrow=2)
    gridExtra::grid.arrange(p1, nrow=1)

}


addFracImpossibleAlleles<-function (x) {
    x$sampleTwoMixtureRatio=1-x$sampleOneMixtureRatio
    x$frac_s1_wrong=x$sampleOneWrongAlleleCount/x$num_umi
    x$frac_s2_wrong=x$sampleTwoWrongAlleleCount/x$num_umi
    return (x)

}

#' Evaluate Doublet Rate (Simple)
#' 
#' Emits statistics about the doublet rate, as well as the fraction of cells coming from unexpected donors 
#' are assigned after various filters.
#' @param df A data frame of cell level metrics generated by getSingletDoubletDF
#' @param fdrThreshold The desired false discovery rate
#' @return A dataframe containing summary info about the fraction of doublets and enrichment.
#' @noRd
evaluateDoubletRateSimple<-function (df, fdrThreshold=0.05) {
	fracDoubletAll=length(which(df$label_simple=="doublet"))/dim (df)[1]
	fracConfidentDoublets=length(which(df$label=="confident_doublet"))/dim(df)[1]
	fracDiffuseContaminationDoublets=length(which(df$label=="diffuse_contamination"))/dim(df)[1]
		
	summaryDF=data.frame(fracDoubletAll=fracDoubletAll, fracConfidentDoublets=fracConfidentDoublets,fracDiffuseContaminationDoublets=fracDiffuseContaminationDoublets)
	impossibleStats=calculateFractionImpossibleSingleDonors(df, fdrThreshold)
	summaryDF=cbind(summaryDF, impossibleStats)
	return (summaryDF)
	
}

########################################################
# Consistent method to label the doublet/single status of each cell.
# 
########################################################
#reads in the singlets and doublets, return the median likelihood and best likelihood scores.  Flag cells as doublets based on the pvalue
#bin data by either some set number of qualtiles or a number of fixed bin sizes.
getSingletDoubletDF<-function (likelihoodSummaryFile, doubletLikelihoodFile, expectedSamplesFile=NULL, 
							   doubletPvalue=0.9, bestPairPvalue=0.9, ignoreDiffuseDoublets=FALSE, 
							   fdrThreshold=0.05, anonymizeDonors=FALSE) {
	a=read.table(likelihoodSummaryFile, header=T, stringsAsFactors = F, sep="\t")
	b=read.table(doubletLikelihoodFile, header=T, stringsAsFactors = F, sep="\t")
	
	#add the fraction impossible alleles up-front.
	b=addFracImpossibleAlleles(b)
	
	#note: the same cell barcodes are in both files
	symdiff=union (setdiff(a$cell, b$cell), setdiff(b$cell, a$cell))
	if (length(symdiff)>0) {
		problemCells=paste(head(symdiff), sep="", collapse=",")
		stop ("Something went wrong, donor assignment and doublets have different numbers of cells.  Check cell barcodes [", problemCells, "]")
	}
	both=intersect(b$cell, a$cell)
	a=a[match(both, a$cell),]
	b=b[match(both, b$cell),]
	
	#starting information from the singlet and doublet files	
	df=data.frame(cell_barcode=a$cell, donor=a$bestSample, singlet_num_inform_umis=a$num_umis, 
				  median_like=a$median_likelihood, best_like=a$bestLikelihood, FDR_pvalue=a$FDR_pvalue, 
				  sampleOneMixtureRatio=b$sampleOneMixtureRatio, sampleTwoMixtureRatio=b$sampleOneMixtureRatio,
				  frac_s1_wrong=b$frac_s1_wrong, frac_s2_wrong=b$frac_s2_wrong,
				  stringsAsFactors=F)
	
	#if the best pair pvalue column exists, add it!
	if ("best_pair_pvalue" %in% colnames (b)) {
		idx=match(df$cell_barcode, b$cell)
		df$doublet_pval=b[idx,]$doublet_pval
		df$best_pair_pvalue=b[idx,]$best_pair_pvalue
		df$doublet_num_inform_umis=b[idx,]$num_inform_umis
	}
	
	#add in expected donors if the file exists
	if (!is.null(expectedSamplesFile)) {
		e=read.table(expectedSamplesFile, header=F, stringsAsFactors=F, sep="\t")
		df$expected_donor=df$donor %in% e$V1
		#making this a factor makes life easier later!
		df$expected_donor=factor(df$expected_donor)
	}
	
	#label each cell with it's singlet/doublet status.
	df=labelDoublets(df, doubletPvalue, bestPairPvalue, ignoreDiffuseDoublets)
	
	#label each singlet with if it passes the current FDR.
	#important later when the FDR threshold for the subset of diffuse doublets that are rescued
	df=labelSinglets(df, fdrThreshold = fdrThreshold)
	
	
	#end early if no need to change donor names.
	if (!anonymizeDonors)
		return (list(cellDF=df, expected_samples=e))
		
	
	#Anonymize donors if requested.
	allDonorNames=unique (c(e$V1, df$donor))
	anonDonorNames <- paste("DONOR", as.numeric(factor(allDonorNames)), sep="_")
	df$donor=anonDonorNames[match(df$donor, allDonorNames)]
	e$V1=anonDonorNames[match(e$V1, allDonorNames)]
	return (list(cellDF=df, expected_samples=e))
}

labelSinglets<-function (x, fdrThreshold=0.05) {
	#NOTE: this is the correct way to do the modification without screwing up when the length is 0.
	#this is because we're accessing a vector of values and modifying them instead of rows of a data table.
	
	x$label[x$label_simple=="singlet" & x$FDR_pvalue>fdrThreshold]="single_fail_FDR"
	x$label[x$label_simple=="singlet" & x$FDR_pvalue<=fdrThreshold]="single_pass_FDR"
	return (x)
}

labelDoublets<-function (x, doubletPvalue=0.9, bestPairPvalue=0.9, ignoreDiffuseDoublets=FALSE) {
	x$label="singlet"
	x$label_simple="singlet"
	
	indexConfDoublet <- which(x$doublet_pval>=doubletPvalue & x$best_pair_pvalue>=bestPairPvalue)
	if (length(indexConfDoublet) > 0) {
		x[indexConfDoublet,]$label="confident_doublet"
	}
	
	if (!ignoreDiffuseDoublets) {
		indexDiffuseDoublet <- which(x$doublet_pval>=doubletPvalue & x$best_pair_pvalue<bestPairPvalue)
		if (length(indexDiffuseDoublet) > 0) {
			x[indexDiffuseDoublet,]$label="diffuse_contamination"
		}	
	}
	
	#################
	# SIMPLE LABEL
	# If diffuse doublets are ignored, only count confident doublets as doublets
	# Otherwise, count any doublet as a doublet
	#################
	
	if (ignoreDiffuseDoublets) {
		#only count confident doublet 
		idxSimpleDoublet=which(x$label=="confident_doublet")
	} else {
		#count confident and diffuse doublets
		idxSimpleDoublet=which(x$label=="diffuse_contamination" | x$label=="confident_doublet")
	}
	if (length(idxSimpleDoublet)>0)
		x[idxSimpleDoublet,]$label_simple="doublet"
	
	
	
	return (x)
}

#FDR threshold only applies to doublets.
calculateFractionImpossibleSingleDonors<-function (df, fdrThreshold=0.05) {
	getFracImpossible<-function (df) {
		expectedDonorCount=length(which(df$expected_donor==T))
		impossibleDonorCount=length(which(df$expected_donor==F))
		fracImpossibleDonors=impossibleDonorCount/(impossibleDonorCount+expectedDonorCount)
		return (fracImpossibleDonors)	
	}
	
	#unexpected donors in entire data set.
	fracAllImpossible=length(which(df$expected_donor==F))/dim(df)[1]
	
	#unexpected donors after filtering for FDR
	#a little more complicated now with modified FDR thresholds after rescue
	#The test set includes all doublets that pass FDR + singlets that pass FDR.
	doubletsPlusFDRSinglets=df[which( (df$label_simple=="doublet" & df$FDR_pvalue<fdrThreshold) | df$label=="single_pass_FDR"),]
	fracFdrFilteredImpossible=getFracImpossible(doubletsPlusFDRSinglets)
	
	#filter all doublets out and filter to only singlets that pass FDR
	fdrSinglets=df[which(df$label=="single_pass_FDR"),]
	fracAllDoubletsFilteredImpossible=getFracImpossible(fdrSinglets)
	
	#what's the largest unexpected donor's contribution?
	largestUnexpected=sort (table(df[df$expected_donor==F,]$singletDonor), decreasing=T)[1]
	
	sumDF=data.frame(fracAllImpossibleDonors=fracAllImpossible, fracFdrFilteredImpossibleDonors=fracFdrFilteredImpossible, fracAllDoubletsFilteredImpossibleDonors=fracAllDoubletsFilteredImpossible, stringsAsFactors = F)
	return (sumDF)	
	
	#x=df[df$label=="diffuse_contamination",]
	#plot (x$num_inform_umis, -log10(x$FDR_pvalue), main="Diffuse Doublets", col=factor(x$expected_donor))
}


########################################################
# Diffuse contamination doublet rescue plots/analysis.
########################################################

#' Plot to evaluate ambient rescue
#'
#' @param df The pre-processed data
#' @param fdrThreshold Threshold for filtering false positive results
#' @import scales
#' @noRd
plotAmbientRescue<-function(df, fdrThreshold=0.05) {
	#TO make R CMD CHECK happy.  Blame ggplot2.
	singlet_num_inform_umis=fdr_p=expected_donor=doublet_num_inform_umis=Var1=Freq<-NULL
	
	x=df[df$label=="diffuse_contamination",]
	
	x$fdr_p=-log10(x$FDR_pvalue)
	
	#how many points rescued?
	xx=x[x$FDR_pvalue<fdrThreshold,]
	#summary=as.data.frame.matrix(t(table (xx$expected_donor)))
	summary=data.frame(table (xx$expected_donor))
	
	
	p1=ggplot(x, aes(x=singlet_num_inform_umis, y=fdr_p, color=expected_donor)) +
		geom_point() +
		xlab("Number of singlet informative UMIs") +
		ylab("Donor assignment FDR [-log10]") +
		ggtitle ("Diffuse contamination rescue")	
	
	p2=ggplot(x, aes(x=doublet_num_inform_umis, y=fdr_p, color=expected_donor)) +
		geom_point() +
		xlab("Number of doublet informative UMIs") +
		ylab("Donor assignment FDR [-log10]")
	
	p3=ggplot(data=summary, aes(x=Var1, y=Freq)) +
		geom_bar(stat="identity", fill=scales::hue_pal()(2)) +
		geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), size=5) +
		xlab("Expected Donor") +
		ggtitle ("Cells with FDR < 0.05")	
	
	
	gridExtra::grid.arrange(p1, p2, p3, nrow = 3, heights=c(0.4, 0.4, 0.2))
	
}

#' Plot to evaluate ambient rescue using donor likelihoods.
#'
#' @param df The pre-processed data
#' @import scales
#' @noRd
plotAmbientRescueDonorLikelihoods<-function (df) {
	#TO make R CMD CHECK happy.  Blame ggplot2.
	median_likelihood=donor_likelihood=expected_donor=ratio=Var1=Freq=NULL
	
	x=df[df$label=="diffuse_contamination",]
	
	#don't -log10 the values, how do the plots look?
	xx=data.frame(cell_barcode=x$cell_barcode, FDR_pvalue=x$FDR_pvalue, median_likelihood=x$median_like, 
				  donor_likelihood=x$best_like,expected_donor=x$expected_donor,stringsAsFactors = F)
	
	p1=ggplot(xx, aes(x=median_likelihood, y=donor_likelihood, color=expected_donor)) +
		geom_point() +
		xlab("Population likelihood [log10]") +
		ylab("Donor best likelihood [log10]") +
		ggtitle ("Diffuse contamination rescue")
	
	#could you capture the ratio of donor likelihood to median, and compute donor cells that were some Z scores away?
	xx$ratio=xx$donor_likelihood-xx$median_likelihood
	
	#get the mean and sd of the unexpected donors
	m=mean (xx[xx$expected_donor==F,]$ratio)
	sd=sd (xx[xx$expected_donor==F,]$ratio)
	numSD=3
	threshold=m+(numSD*sd)
	
	p2= ggplot(xx, aes(x = ratio, colour = expected_donor)) +
		geom_density() + 
		xlab("Ratio of donor likelihood to median likelihood") +
		geom_vline(xintercept = threshold, linetype='dotted', color = "black", size=1) +
		ggtitle (paste("Diffuse contamination rescue by median likelihood \nMean + ",numSD, "SD from unexpected donor cells"))		
	
	
	summary=data.frame(table (xx[xx$ratio>threshold,]$expected_donor))
	
	p3=ggplot(data=summary, aes(x=Var1, y=Freq)) +
		geom_bar(stat="identity", fill=scales::hue_pal()(2)) +
		geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), size=5) +
		xlab("Expected Donor") +
		ggtitle (paste("Cells with donor likelihood- median likelihood >= ", round (threshold,3)))	
		
	gridExtra::grid.arrange(p1, p2, p3, nrow = 3, heights=c(0.4, 0.4, 0.2))
	
}

#' Plot to evaluate ambient rescue using both FDR filtering and donor likelihoods.
#'
#' @param df The pre-processed data
#' @param fdrThreshold Threshold for filtering false positive results
#' @import scales
#' @noRd
plotAmbientRescueFDRPlusDonorLike<-function (df, fdrThreshold=0.05) {
	#TO make R CMD CHECK happy.  Blame ggplot2.
	ratio=fdr_p=expected_donor=Var1=Freq=NULL
	
	x=df[df$label=="diffuse_contamination",]
	numCells=dim (x)[1]
	
	xx=data.frame(cell_barcode=x$cell_barcode, FDR_pvalue=x$FDR_pvalue, median_likelihood=x$median_like, 
				  donor_likelihood=x$best_like,expected_donor=x$expected_donor,stringsAsFactors = F)
	
	xx$fdr_p=-log10(xx$FDR_pvalue)
	
	#could you capture the ratio of donor likelihood to median, and compute donor cells that were some Z scores away?
	xx$ratio=xx$donor_likelihood-xx$median_likelihood
	
	#get the mean and sd of the unexpected donors
	m=mean (xx[xx$expected_donor==F,]$ratio)
	sd=sd (xx[xx$expected_donor==F,]$ratio)
	numSD=2
	threshold=m+(numSD*sd)
	
	strTitle=paste("Mean + ",numSD, "SD from unexpected donor cells [", round (threshold,3), "]")
	
	p1=ggplot(xx, aes(x=ratio, y=fdr_p, color=expected_donor)) +
		geom_point() +
		geom_vline(xintercept = threshold, linetype='dotted', color = "black", size=1) +
		geom_hline(yintercept = -log10(fdrThreshold), linetype='dotted', color = "black", size=1) +
		xlab("Ratio of donor likelihood to median likelihood") +
		ylab("Donor assignment FDR [-log10]") +
		ggtitle (paste("Diffuse contamination rescue\n",strTitle, sep=""))
	
	idx=which(xx$ratio>threshold & xx$FDR_pvalue<fdrThreshold)
	summary=data.frame(table (xx[idx,]$expected_donor))
	
	strTitleSuffix=paste("[", length(idx), "] of [", numCells, "] cells rescued [", round (length(idx)/numCells*100, 2), "%]", sep="")
	p2=ggplot(data=summary, aes(x=Var1, y=Freq)) +
		geom_bar(stat="identity", fill=scales::hue_pal()(2)) +
		geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), size=5) +
		xlab("Expected Donor") +
		ggtitle (paste("Combined diffuse doublet rescue\n",strTitleSuffix, sep=""))	
	
	gridExtra::grid.arrange(p1, p2, nrow = 2, heights=c(0.6, 0.4))

	#these diffuse doublets become singlets that pass FDR
	if (length(idx)>0) {
		cellsRescued=xx[idx,]$cell_barcode
		df[match(cellsRescued, df$cell_barcode),]$label="single_pass_FDR"
		df[match(cellsRescued, df$cell_barcode),]$label_simple="singlet"
	}
	
	return (df)
	
}

#' Plot to evaluate ambient rescue using FDR and donor likelihoods,
#' where the final thresholds are optimized to minimize the number
#' of false positives to be under the given FDR threshold.
#'
#' @param df The pre-processed data
#' @param fdrThreshold Threshold for filtering false positive results
#' @import scales
#' @noRd
plotOptimizedAmbientDonor<-function (df, fdrThreshold=0.05) {
	#TO make R CMD CHECK happy.  Blame ggplot2.
	ratio=fdr_p=expected_donor=Var1=Freq=NULL
	
	x=df[df$label=="diffuse_contamination",]
	numCells=dim (x)[1]
	
	#don't -log10 the values, how do the plots look?  Resonable!
	xx=data.frame(cell_barcode=x$cell_barcode, FDR_pvalue=x$FDR_pvalue, median_likelihood=x$median_like, 
				  donor_likelihood=x$best_like,expected_donor=x$expected_donor,stringsAsFactors = F)
	
	xx$fdr_p=-log10(xx$FDR_pvalue)
	
	#could you capture the ratio of donor likelihood to median, and compute donor cells that were some Z scores away?
	xx$ratio=xx$donor_likelihood-xx$median_likelihood
	
	#get the mean and sd of the unexpected donors
	m=mean (xx[xx$expected_donor==F,]$ratio)
	sd=sd (xx[xx$expected_donor==F,]$ratio)
	#this threshold is a starting spot for the optimization
	numSD=2
	threshold=m+(numSD*sd)
	
	###########################################
	# FOR FUN-ISH, optimize the two parameters
	###########################################
	
	scoreThresholds<-function (param, fdrThreshold, xx) {
		#each included expected donor scores a point, each included unexpected donor subtracts 1/fdrThreshold points
		t1=param[1]
		t2=param[2]
		idx=which(xx$ratio>=t1 & xx$fdr_p>=t2)
		#short circuit
		if (length(idx)==0) return (0)
		pos=length(which(xx[idx,]$expected_donor==T))
		neg=length(which(xx[idx,]$expected_donor==F))
		score=pos-(neg*1/fdrThreshold)
		return(score)
	}
	
	opt=optim(par=c(threshold,-log10(fdrThreshold)), fn=scoreThresholds, fdrThreshold=fdrThreshold, xx=xx, control=list(fnscale=-0.01), 
			  method="Nelder-Mead") 
	
	ratioThreshold=opt$par[1]
	fdr_pThreshold=opt$par[2]
	#updated fdr threshold from optimizaton
	fdrThreshold=10^-fdr_pThreshold
	
	idx=which(xx$ratio>=ratioThreshold & xx$fdr_p>=fdr_pThreshold)
	summary=data.frame(table (xx[idx,]$expected_donor))
	
	
	strTitle=paste("ratio [",round(ratioThreshold,3),"] FDR [", round (fdrThreshold,3), "]")
	p1=ggplot(xx, aes(x=ratio, y=fdr_p, color=expected_donor)) +
		geom_point() +
		geom_vline(xintercept = ratioThreshold, linetype='dotted', color = "black", size=1) +
		geom_hline(yintercept = fdr_pThreshold, linetype='dotted', color = "black", size=1) +
		xlab("Ratio of donor likelihood to median likelihood") +
		ylab("Donor assignment FDR [-log10]") +
		ggtitle (paste("Diffuse contamination rescue\nOptimized thresholds", strTitle))
	
	strTitleSuffix=paste("[", length(idx), "] of [", numCells, "] cells rescued [", round (length(idx)/numCells*100, 2), "%]")
	
	p2= ggplot(data=summary, aes(x=Var1, y=Freq)) +
		geom_bar(stat="identity", fill=scales::hue_pal()(2)) +
		geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), size=5) +
		xlab("Expected Donor") +
		ggtitle (paste("Optimized Threshold Combined diffuse doublet rescue\n",strTitleSuffix))	
	
	gridExtra::grid.arrange(p1, p2, nrow = 2, heights=c(0.6, 0.4))
	
	#these diffuse doublets become singlets that pass FDR
	if (length(idx)>0) {
		cellsRescued=xx[idx,]$cell_barcode
		df[match(cellsRescued, df$cell_barcode),]$label="single_pass_FDR"
		df[match(cellsRescued, df$cell_barcode),]$label_simple="singlet"
	}
	
	return (df)
}

