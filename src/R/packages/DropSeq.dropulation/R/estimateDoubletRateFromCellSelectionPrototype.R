# library (logger)
# library (ggplot2)

#' Estimate the expected number of doublets from the count of empty/nonempty cell barcodes
#' 
#' Given the count of empty and non-empty cell barcodes, fit data to a poisson distribution
#' where the empty count matches the number given, and the sum of counts >= 1 equals the non-empty count.
#' 
#' The number of nonEmpty cell barcodes is estimated by any cell selection software you'd like to use - 10x
#' and other scRNASeq methods will generate this number
#' 
#' The number of empty cell barcodes is slightly harder to estimate - it can be inferred using
#' cell level metrics of the number of UMIs and %intronic to find cell barcodes that have fewer UMIs and 
#' (in the case of nuclei) lower pct inronic.  See \code{\link{estimateEmptyCount}}
#' 
#' This estimates a lambda parameter for a Poisson distribution.  Assuming this distribution is reasonable,
#' it is them possible to estimate the doublet rate as the number of single observations divided by the number of
#' observations of one or more counts.
#' 
#' @param emptyCount The number of cell barcodes that do not encapsulate a cell in the experiment
#' @param nonEmptyCount The numbeYr of cell barcodes that encapsulate a cell in the experiment
#' @param showPlot Plot the poisson counts distribution using the inferred lambda parameter for the total number of cell barcodes.
#' 
#' @return A dataframe contain the input parameters, and the esimated doublet rate and lambda parameter.
#' @export
#' @import ggplot2
#' @examples estimateDoubletRatePoissonDistribution(10000,1000)
estimateDoubletRatePoissonDistribution<-function (emptyCount, nonEmptyCount, showPlot=TRUE) {
	
	#estimate the lambda to minimize difference between the observed counts of 0 and not 0 vs the counts from the poisson
	estimatePoisionGreaterThanZero<-function (lambda, fracNonZero) {
		#what is the cumulative probability of points with a count > 0?  IE: not empty
		estimatedFracNonZero=ppois(0, lambda, lower.tail = FALSE)
		#difference in observed and estimated, minimize the absolute difference.
		score=abs(fracNonZero-estimatedFracNonZero)
		#log_info(paste("Lambda [", lambda, "] score [", score, "]"), sep="")
		return(score)
	}
	
	total=emptyCount+nonEmptyCount
	fracNonZero=nonEmptyCount/total
	
	#an upper bound of 20 would have ~ 0 empty droplets.
	o=optim(par=1, fn=estimatePoisionGreaterThanZero, fracNonZero=fracNonZero, method="Brent", lower=0, upper=20)
	lambda=o$par
		
	#let's look at a maximum of 10 cells in a single droplet.
	success <- 0:10
	d=round (dpois(success, lambda=lambda)*total)
	df=data.frame(success=success, counts=d)
	fracDoublet=round (sum (df[df$success>1,]$counts)/sum (df[df$success>0,]$counts),3)
	summaryDF=data.frame(emptyCount=emptyCount, nonEmptyCount=nonEmptyCount, total=total, lambda=lambda, doubletRate=fracDoublet)
	
	#I hate making ggplot2 / R CMD check happy.
	counts<-NULL
	
	if (showPlot) {
		p=ggplot(df, aes(x=success, y=counts)) + 
			geom_bar(stat = "identity", fill="light blue") +
			geom_text(aes(label = counts), vjust = -0.5) +
			ggtitle(paste("Estimated Lambda [", round (lambda,3), "], Expected Fraction Doublets [", fracDoublet, "]")) +
			scale_x_continuous("", labels = as.character(success), breaks = success)
		print (p)
	}
	
	return (summaryDF)	
}

#' Estimate the number of empty cell barcodes in an experiment
#'
#' Manually draw thresholds around the ambient peak, and count the number of cells
#' 
#' @param cellMetricsFile A file containing one entry per cell.  The file is tab separated,
#' and the required columns are:
#' 1. cell_barcode - the string identifying the cell in the experiment.  Usually 12-16 bases long.
#' 2. num_transcripts - the total number of UMIs assigned to each cell.
#' 3. pct_intronic - the fraction of reads or UMIs that are intronic.
#' @param minPctIntronic The minimum PCT intronic for the ambient cell barcodes
#' @param maxPctIntronic The maximum PCT intronic for the ambient cell barcodes
#' @param minUMIs The minimum number of UMIs (in log10) for the ambient cell barcodes
#' @param maxUMIs The maximum number of UMIs (in log10) for the ambient cell barcodes
#'
#' @return The number of empty barcodes captured by the bounding box.
#' @export
estimateEmptyCount<-function (cellMetricsFile, minPctIntronic=0.28, maxPctIntronic=0.5, minUMIs=2.8, maxUMIs=3.3) {
	a=read.table(cellMetricsFile, header=T, stringsAsFactors = F, sep="\t")
	numEmpty=length(which(log10(a$num_transcripts)>minUMIs & log10(a$num_transcripts)<maxUMIs &
						  	a$pct_intronic<maxPctIntronic & a$pct_intronic>minPctIntronic))
	smoothScatter (log10(a$num_transcripts), a$pct_intronic, main=paste("Num empty cell barcodes [boxed selection] ", numEmpty))
	rect(xleft =minUMIs, ybottom= minPctIntronic, xright = maxUMIs, ytop = maxPctIntronic, col = "NA", lwd=2, border='red')  
	return (numEmpty)
}
	
	
