#metaCellFile="/downloads/meta_cell/d21NGN2Live.meta_cell.expression.txt"

#' Downsample the UMIs of donors by some fraction of their UMIs.
#' @param metaCellFile A meta-cell file
#' @param fractionUMIs the fraction of UMIs to retain
#' @param outFile The output file to write results to.  If null return the data frame.
#' @return The meta cell matrix after downsampling.
#' @export
downsampleMetaCells<-function (metaCellFile, fractionUMIs, outFile=NULL) {
    a=read.table(metaCellFile, header=T, stringsAsFactors = F, sep="\t", check.names = F)
    donors=colnames(a)[-1]

    #downsamples some number of observations by the fraction to keep.
    downsampleUMI<-function (count, fractionToKeep) {
        idx=runif(count)
        return(length(which(idx<=fractionToKeep)))
    }

    downsampleDonor<-function (donor, fractionToKeep,a) {
        x=a[,donor]
        numUMIsToKeep=sum(x)*fractionToKeep
        if (numUMIsToKeep>=sum(x)) return (x)
        xD=sapply(x, downsampleUMI, fractionToKeep)
        return (xD)
    }

    downsampleDonors<-function (donors, fractionToKeep, a) {
        aD=lapply(donors, downsampleDonor, fractionToKeep, a)
        aD=do.call(cbind, aD)
        rownames(aD)=rownames(a)
        colnames(aD)=colnames(a)[-1]
        result=data.frame(GENE=a$GENE, aD, stringsAsFactors = F)
        return (result)
    }
    #a, downsampled.
    aD=downsampleDonors(donors, fractionUMIs, a)
    if (!is.null(outFile)) write.table(aD, outFile, row.names = F, col.names = T, quote=F, sep="\t")
}
