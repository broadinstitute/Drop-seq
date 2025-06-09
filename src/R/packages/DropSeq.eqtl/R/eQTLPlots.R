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

#' Plot diversity of UMIs/Donors in an eQTL experiment.
#' @param metaCellFile A single file of metacells.
#' @param outPDF The PDF file to write to, or NULL.
#' @import vegan
#' @export
plotDiversity <-function (metaCellFile, outPDF=NULL) {
    a=read.table(metaCellFile, header=T, stringsAsFactors = F, check.names=F)
    d=colSums(a[,-1])
    d=d[order(d, decreasing = T)]
    div=vegan::diversity(d)
    expName=sub (".meta_cell.expression.txt", "", basename (metaCellFile))
    strTitle=paste(expName, "\nDistribution of UMIs across donors\n",length(d), " donors; ", round (sum (d)/1e6,1), "M UMIs; SW Div: ", sprintf("%.2f",round(div,2)), "; SW Eq: ", sprintf("%.2f", round (div/log(length(d)),2)), sep="")
    if (!is.null(outPDF)) pdf(outPDF, width=8, height=6)
    par(mar=c(8,4,4,2))
    barplot ((d/1e6), las=3, cex.axis=0.5, cex.names=0.75, cex.lab=1, ylab="Total UMIs [millions]", names.arg = names(d), main=strTitle, col="light blue")
    par(mar=c(4,4,4,4))
    if (!is.null(outPDF)) dev.off()
}
