
inMetaCellFile="/broad/mccarroll/nemesh/eQTL_Runs/05-29-2020/d5QueenBasic.meta_cell.expression.txt"
#inKinshipFile="/broad/mccarroll/dropulation_census_metadata/COVARIATES/CIRMw1w2w3.no_twins.kin"
inKinshipFile="/broad/mccarroll/nemesh/dropulation_paper/kinship/CIRM.kin"

outKinshipMatrix="/broad/mccarroll/dropulation_census_metadata/COVARIATES/CIRMw1w2w3.kinship_matrix.txt"

getMatrix<-function (inMetaCellFile, inKinshipFile, outKinshipMatrix) {
    a=read.table(inMetaCellFile, header=T, stringsAsFactors = F, sep="\t", check.names = F, nrows = 1)
    b=read.table(inKinshipFile, header=T, stringsAsFactors = F, sep="\t", check.names = F)

    donorOrder=colnames(a)[-1]

    getOne<-function (donor1, donor2, b) {
        idx=c(which(b$ID1==donor1 & b$ID2==donor2), which(b$ID1==donor2 & b$ID2==donor1))
        #donor pair doesn't exist.
        if (length(idx)==0) return (0)
        kinship=b[idx,]$Kinship
        if (kinship<0) kinship=0
        return (kinship)
    }
    #initialize to 0.
    outMatrix=matrix(0, nrow=length(donorOrder), ncol=length(donorOrder), dimnames=list(donorOrder, donorOrder))
    #rowIdx=1;colIdx=2
    for (rowIdx in 1:dim(outMatrix)[1]) {
        donor1=rownames(outMatrix)[rowIdx]
        for (colIdx in 1:dim(outMatrix)[2]) {
            donor2=colnames(outMatrix)[colIdx]
            outMatrix[rowIdx,colIdx]=getOne(donor1, donor2, b)
        }
    }

    #set diagonal to 1.  Self testing doesn't occur in kinship file.
    diag(outMatrix)=0.5

    #matrix should be 2*kinship
    outMatrix=outMatrix*2

    heatmap(outMatrix, symm=T, cexCol=0.0001)
    write.table(outMatrix, outKinshipMatrix, row.names = T, col.names = T, quote=F, sep="\t")


}

compareTwoKingRuns<-function () {

    getTruncatedID<-function (x) {
        strsplit(x, "_", fixed=T)[[1]][1]
    }
    a=read.table("/broad/mccarroll/nemesh/dropulation_paper/kinship/CIRM.kin", header=T, stringsAsFactors = F, sep="\t", check.names = F)
    b=read.table("/broad/mccarroll/dropulation_census_metadata/COVARIATES/cirm_core24v10_core24v12.GrCh37.merged.1507.kin", header=T, stringsAsFactors = F, sep="\t", check.names = F)

    a$ID1=sapply(a$ID1, getTruncatedID)
    a$ID2=sapply(a$ID1, getTruncatedID)

    b$ID1=sapply(b$ID1, getTruncatedID)
    b$ID2=sapply(b$ID1, getTruncatedID)



}
