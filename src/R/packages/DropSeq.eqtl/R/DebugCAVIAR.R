#Debug eQTL CAVIAR AKAP10

# inSNPFile="/downloads/GWAS_correlation/CAVIAR_FILES/snpBlocks_daner_NGN2/28.Zscores.txt.gz"
# inCAVIARSNPFile="/downloads/GWAS_correlation/CAVIAR_FILES/snpBlocks_daner_NGN2/out/28.eQTL.Z.txt"
# eQTLLargeWindowFile="/downloads/GWAS_correlation/CAVIAR_FILES/d14-42_NGN2.maf_0.05_cisDist_250kb.eQTL_results.txt"
#
# caviar95File="/downloads/GWAS_correlation/CAVIAR_FILES/snpBlocks_daner_NGN2/out/28.eQTL_set"
# eQTL_LDFile="/downloads/GWAS_correlation/CAVIAR_FILES/snpBlocks_daner_NGN2/out/28.eQTL.LD.txt"
# eQTLPosteriors="/downloads/GWAS_correlation/CAVIAR_FILES/snpBlocks_daner_NGN2/out/28.eQTL_post"

# testSNPDFEQTL<-function (inSNPFile, caviar95File, eQTLPosteriors) {
#     df=read.table(inSNPFile, header=T, stringsAsFactors = F, sep="\t")
#     df=df[!is.na(df$eQTL_P),]
#     geneName=unique (df$gene)
#     a=fread(eQTLLargeWindowFile, header=T, stringsAsFactors = F, sep="\t")
#     a=a[a$gene==geneName,]
#
#     idx=match(df$snp, a$SNP)
#     which(is.na(idx))
#     which (df$eQTL_P!=a[idx,]$`p-value`)
#     which (df$eQTL_beta!=a[idx,]$beta)
#
#     #p-value and beta not mangled.
#     #p-value and z-score scale appropriately.
#     plot (-log10(df$eQTL_P), abs(df$eQTL_Z))
#
#     #are there snps missing in the input to caviar?
#     inCav=read.table(inCAVIARSNPFile, header=F, stringsAsFactors = F)
#     dfMiss=df[match(setdiff(df$snp, inCav$V1), df$snp),]
#     #yes, because there's no LD score in 1kg.
#     ldBlock=read.table(eQTL_LDFile, header=F, stringsAsFactors=F)
#     intersect (colnames(ldBlock), dfMiss$snp)
#
#     intersect (inCav$V1, colnames(ldBlock))
#
#     #test the caviar SNP input
#     b=read.table(inCAVIARSNPFile, header=F, stringsAsFactors = F)
#     idx=match(b$V1, df$snp)
#     df=df[idx,]
#     b$V2==df$eQTL_Z
#     setdiff(df$snp, b$V1)
#
#     #test the posteriors probs...
#     z=read.table(eQTLPosteriors, header=T)
#     idx=match(z$SNP_ID, df$snp)
#     z$pos=df[idx,]$pos
#     par (mfrow=c(2,1), mar=c(2,4,0,2))
#     idx2=match(z[z$Causal_Post._Prob.>0.1,]$SNP_ID, df$snp)
#     plot (df$pos, -log10(df$eQTL_P))
#     points(df[idx2,]$pos, -log10(df[idx2,]$eQTL_P), col="red", pch=16)
#     plot (z$pos, z$Causal_Post._Prob.)
#
#     #the results are the same.
# }

# testCAVIAReQTLResult<-function (inSNPFile, caviar95File) {
#     df=read.table(inSNPFile, header=T, stringsAsFactors = F, sep="\t")
#     df=df[!is.na(df$eQTL_P),]
#     a=read.table(caviar95File, header=F,stringsAsFactors = F)
#     idx=match(a$V1, df$snp)
#     plot (df$pos, -log10(df$eQTL_P))
#     points (df[idx,]$pos, -log10(df[idx,]$eQTL_P), col='red')
#
#     z=read.table(eQTLPosteriors, header=T)
#     idx=match(df$snp, z$SNP_ID)
#     z=z[idx,]
#     plot (df$pos, z$Causal_Post._Prob.)
#
#     library(LDheatmap)
#     ldBlock=read.table(eQTL_LDFile, header=F, stringsAsFactors=F)
#     LDheatmap(gdat=abs(as.matrix(ldBlock[1:100,1:100])), flip=T, add.map=F, newpage=F, title="")
#     LDheatmap(gdat=abs(as.matrix(ldBlock)), flip=T, add.map=F, newpage=F, title="")
#     LDheatmap(gdat=as.matrix(ldBlock), flip=T, add.map=F, newpage=F, title="")
#
#     hist (as.matrix(ldBlock), xlab="LD", main ="Genotype correlation")
#     hist (as.matrix(abs(ldBlock)), xlab="LD (absolute value)", main ="Genotype correlation")
#
#
# }
