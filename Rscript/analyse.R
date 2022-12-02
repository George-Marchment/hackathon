library(DESeq2)
library(ggplot2)

#recup des fichiers
args = commandArgs(trailingOnly=TRUE)  
countingReads = args[1]
metadata = args[2]

#lecture
countData = read.table(header=TRUE, row.names = 1, countingReads)
countData = countData[,6:13]
metadata  = read.table(header=TRUE, sep="\t", metadata)

#initiaiton pour dire si mutant ou non
cond = c()
for (ind in colnames(countData)){
  cond = c(cond, metadata[which(metadata[,1]==ind),2] == 1)
}
cond = factor(cond)

#on enleve les lignes totalement egales à 0
countData = countData[-which(rowSums(countData) == 0),]

###########################################################################

#premiere analyse before cleaning
dds_before = DESeqDataSetFromMatrix(countData=countData, colData=DataFrame(cond), design=~cond)
dds_before = DESeq(dds_before)
res_before = results(dds_before)

#save in a dpf all before cleaning
pdf('beforeCleaning.pdf')
plotDispEsts(dds_before, main="PlotDispEsts before cleaning")

hist(res_before$pvalue, main='Histogram of pvalues before cleaning')

rld_before = rlog(dds_before)
plotPCA(rld_before, intgroup=c("cond")) + geom_text(aes(label=name),vjust=2) + ggtitle("PCA (log2 scale) before cleaning with mutant = TRUE")

plotMA(res_before, main='PlotMa before cleaning and with alpha = 0.05',alpha=0.05)

title = paste("PlotCounts - Min padj : ", row.names(res_before[which.min(res_before$padj),]), " : ", res_before[which.min(res_before$padj),]$padj, '\n')
plotCounts(dds_before,which.min(res_before$padj),intgroup = "cond",  normalized = TRUE, main= title)
title = paste("PlotCounts - Max padj : ", row.names(res_before[which.max(res_before$padj),]), " : ", res_before[which.max(res_before$padj),]$padj, '\n')
plotCounts(dds_before,which.max(res_before$padj),intgroup = "cond",  normalized = TRUE, main = title)
dev.off()

DEgenes_before = res_before[which(res_before$padj<0.05),]
write.table(DEgenes_before, "DEgenes_beforeCleaning.txt", sep="\t")

###########################################################################

##Cleaning 
mutant = metadata[which(metadata[,2]==1),1]
nonMutant = metadata[-which(metadata[,2]==1),1]
for (indMutant in mutant){
  for (indNonMutant in nonMutant) {
    countData = countData[countData[,indMutant]!=0 | countData[,indNonMutant]!=0,]
  }
}

# Analyse différentielle
dds = DESeqDataSetFromMatrix(countData=countData, colData=DataFrame(cond), design=~cond)
dds = DESeq(dds)
res = results(dds)

pdf("afterCleaning.pdf")
plotDispEsts(dds, main="PlotDispEsts")

hist(res$pvalue, main='Histogram of pvalues')

rld = rlog(dds)
plotPCA(rld, intgroup=c("cond")) + geom_text(aes(label=name),vjust=2)  + ggtitle("PCA (log2 scale) with mutant = TRUE")

plotMA(res, main='PlotMa with alpha = 0.05',alpha=0.05)

title = paste("PlotCounts - Min padj : ", row.names(res[which.min(res$padj),]), " : ", res[which.min(res$padj),]$padj, '\n')
plotCounts(dds,which.min(res$padj),intgroup = "cond",  normalized = TRUE, main=title)
title = paste("PlotCounts - Min padj : ", row.names(res[which.max(res$padj),]), " : ", res[which.max(res$padj),]$padj, '\n')
plotCounts(dds,which.max(res$padj),intgroup = "cond",  normalized = TRUE, main=title)
dev.off()

DEgenes = res[which(res$padj<0.05),]
write.table(DEgenes, "DEgenes.txt", sep="\t")

pdf("SF3B1.pdf")
plotCounts(dds_before,"ENSG00000115524",intgroup = "cond",  normalized = TRUE, main= 'PlotCounts ENSG00000115524 before cleaning')
plotCounts(dds,"ENSG00000115524",intgroup = "cond",  normalized = TRUE, main= 'PlotCounts ENSG00000115524 after cleaning')
dev.off()

#save stat in a txt
sink("stat.txt")
cat("The results were obtained with alpha = 0.05 and used the adjusted p-value (padj)")
cat("\n\n\n")
before = paste ("Numbers of DE genes before cleaning : ", dim(DEgenes_before)[1])
cat(before)
cat("\n")
mini = paste("Min padj : ", row.names(res_before[which.min(res_before$padj),]) , " : ", res_before[which.min(res_before$padj),]$padj, '\n')
cat(mini)
maxi = paste("Max padj : ", row.names(res_before[which.max(res_before$padj),]) , " : ", res_before[which.max(res_before$padj),]$padj, '\n')
cat(maxi)
cat("\n")
after = paste ("Numbers of DE genes after cleaning : ", dim(DEgenes)[1])
cat(after)
cat("\n")
mini = paste("Min padj : ", row.names(res[which.min(res$padj),]) , " : ", res[which.min(res$padj),]$padj, '\n')
cat(mini)
maxi = paste("Max padj : ", row.names(res[which.max(res$padj),]) , " : ", res[which.max(res$padj),]$padj, '\n')
cat(maxi)
cat("\n\n")
cat("Gene SF3B1, ENSG00000115524 : \n")
sf3b1 = paste("padj before cleaning : ", res_before['ENSG00000227232',]$padj)
cat(sf3b1)
cat("\n")
sf3b1 = paste("padj after cleaning : ", res['ENSG00000227232',]$padj)
cat(sf3b1)
sink()
