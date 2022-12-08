#   Nextflow workflow performing an ARN-seq Analysis created by :
# - Ambre Baumann (https://github.com/ambrebaumann)
# - Lindsay Goulet (https://github.com/Lindsay-Goulet)
# - George Marchment (https://github.com/George-Marchment)
# - Clémence Sebe (https://github.com/ClemenceS)

#This workflow was created in the context of the 2022 Hackathon project (master M2 AMI2B) 
#This project was surpervised by Frédéric Lemoine & Thomas Cokelaer

# -------------------------------------- 
# Differential Analyse

# -------------------------------------- 
library(DESeq2)
library(ggplot2)
# -------------------------------------- 


#Load the files (countingReads and metadata (give for each read if it is a mutant or not))
args = commandArgs(trailingOnly=TRUE)  
countingReads = args[1]
metadata = args[2]

countData = read.table(header=TRUE, row.names = 1, countingReads)
countData = countData[,6:13]
metadata  = read.table(header=TRUE, sep="\t", metadata)

# -------------------------------------- 

#Initialisation of a vector which contains the information if a read is a mutant 
cond = c()
for (ind in colnames(countData)){
  cond = c(cond, metadata[which(metadata[,1]==ind),2] == 1)
}
cond = factor(cond)

#Revome the gene's lines with all its counts egal to 0 
countData = countData[-which(rowSums(countData) == 0),]

# -------------------------------------- 

# ------ First analyse : Before Cleaning ------ 

#Initialisation of the matrix in the 'good' format for running the DESeq analyse
dds_before = DESeqDataSetFromMatrix(countData=countData, colData=DataFrame(cond), design=~cond)
#Run the analyse
dds_before = DESeq(dds_before)
#Get the results
res_before = results(dds_before)

# -- Plot save in a pdf --
pdf('beforeCleaning.pdf')

#Dispersion plot 
plotDispEsts(dds_before, main="PlotDispEsts before cleaning")

#Histogram of the p-values
hist(res_before$pvalue, main='Histogram of pvalues before cleaning')

#Perform a log transformation to the data
rld_before = rlog(dds_before)
#Plot the PCA 
plotPCA(rld_before, intgroup=c("cond")) + geom_text(aes(label=name),vjust=2) + ggtitle("PCA (log2 scale) before cleaning with mutant = TRUE")

#Plot MA
plotMA(res_before, main='PlotMa before cleaning and with alpha = 0.05',alpha=0.05)

#PlotCounts for the gene with the smallest p-value
title = paste("PlotCounts - Min padj : ", row.names(res_before[which.min(res_before$padj),]), " : ", res_before[which.min(res_before$padj),]$padj, '\n')
plotCounts(dds_before,which.min(res_before$padj),intgroup = "cond",  normalized = TRUE, main= title)

#PlotCounts for the gene with the biggest p-value
title = paste("PlotCounts - Max padj : ", row.names(res_before[which.max(res_before$padj),]), " : ", res_before[which.max(res_before$padj),]$padj, '\n')
plotCounts(dds_before,which.max(res_before$padj),intgroup = "cond",  normalized = TRUE, main = title)

dev.off()

# -- Save in a txt file the genes with a p-value ajusted less than 0.05 --
DEgenes_before = res_before[which(res_before$padj<0.05),]
write.table(DEgenes_before, "DEgenes_beforeCleaning.txt", sep="\t")

# -------------------------------------- 
# Apply a Cleaning 'function' : Keep all genes with at least 5 reads in both conditions
# -------------------------------------- 
countData = countData[which(rowSums(countData) > 5),]

# ------ Second analyse : After Cleaning ------ 

#Initialisation of the matrix in the 'good' format for running the DESeq analyse
dds = DESeqDataSetFromMatrix(countData=countData, colData=DataFrame(cond), design=~cond)
#Run the analyse
dds = DESeq(dds)
#Get the results
res = results(dds)

# -- Plot save in a pdf --
pdf("afterCleaning.pdf")

#Dispersion plot 
plotDispEsts(dds, main="PlotDispEsts")

#Histogram of the p-values
hist(res$pvalue, main='Histogram of pvalues')

#Perform a log transformation to the data
rld = rlog(dds)

#Plot the PCA 
plotPCA(rld, intgroup=c("cond")) + geom_text(aes(label=name),vjust=2)  + ggtitle("PCA (log2 scale) with mutant = TRUE")


#Plot MA
plotMA(res, main='PlotMa with alpha = 0.05',alpha=0.05)


#PlotCounts for the gene with the smallest p-value
title = paste("PlotCounts - Min padj : ", row.names(res[which.min(res$padj),]), " : ", res[which.min(res$padj),]$padj, '\n')
plotCounts(dds,which.min(res$padj),intgroup = "cond",  normalized = TRUE, main=title)

#PlotCounts for the gene with the biggest p-value
title = paste("PlotCounts - Max padj : ", row.names(res[which.max(res$padj),]), " : ", res[which.max(res$padj),]$padj, '\n')
plotCounts(dds,which.max(res$padj),intgroup = "cond",  normalized = TRUE, main=title)

dev.off()

# -- Save in a txt file the genes with a p-value ajusted less than 0.05 --
DEgenes = res[which(res$padj<0.05),]
write.table(DEgenes, "DEgenes.txt", sep="\t")


# ------ Focus on the SF3B1 gene ------ 

# -- Plot save in a pdf --
pdf("SF3B1.pdf")

#PlotCounts for SF3B1 before cleaning
plotCounts(dds_before,"ENSG00000115524",intgroup = "cond",  normalized = TRUE, main= 'PlotCounts ENSG00000115524 before cleaning')

#PlotCounts for SF3B1 after cleaning
plotCounts(dds,"ENSG00000115524",intgroup = "cond",  normalized = TRUE, main= 'PlotCounts ENSG00000115524 after cleaning')

dev.off()

# -- Save in a txt file some statistics --
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