options(stringsAsFactors = FALSE)

library(openxlsx)
library(edgeR)
library(limma)
library(plyr)
library(pheatmap)
library(reshape2)

source("helper_fxn.R")

# Load up workspace containing viral counts and metadata ------------------
load("data/MSBB_RNA_workspace.RData")

tissues <- c("BM_22", "BM_36", "BM_10", "BM_44")


while(dev.cur() != 1) dev.off()
pdf("results/zero_histogram.pdf", width = 12, height = 6)
layout(
  mat = matrix(c(9,10,1,2,3,4,5,6,7,8), nrow = 2, ncol = 5, byrow = F),
  heights = c(4,4),
  widths = c(1,4,4,4,4)
)

## Differential expression of virus level abundance
virus_level_DE_per_region <- vector("list", length(tissues))
names(virus_level_DE_per_region) <- tissues


for(tissue_i in 1:length(tissues)){
  
  featureDEPerPathologySubset <- vector("list", length(MSBB_RNA_workspace$pathologySampleSets))
  names(featureDEPerPathologySubset) <- names(MSBB_RNA_workspace$pathologySampleSets)
  
  for(subset_i in 3:length(MSBB_RNA_workspace$pathologySampleSets)){
    
    caseMat <- MSBB_RNA_workspace$virusLevelCounts[,intersect(MSBB_RNA_workspace$pathologySampleSets[[subset_i]], MSBB_RNA_workspace$metadata$Sample_ID[which(MSBB_RNA_workspace$metadata$Region == tissues[tissue_i])])]
    controlMat <- MSBB_RNA_workspace$virusLevelCounts[,intersect(MSBB_RNA_workspace$controlSampleSets, MSBB_RNA_workspace$metadata$Sample_ID[which(MSBB_RNA_workspace$metadata$Region == tissues[tissue_i])])]
    
    status <- factor(x = c(rep("Control", ncol(controlMat)), rep("Case", ncol(caseMat))), levels = c("Control","Case"))
    
    otherCovariates <- MSBB_RNA_workspace$metadata[match(c(colnames(controlMat),colnames(caseMat)), MSBB_RNA_workspace$metadata$Sample_ID),c("RIN", "AOD", "Race", "Sex", "Batch", "PMI")]
    otherCovariates$AOD <- gsub(x = otherCovariates$AOD, pattern = "+", replacement = "", fixed = TRUE)
    
    otherCovariates$AOD <- as.numeric(otherCovariates$AOD)
    otherCovariates$RIN <- as.numeric(otherCovariates$RIN)
    otherCovariates$PMI <- as.numeric(otherCovariates$PMI)
    
    otherCovariates$Race <- as.factor(otherCovariates$Race)
    otherCovariates$Sex <- as.factor(otherCovariates$Sex)
    otherCovariates$Batch <- as.factor(otherCovariates$Batch)
    
    comparison_annot <- data.frame(status = status, otherCovariates)
    
    designMat = model.matrix( ~ 0 + status + AOD + Race + RIN + Sex + Batch + PMI ,data = comparison_annot)
    
    sign_threshold_population <- 10
    virMat <- cbind(controlMat,caseMat)
    virMat <- virMat[rowSums(sign(virMat)) > 0,]
    retainVir <- rowSums(virMat >= 2) >= sign_threshold_population
    zero_prop <- sum(virMat==0) / sum(virMat>=0) * 100
    h = hist(unlist(log10(virMat+0.1)), plot = F)
    #h = hist(unlist(virMat), breaks = 1000, plot = F)
    h$density = h$counts/sum(h$counts) * 100
    cuts <- cut(h$breaks, c(-Inf,-1,Inf))
    
    
    main_str = sprintf("(%d viruses, %.0f%% zero)", nrow(virMat), zero_prop)
    par(mar = c(4, 4, 1, 2))
    
    plot(h, freq=F, ylim=c(0,100), ylab = "Percentage", xlab = "log10(count+0.1)", col = c("red", "grey")[cuts], main = NULL, xlim=c(-1,3))
    text(x = -0.3, y = 85, main_str, cex = 1.2, adj = 0, font = 2)
    text(x = 1.5, y = 95, tissues[tissue_i], adj = 1, font = 2, cex = 2 ) 

    virMat <- virMat[retainVir > 0,]
    
    # print(sum(virMat==0) / sum(virMat>=0))
    h = hist(unlist(log10(virMat+0.1)), plot = F)
    #h = hist(unlist(virMat), breaks = 1000, plot = F)
    h$density = h$counts/sum(h$counts) * 100
    
    cuts <- cut(h$breaks, c(-Inf,-1,Inf))
    zero_prop <- sum(virMat==0) / sum(virMat>=0) * 100
    main_str = sprintf("(%d viruses, %.0f%% zero)", nrow(virMat), zero_prop)
    #par(mar = c(5, 4, 0, 0))
    plot(h, freq=F, ylim=c(0,100), ylab = "Percentage", xlab = "log10(count+0.1)", col = c("red", "grey")[cuts], main = NULL, xlim=c(-1,3))
    text(x = -0.3, y = 85, main_str, cex = 1.2, adj = 0, font = 2)
    text(x = 1.5, y = 95, tissues[tissue_i], adj = 1, font = 2, cex = 2 ) 
    
  }
  
}


plot_title <- function(main) {
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, main,
       cex = 2, adj=0.5, srt = 90)
}

plot_title("Before Filtering")
plot_title("After Filtering")

dev.off()
