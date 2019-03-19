options(stringsAsFactors = FALSE)

library(openxlsx)
library(edgeR)
library(limma)
library(plyr)
library(pheatmap)
library(reshape2)

# Load up workspace containing viral counts and metadata ------------------
load("data/MSBB_RNA_workspace.RData")
tissues <- c("BM_22", "BM_36", "BM_10", "BM_44")


## Differential expression of virus level abundance
virus_level_DE_per_region <- vector("list", length(tissues))
names(virus_level_DE_per_region) <- tissues

for(tissue_i in 1:length(tissues)){
  
  featureDEPerPathologySubset <- vector("list", length(MSBB_RNA_workspace$pathologySampleSets))
  names(featureDEPerPathologySubset) <- names(MSBB_RNA_workspace$pathologySampleSets)
  
  for(subset_i in 1:length(MSBB_RNA_workspace$pathologySampleSets)){
    
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
    
    virMat <- cbind(controlMat,caseMat)
    
    
    contrast.matrix<-makeContrasts(statusCase-statusControl,levels=designMat)##
    
    fullReadsPerSample <- as.numeric(MSBB_RNA_workspace$metadata$TotalReads[match(colnames(virMat), MSBB_RNA_workspace$metadata$Sample_ID)])# - as.numeric(U01_workspace$metadata$Mapped[match(colnames(virMat), U01_workspace$metadata$mapped_column_names)])
    
    sign_threshold_population <- 10
    
    virMat <- virMat[rowSums(sign(virMat)) > 0,]
    
    retainVir <- rowSums(virMat >= 2) >= sign_threshold_population
    
    dge <- DGEList(counts=virMat, lib.size = fullReadsPerSample)
    v <- voom(dge,designMat, lib.size = fullReadsPerSample, normalize.method = "quantile")
    fit <- lmFit(object = v[retainVir,], design = designMat)
    fit<- contrasts.fit(fit,contrast.matrix)
    fit <- eBayes(fit, robust = TRUE)
    
    viral_tT <- topTable(fit,number=Inf)
    
    viral_tT <- data.frame(sequence = MSBB_RNA_workspace$virus_name_2_accession_ID$VirusName[match(sapply(strsplit(rownames(viral_tT), split = "_"), function(x) paste(x[1:2], collapse = "_")), MSBB_RNA_workspace$virus_name_2_accession_ID$Accession)], name = rownames(viral_tT),viral_tT, row.names = NULL)
 
     featureDEPerPathologySubset[[subset_i]] <- viral_tT
    
  }
  
  virus_level_DE_per_region[[tissue_i]] <- featureDEPerPathologySubset
  
}

virus_level_DE_per_region_origin <- virus_level_DE_per_region
save(virus_level_DE_per_region_origin, file = "results/tmp.Rdata")
