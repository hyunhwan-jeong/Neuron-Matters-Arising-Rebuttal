library(DESeq2)
library(edgeR)
library(tidyverse)
source("../scripts/R/helper_fxn.R")
load("data/MSBB_RNA_workspace.RData")

tissues <- c("BM_22", "BM_36", "BM_10", "BM_44")

virus_level_DE_per_region <- vector("list", length(tissues))
names(virus_level_DE_per_region) <- tissues
interests <- c("NC_001716.2_region_1_153080__ID=id0", "NC_001664.2_region_1_159322__ID=id0")

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
    
    df_tmp <- data.frame()    
    viral_tT <- run_deseq2(virMat, fullReadsPerSample, comparison_annot)
    viral_tT <- data.frame(sequence = MSBB_RNA_workspace$virus_name_2_accession_ID$VirusName[match(sapply(strsplit(rownames(viral_tT), split = "_"), function(x) paste(x[1:2], collapse = "_")), MSBB_RNA_workspace$virus_name_2_accession_ID$Accession)], name = rownames(viral_tT),viral_tT, row.names = NULL)
    viral_tT <- viral_tT %>% 
      select(sequence, name, logFC = log2FoldChange, FDR = padj) %>% 
      mutate(method = "DESeq2noFilter")
    df_tmp <- rbind(df_tmp, viral_tT)
    
    viral_tT <- run_deseq2(virMat, fullReadsPerSample, comparison_annot, filter = F)
    viral_tT <- data.frame(sequence = MSBB_RNA_workspace$virus_name_2_accession_ID$VirusName[match(sapply(strsplit(rownames(viral_tT), split = "_"), function(x) paste(x[1:2], collapse = "_")), MSBB_RNA_workspace$virus_name_2_accession_ID$Accession)], name = rownames(viral_tT),viral_tT, row.names = NULL)
    viral_tT <- viral_tT %>% 
      select(sequence, name, logFC = log2FoldChange, FDR = padj) %>% 
      mutate(method = "DESeq2wtFilter")
    df_tmp <- rbind(df_tmp, viral_tT)
    
    viral_tT <- run_voom(virMat, fullReadsPerSample, comparison_annot)
    viral_tT <- data.frame(sequence = MSBB_RNA_workspace$virus_name_2_accession_ID$VirusName[match(sapply(strsplit(rownames(viral_tT), split = "_"), function(x) paste(x[1:2], collapse = "_")), MSBB_RNA_workspace$virus_name_2_accession_ID$Accession)], name = rownames(viral_tT),viral_tT, row.names = NULL)
    viral_tT <- viral_tT %>% 
      select(sequence, name, logFC = logFC, FDR = adj.P.Val) %>% 
      mutate(method = "VoomLimmanoFilter")
    df_tmp <- rbind(df_tmp, viral_tT)
    
    viral_tT <- run_voom(virMat, fullReadsPerSample, comparison_annot, filter = F)
    viral_tT <- data.frame(sequence = MSBB_RNA_workspace$virus_name_2_accession_ID$VirusName[match(sapply(strsplit(rownames(viral_tT), split = "_"), function(x) paste(x[1:2], collapse = "_")), MSBB_RNA_workspace$virus_name_2_accession_ID$Accession)], name = rownames(viral_tT),viral_tT, row.names = NULL)
    viral_tT <- viral_tT %>% 
      select(sequence, name, logFC = logFC, FDR = adj.P.Val) %>% 
      mutate(method = "VoomLimmawtFilter")
    df_tmp <- rbind(df_tmp, viral_tT)

    viral_tT <- run_edgeR(virMat, fullReadsPerSample, comparison_annot)
    viral_tT <- data.frame(sequence = MSBB_RNA_workspace$virus_name_2_accession_ID$VirusName[match(sapply(strsplit(rownames(viral_tT), split = "_"), function(x) paste(x[1:2], collapse = "_")), MSBB_RNA_workspace$virus_name_2_accession_ID$Accession)], name = rownames(viral_tT),viral_tT, row.names = NULL)
    viral_tT <- viral_tT %>% 
      select(sequence, name, logFC, FDR) %>% 
      mutate(method = "EdgeRnoFilter")
    df_tmp <- rbind(df_tmp, viral_tT)
    
    viral_tT <- run_edgeR(virMat, fullReadsPerSample, comparison_annot, filter = F)
    viral_tT <- data.frame(sequence = MSBB_RNA_workspace$virus_name_2_accession_ID$VirusName[match(sapply(strsplit(rownames(viral_tT), split = "_"), function(x) paste(x[1:2], collapse = "_")), MSBB_RNA_workspace$virus_name_2_accession_ID$Accession)], name = rownames(viral_tT),viral_tT, row.names = NULL)
    viral_tT <- viral_tT %>% 
      select(sequence, name, logFC, FDR) %>% 
      mutate(method = "EdgeRwithFilter")
    df_tmp <- rbind(df_tmp, viral_tT)
    
    featureDEPerPathologySubset[[subset_i]] <- list(viral_DE = df_tmp)
    
  }
  
  virus_level_DE_per_region[[tissue_i]] <- featureDEPerPathologySubset
  
}

df_ret <- data.frame()

for(i in tissues) {
  traits <- names(virus_level_DE_per_region[[i]])
  for(j in traits) {
    df_ij <- virus_level_DE_per_region[[i]][[j]]$viral_DE
    df_interest <- df_ij %>% dplyr::filter(name %in% interests) 
    df_interest$tissue <- i
    df_interest$trait <- j
    df_ret <- rbind(df_ret, df_interest)
  }
}


p1 <- df_ret %>%
  ggplot(aes(x=trait, y=-log10(FDR))) +
  #geom_line(stat = "identity", aes(group=method, color=method)) + 
  geom_bar(stat = "identity", aes(fill=method), width=0.5, position = "dodge") + 
  geom_hline(yintercept = -log10(0.1)) +
  facet_grid(tissue~sequence) 


p2 <- df_ret %>% 
  ggplot(aes(x=trait, y=logFC)) +
  geom_bar(stat = "identity", aes(fill=method), width=0.5, position = "dodge") + 
  geom_hline(yintercept = 0) +
  facet_grid(tissue~sequence)

p1
p2
