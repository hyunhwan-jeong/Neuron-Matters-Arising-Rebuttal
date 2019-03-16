library(edgeR)
library(tidyverse)
source("../scripts/R/helper_fxn.R")
load("data/MSBB_RNA_workspace.RData")

tissues <- c("BM_22", "BM_36", "BM_10", "BM_44")

virus_level_DE_per_region <- vector("list", length(tissues))
names(virus_level_DE_per_region) <- tissues
interests <- c("NC_001716.2_region_1_153080__ID=id0", "NC_001664.2_region_1_159322__ID=id0")

df_log2CPM <- data.frame()
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
    
    log2CPM <- cpm(virMat, lib.size = fullReadsPerSample, log=T, prior.count = 0.1)
    
    HHV6A <- data.frame(id = "HHV6A", log2CPM = log2CPM[interests[2],], comparison_annot) %>% 
      mutate(tissue = tissues[tissue_i])
    HHV7 <- data.frame(id = "HHV7", log2CPM = log2CPM[interests[1],], comparison_annot) %>% 
      mutate(tissue = tissues[tissue_i])
    cat(tissues[tissue_i], table(comparison_annot$status), "\n")
    HHV6A_t <- HHV6A %>% filter(status == "Case") %>% mutate(status = names(MSBB_RNA_workspace$pathologySampleSets)[subset_i])
    HHV7_t <- HHV7 %>% filter(status == "Case")%>% mutate(status = names(MSBB_RNA_workspace$pathologySampleSets)[subset_i])
    print(nrow(HHV6A_t))
    print(nrow(HHV7_t))
    df_log2CPM <- rbind(df_log2CPM, HHV6A_t, HHV7_t)
    if(subset_i == 1) {
      df_log2CPM <- rbind(df_log2CPM, 
                          HHV6A %>% filter(status == "Control"),
                          HHV7 %>% filter(status == "Control"))
    }
    
  }

}

# https://github.com/tidyverse/ggplot2/issues/2499

df_log2CPM <- as_tibble(df_log2CPM)

df_log2CPM %>% filter(id == "HHV6A", tissue=="BM_10", status=="Control") %>% 
  ggplot(aes(x=log2CPM, y=stat(width*density))) +
  geom_histogram() 

p1 <- df_log2CPM %>% filter(id == "HHV6A") %>% 
  ggplot(aes(x=log2CPM, y=stat(width*density))) +
  geom_histogram() +
  facet_grid(status~tissue) + ggtitle("HHV6A") + theme_bw()

p2 <- df_log2CPM %>% filter(id == "HHV7") %>% 
  ggplot(aes(x=log2CPM, y=stat(width*density))) +
  geom_histogram() +
  facet_grid(status~tissue) + ggtitle("HHV7") + theme_bw()

library(cowplot)

plot_grid(p1, p2)
