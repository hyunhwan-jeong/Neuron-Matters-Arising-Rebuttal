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
    
    sign_threshold_population <- 10
    
    virMat <- virMat[rowSums(sign(virMat)) > 0,]
    
    retainVir <- rowSums(virMat >= 2) >= sign_threshold_population
    
    dge <- DGEList(counts=virMat[retainVir,], lib.size = fullReadsPerSample)
    v <- voom(dge,designMat, lib.size = fullReadsPerSample, normalize.method = "quantile")
    fit <- lmFit(object = v, design = designMat)
    fit<- contrasts.fit(fit,contrast.matrix)
    fit <- eBayes(fit, robust = TRUE)
    lh <- limma.one.sided(fit, lower = TRUE)
    rh = 1 - lh
    
    viral_tT <- topTable(fit,number=Inf)
    
    viral_tT <- data.frame(sequence = MSBB_RNA_workspace$virus_name_2_accession_ID$VirusName[match(sapply(strsplit(rownames(viral_tT), split = "_"), function(x) paste(x[1:2], collapse = "_")), MSBB_RNA_workspace$virus_name_2_accession_ID$Accession)], name = rownames(viral_tT),viral_tT, row.names = NULL)
    
    viral_tT$downregulated_pvalue <- unname(lh)[match(viral_tT$name, names(lh))]
    viral_tT$upregulated_pvalue <- unname(rh)[match(viral_tT$name, names(rh))]
    
    featureDEPerPathologySubset[[subset_i]] <- list(viral_DE = viral_tT)
    
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

df_ret %>% select(tissue, trait, sequence, logFC, AveExpr, P.Value, adj.P.Val) %>% View
df_ret %>% select(tissue, trait, sequence, logFC, AveExpr, P.Value, adj.P.Val) %>% write_csv("results/voom_limma_prenorm.csv")

df_ret %>% select(tissue, trait, sequence, logFC, AveExpr, P.Value, adj.P.Val) %>%
  ggplot(aes(x=trait, y=-log10(adj.P.Val))) +
  geom_bar(stat = "identity", aes(fill=trait), width=0.5) + 
  geom_hline(yintercept = -log10(0.1)) +
  facet_grid(tissue~sequence)

df_ret %>% select(tissue, trait, sequence, logFC, AveExpr, P.Value, adj.P.Val) %>%
  ggplot(aes(x=trait, y=logFC)) +
  geom_bar(stat = "identity", aes(fill=trait), width=0.5) + 
  geom_hline(yintercept = 0) +
  facet_grid(tissue~sequence)
