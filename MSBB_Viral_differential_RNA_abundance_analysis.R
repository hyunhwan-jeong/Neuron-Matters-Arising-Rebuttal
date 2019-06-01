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


## Differential expression and presence of viral gene feature level counts
viral_genomic_feature_DE_per_region <- vector("list", length(tissues))
names(viral_genomic_feature_DE_per_region) <- tissues

for(tissue_i in 1:length(tissues)){
  
  featureDEPerPathologySubset <- vector("list", length(MSBB_RNA_workspace$pathologySampleSets))
  names(featureDEPerPathologySubset) <- names(MSBB_RNA_workspace$pathologySampleSets)
  
  for(subset_i in 1:length(MSBB_RNA_workspace$pathologySampleSets)){
    
    caseMat <- MSBB_RNA_workspace$virusGenomicFeatureCounts[,intersect(MSBB_RNA_workspace$pathologySampleSets[[subset_i]], MSBB_RNA_workspace$metadata$Sample_ID[which(MSBB_RNA_workspace$metadata$Region == tissues[tissue_i])])]
    controlMat <- MSBB_RNA_workspace$virusGenomicFeatureCounts[,intersect(MSBB_RNA_workspace$controlSampleSets, MSBB_RNA_workspace$metadata$Sample_ID[which(MSBB_RNA_workspace$metadata$Region == tissues[tissue_i])])]
    
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
    lh <- limma.one.sided(fit, lower = TRUE)
    rh = 1 - lh
    
    viral_tT <- topTable(fit,number=Inf)
    
    viral_tT <- data.frame(sequence = MSBB_RNA_workspace$virus_name_2_accession_ID$VirusName[match(sapply(strsplit(rownames(viral_tT), split = "_"), function(x) paste(x[1:2], collapse = "_")), MSBB_RNA_workspace$virus_name_2_accession_ID$Accession)], name = rownames(viral_tT),viral_tT, row.names = NULL)
    
    viral_tT$downregulated_pvalue <- unname(lh)[match(viral_tT$name, names(lh))]
    viral_tT$upregulated_pvalue <- unname(rh)[match(viral_tT$name, names(rh))]
    
    featureDEPerPathologySubset[[subset_i]] <- list(viral_DE = viral_tT)
    
  }
  
  viral_genomic_feature_DE_per_region[[tissue_i]] <- featureDEPerPathologySubset
  
}



#Organize and annotate and output

# Aggregated --------------------------------------------------------------

multiregional_viral_level_DE <- ldply(lapply(virus_level_DE_per_region, function(y) ldply(lapply(y, function(x) x$viral_DE), rbind, .id = "Comparison")), rbind, .id = "Region") 
multiregional_viral_level_DE <- multiregional_viral_level_DE[order(multiregional_viral_level_DE$adj.P.Val, decreasing = FALSE),]

# Viral gene --------------------------------------------------------------

multiregional_viral_genomic_feature_level_DE <- ldply(lapply(viral_genomic_feature_DE_per_region, function(y) ldply(lapply(y, function(x) x$viral_DE), rbind, .id = "Comparison")), rbind, .id = "Region") 
multiregional_viral_genomic_feature_level_DE <- multiregional_viral_genomic_feature_level_DE[order(multiregional_viral_genomic_feature_level_DE$adj.P.Val, decreasing = FALSE),]

MSBB_RNA_DE_summary <- list(viral_level_DE = multiregional_viral_level_DE, 
                            viral_genomic_feature_level_DE = multiregional_viral_genomic_feature_level_DE
)

# Output results as heatmaps ----------------------------------------------
# Viral level heatmap -----------------------------------------------------
viralLevel_DE <- MSBB_RNA_DE_summary$viral_level_DE
viralLevel_DE$sequence<- gsub(viralLevel_DE$sequence, pattern = ", complete genome", replacement = "")

viralLevel_DE$tissue_comparison_compound_id <- paste(viralLevel_DE$Region, viralLevel_DE$Comparison, sep = "_")

viralLevel_DE$DE <- -log10(viralLevel_DE$P.Value) * sign(viralLevel_DE$logFC)
DE_expansive_matrix <- dcast(viralLevel_DE, sequence ~ tissue_comparison_compound_id, value.var = "DE", fill = 0)

rownames(DE_expansive_matrix) <- DE_expansive_matrix$sequence
DE_expansive_matrix <- DE_expansive_matrix[,-1] 
FDR_expansive_matrix <- dcast(viralLevel_DE, sequence ~ tissue_comparison_compound_id, value.var = "adj.P.Val", fill = NA)

rownames(FDR_expansive_matrix) <- FDR_expansive_matrix$sequence
FDR_expansive_matrix <- FDR_expansive_matrix[,-1]

Pval_expansive_matrix <- dcast(viralLevel_DE, sequence ~ tissue_comparison_compound_id, value.var = "P.Value", fill = NA)
rownames(Pval_expansive_matrix) <- Pval_expansive_matrix$sequence
Pval_expansive_matrix <- Pval_expansive_matrix[,-1]


#Collectively Significant
