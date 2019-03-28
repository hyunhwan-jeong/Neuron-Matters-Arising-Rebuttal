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
  
  Mat <- MSBB_RNA_workspace$virusLevelCounts[,MSBB_RNA_workspace$metadata$Sample_ID[which(MSBB_RNA_workspace$metadata$Region == tissues[tissue_i])]]
  
  
  otherCovariates <- MSBB_RNA_workspace$metadata[match(colnames(Mat), MSBB_RNA_workspace$metadata$Sample_ID),c("RIN", "AOD", "Race", "Sex", "Batch", "PMI", "NP.1")]
  otherCovariates$AOD <- gsub(x = otherCovariates$AOD, pattern = "+", replacement = "", fixed = TRUE)
  
  otherCovariates$AOD <- as.numeric(otherCovariates$AOD)
  otherCovariates$RIN <- as.numeric(otherCovariates$RIN)
  otherCovariates$PMI <- as.numeric(otherCovariates$PMI)
  
  otherCovariates$Race <- as.factor(otherCovariates$Race)
  otherCovariates$Sex <- as.factor(otherCovariates$Sex)
  otherCovariates$Batch <- as.factor(otherCovariates$Batch)
  otherCovariates$NP.1 <- as.factor(otherCovariates$NP.1)
  
  comparison_annot <- data.frame(otherCovariates)
  
  designMat = model.matrix( ~ 0 + NP.1 + AOD + Race + RIN + Sex + Batch + PMI ,data = comparison_annot)
  
  virMat <- Mat
  
  contrast.matrix <- makeContrasts(AD_Definite = NP.12 - NP.11,
                                   AD_Likely = (NP.12+NP.13)/2 - NP.11,
                                   AD_Possible = (NP.12+NP.13+NP.14)/3 - NP.11,
                                   levels=designMat)
  
  
  fullReadsPerSample <- as.numeric(MSBB_RNA_workspace$metadata$TotalReads[match(colnames(virMat), MSBB_RNA_workspace$metadata$Sample_ID)])
  
  sign_threshold_population <- 10
  
  virMat <- virMat[rowSums(sign(virMat)) > 0,]
  
  retainVir <- rowSums(virMat >= 2) >= sign_threshold_population
  
  dge <- DGEList(counts=virMat, lib.size = fullReadsPerSample)
  v <- voom(dge,designMat, lib.size = fullReadsPerSample, normalize.method = "quantile")
  fit <- lmFit(object = v[retainVir,], design = designMat)
  fit<- contrasts.fit(fit,contrast.matrix)
  fit <- eBayes(fit, robust = TRUE)
  
  paths <- names(MSBB_RNA_workspace$pathologySampleSets)
  tests <- colnames(contrast.matrix)
  for(i in 1:3) {
    viral_tT <- topTable(fit,number=Inf, tests[i])
    viral_tT <- data.frame(sequence = MSBB_RNA_workspace$virus_name_2_accession_ID$VirusName[match(sapply(strsplit(rownames(viral_tT), split = "_"), function(x) paste(x[1:2], collapse = "_")), MSBB_RNA_workspace$virus_name_2_accession_ID$Accession)], name = rownames(viral_tT),viral_tT, row.names = NULL)
    featureDEPerPathologySubset[[paths[i]]] <- viral_tT
  }
  virus_level_DE_per_region[[tissues[tissue_i]]] <- featureDEPerPathologySubset
  
}


library(tidyverse)
load('results/tmp.Rdata')

orig <- virus_level_DE_per_region_origin
cont <- virus_level_DE_per_region

all(names(orig) == names(cont))

for(i in names(orig)) {
  print(all(names(orig[[i]])==names(cont[[i]]))) 
}

df_mg <- tibble()
for(i in names(orig)) {
  for(j in names(orig[[i]])) {
    df_o <- orig[[i]][[j]]
    df_c <- cont[[i]][[j]]
    df_mg <- bind_rows(
      df_mg,
      full_join(df_o, df_c, by = "name", suffix = c(".orig", ".cont")) %>% 
        mutate(tissue = i) %>% 
        mutate(test = j)
    )
    
  }
}

interests <- c("NC_001716.2_region_1_153080__ID=id0", "NC_001664.2_region_1_159322__ID=id0")

(df_mg %>% select(tissue, test, name, starts_with("logFC"), starts_with("P.Value"))) %>%
  mutate(name = ifelse(name %in% interests, substr(name, 1, 9), "")) %>% 
  ggplot(aes(x=P.Value.orig, y=P.Value.cont)) +
  geom_point() + theme_bw() + facet_grid(tissue~test) + geom_smooth(method = "lm") + theme(aspect.ratio=1) + geom_text(aes(label=name))

(df_mg %>% select(tissue, test, name, starts_with("logFC"), starts_with("P.Value"))) %>% 
  mutate(name = ifelse(name %in% interests, substr(name, 1, 9), "")) %>% 
  ggplot(aes(x=logFC.orig, y=logFC.cont)) +
  geom_point() + theme_bw() + facet_wrap(tissue~test, ncol=3, scales = "free") + geom_abline() + theme(aspect.ratio = 1) + geom_text(aes(label=name))
