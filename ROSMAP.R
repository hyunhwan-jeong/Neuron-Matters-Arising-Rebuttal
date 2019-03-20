options(stringsAsFactors = FALSE)

library(openxlsx)
library(edgeR)
library(limma)
library(plyr)
library(pheatmap)
library(reshape2)
library(DESeq2)
# Load up workspace containing viral counts and metadata ------------------

run_voom <- function(virMat, fullReadsPerSample, otherCovariates, ROSMAP_RNA_workspace, filter = F) {
  otherCovariates$status <- ifelse(otherCovariates$CeradScore == 4, "Control", "Case")
  
  designMat = model.matrix( ~ 0 + status + Education + AOD  + MSex + Race + PMI + RIN + Batch,
                            data = droplevels.data.frame(otherCovariates))
  
  sign_threshold_population <- 10
  virMat <- virMat[rowSums(sign(virMat)) > 0,]
  retainVir <- rowSums(virMat >= 2) >= sign_threshold_population
  print(sum(retainVir))
  
  contrast.matrix <- makeContrasts(statusCase - statusControl,
                                   levels=designMat)

  if(!filter) {
    dge <- DGEList(counts=virMat, lib.size = fullReadsPerSample) 
    v <- voom(dge,designMat, lib.size = fullReadsPerSample, normalize.method = "quantile")
    fit <- lmFit(object = v[retainVir,], design = designMat)
    fit<- contrasts.fit(fit,contrast.matrix)
    fit <- eBayes(fit, robust = TRUE)
  }
  else {
    dge <- DGEList(counts=virMat[retainVir,], lib.size = fullReadsPerSample) 
    v <- voom(dge, designMat, lib.size = fullReadsPerSample, normalize.method = "quantile")
    fit <- lmFit(object = v, design = designMat)
    fit<- contrasts.fit(fit,contrast.matrix)
    fit <- eBayes(fit, robust = TRUE)
  }
  
  viral_tT <- topTable(fit,number=Inf)
  viral_tT <- data.frame(sequence = ROSMAP_RNA_workspace$virus_name_2_accession_ID$VirusName[match(sapply(strsplit(rownames(viral_tT), split = "_"), function(x) paste(x[1:2], collapse = "_")), ROSMAP_RNA_workspace$virus_name_2_accession_ID$Accession)], name = rownames(viral_tT),viral_tT, row.names = NULL)
}

run_edgeR <- function(virMat, fullReadsPerSample, otherCovariates, ROSMAP_RNA_workspace, filter = F) {
  otherCovariates$status <- ifelse(otherCovariates$CeradScore == 4, "Control", "Case")
  
  designMat = model.matrix( ~ 0 + status + Education + AOD  + MSex + Race + PMI + RIN + Batch,
                            data = droplevels.data.frame(otherCovariates))
  
  sign_threshold_population <- 10
  virMat <- virMat[rowSums(sign(virMat)) > 0,]
  retainVir <- rowSums(virMat >= 2) >= sign_threshold_population
  print(sum(retainVir))
  
  contrast.matrix <- makeContrasts(statusCase - statusControl,
                                   levels=designMat)
  
  if(!filter){
    dge <- DGEList(counts=virMat, lib.size = fullReadsPerSample)
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, designMat)
    fit <- glmQLFit(dge[retainVir,], designMat)
    glf <- glmQLFTest(fit,contrast = contrast.matrix)
  }
  else {
    dge <- DGEList(counts=virMat[retainVir,], lib.size = fullReadsPerSample) 
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, designMat)
    fit <- glmQLFit(dge, designMat)
    glf <- glmQLFTest(fit,contrast = contrast.matrix)
    as.data.frame(topTags(glf, n = Inf))
  }
  
  viral_tT <- as.data.frame(topTags(glf, n = Inf))
  viral_tT <- data.frame(sequence = ROSMAP_RNA_workspace$virus_name_2_accession_ID$VirusName[match(sapply(strsplit(rownames(viral_tT), split = "_"), function(x) paste(x[1:2], collapse = "_")), ROSMAP_RNA_workspace$virus_name_2_accession_ID$Accession)], name = rownames(viral_tT),viral_tT, row.names = NULL)
}

run_deseq2 <- function(virMat, fullReadsPerSample, otherCovariates, ROSMAP_RNA_workspace, filter = F) {
  otherCovariates$status <- ifelse(otherCovariates$CeradScore == 4, "ctrl", "trt")
  deseq2_coldata <- otherCovariates
  rownames(deseq2_coldata) <- colnames(virMat)
  sign_threshold_population <- 10
  virMat <- virMat[rowSums(sign(virMat)) > 0,]
  retainVir <- rowSums(virMat >= 2) >= sign_threshold_population
  
  # perform DESeq before the filteration
  if(!filter) {
    dds <- DESeqDataSetFromMatrix(countData = virMat, 
                                  colData = deseq2_coldata, 
                                  design = ~ 0 + Education + AOD  + MSex + Race + PMI + RIN + Batch + status)
    dds <- DESeq(dds[retainVir,])
    as.data.frame(results(dds))
  } else {
    dds <- DESeqDataSetFromMatrix(countData = virMat[retainVir,], 
                                  colData = deseq2_coldata, 
                                  design =  ~ 0 + Education + AOD  + MSex + Race + PMI + RIN + Batch + status)
    
    print(summary(dds))
    dds <- DESeq(dds)
  }
  viral_tT <- as.data.frame(results(dds))
  viral_tT <- data.frame(sequence = ROSMAP_RNA_workspace$virus_name_2_accession_ID$VirusName[match(sapply(strsplit(rownames(viral_tT), split = "_"), function(x) paste(x[1:2], collapse = "_")), ROSMAP_RNA_workspace$virus_name_2_accession_ID$Accession)], name = rownames(viral_tT),viral_tT, row.names = NULL)
}


run_rosmap <- function(study_id = "ROS") {
  load("data/ROSMAP_RNA_workspace.RData")
  ROSMAP_RNA_workspace$metadata <- ROSMAP_RNA_workspace$metadata[which(ROSMAP_RNA_workspace$metadata$Study == study_id),]
  
  # Viral level -------------------------------------------------------------
  
  path_level <- c("AD_Definite", "AD_Likely", "AD_Possible")
  
  df_ret <- list()
  for(cerad_i in 1:3) {
    expr <- ROSMAP_RNA_workspace$virusLevelCounts[,ROSMAP_RNA_workspace$metadata$Sample_ID]
    
    otherCovariates <- ROSMAP_RNA_workspace$metadata[match(colnames(expr), ROSMAP_RNA_workspace$metadata$Sample_ID),
                                                     c("Study", "AOD", "PMI","MSex", "Race", "RIN", "Batch", "CeradScore", "Education")]
    
    otherCovariates$MSex <- as.factor(otherCovariates$MSex)
    otherCovariates$Race <- as.factor(otherCovariates$Race)
    otherCovariates$Study <- as.factor(otherCovariates$Study)
    otherCovariates$Batch <- as.factor(otherCovariates$Batch)
    
    
    otherCovariates$AOD <- as.numeric(gsub(x = otherCovariates$AOD, pattern = "+", replacement = "", fixed = TRUE))
    otherCovariates$PMI <- as.numeric(otherCovariates$PMI)
    otherCovariates$RIN <- as.numeric(otherCovariates$RIN)
    otherCovariates$Education <- as.numeric(otherCovariates$Education)
    
    virMat <- expr
    fullReadsPerSample <- as.numeric(ROSMAP_RNA_workspace$metadata$TotalReads[match(colnames(virMat), ROSMAP_RNA_workspace$metadata$Sample_ID)])# 
    
    interests <- otherCovariates$CeradScore == 4 | otherCovariates$CeradScore <= cerad_i
    virMat <- virMat[,interests]
    fullReadsPerSample <- fullReadsPerSample[interests]
    otherCovariates <- otherCovariates[interests,]
    
    df_ret[[path_level[cerad_i]]] <- list()
    df_ret[[path_level[cerad_i]]][["VoomLimmaNoFilter"]] <- run_voom(virMat, fullReadsPerSample, otherCovariates, ROSMAP_RNA_workspace, filter = F)
    df_ret[[path_level[cerad_i]]][["VoomLimmaWithFilter"]] <- run_voom(virMat, fullReadsPerSample, otherCovariates, ROSMAP_RNA_workspace, filter = T)
    df_ret[[path_level[cerad_i]]][["edgeRNoFilter"]] <- run_edgeR(virMat, fullReadsPerSample, otherCovariates, ROSMAP_RNA_workspace, filter = F)
    df_ret[[path_level[cerad_i]]][["edgeRWithFilter"]] <- run_edgeR(virMat, fullReadsPerSample, otherCovariates, ROSMAP_RNA_workspace, filter = F)
    df_ret[[path_level[cerad_i]]][["DESeq2"]] <- run_deseq2(virMat, fullReadsPerSample, otherCovariates, ROSMAP_RNA_workspace)
    
  }
  df_ret
  
}

ret <- list()
ret[["ROS"]] <- run_rosmap("ROS")
ret[["MAP"]] <- run_rosmap("MAP")
library(tidyverse)

interests <- c("NC_001716.2_region_1_153080__ID=id0", "NC_001664.2_region_1_159322__ID=id0")

df_all <- data.frame()
for(study_id in c("ROS", "MAP")) {
  df_cur <- ret[[study_id]]
  for(path_level in names(df_cur)) {
    for(met in names(df_cur[[path_level]])) {
      df_ft <- df_cur[[path_level]][[met]] %>% filter(name %in% interests) %>% select(sequence, matches("logFC|log2FoldChange"), matches("padj|FDR|adj.P.Val")) 
      df_all <- rbind(df_all, 
                      df_ft %>% select(name = 1, logFC = 2, FDR = 3) %>% 
                        mutate(study = study_id) %>% mutate(method = met) %>% 
                        mutate(AD_level = path_level))
    }
  }
}

ggplot(df_all, aes(x=AD_level, y=-log10(FDR))) +
  geom_bar(aes(fill=method), stat = "identity", position = "dodge") +
  geom_hline(yintercept = -log10(0.1)) +
  facet_grid(study~name)

ggplot(df_all, aes(x=AD_level, y=logFC)) +
  geom_bar(aes(fill=method), stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0) +
  facet_grid(study~name)


draw_cpm_hist <- function(study_id) {
  load("data/ROSMAP_RNA_workspace.RData")
  ROSMAP_RNA_workspace$metadata <- ROSMAP_RNA_workspace$metadata[which(ROSMAP_RNA_workspace$metadata$Study == study_id),]
  expr <- ROSMAP_RNA_workspace$virusLevelCounts[,ROSMAP_RNA_workspace$metadata$Sample_ID]
  
  path_level <- c("AD_Definite", "AD_Likely", "AD_Possible")
  
  df_count <- data.frame()
  for(cerad_i in 1:3) {
    
    otherCovariates <- ROSMAP_RNA_workspace$metadata[match(colnames(expr), ROSMAP_RNA_workspace$metadata$Sample_ID),
                                                     c("Study", "AOD", "PMI","MSex", "Race", "RIN", "Batch", "CeradScore", "Education")]
    
    virMat <- expr
    fullReadsPerSample <- as.numeric(ROSMAP_RNA_workspace$metadata$TotalReads[match(colnames(virMat), ROSMAP_RNA_workspace$metadata$Sample_ID)])# 
    will_keep <- otherCovariates$CeradScore == 4 | otherCovariates$CeradScore <= cerad_i
    
    virMat <- virMat[,will_keep]
    fullReadsPerSample <- fullReadsPerSample[will_keep]
    otherCovariates <- otherCovariates[will_keep,]  
    
    log2CPM <- cpm(virMat, lib.size = fullReadsPerSample, log=T, prior.count = 0.1)
    
    hhv7 <- data.frame(
      cpm = unlist(log2CPM[interests[1],]),
      name = "HHV7",
      status = ifelse(otherCovariates$CeradScore == 4, "Control", path_level[cerad_i])
    )
    hhv6a <- data.frame(
      cpm = unlist(log2CPM[interests[2],]),
      name = "HHV6A",
      status = ifelse(otherCovariates$CeradScore == 4, "Control", path_level[cerad_i])
    )
    
    if(cerad_i != 1) {
      hhv7 <- hhv7 %>% filter(status != "Control")
      hhv6a <- hhv6a %>% filter(status != "Control")
    }
    df_count <- rbind(df_count, hhv7, hhv6a)
  }
  
  ggplot(df_count, aes(x=cpm)) +
    geom_histogram(aes(y=stat(width*density))) + facet_grid(status~name) + ggtitle(study_id) + theme_bw() +
    ylab("Frequency (%)") + xlab("log2CPM")
}

library(cowplot)
plot_grid(
  draw_cpm_hist("ROS"),
  draw_cpm_hist("MAP")  
)

