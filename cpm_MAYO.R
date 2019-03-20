options(stringsAsFactors = FALSE)

library(openxlsx)
library(edgeR)
library(limma)
library(plyr)
library(pheatmap)
library(reshape2)
library(tidyverse)


# Load up workspace containing viral counts and metadata ------------------


draw_cpm_hist <- function() {
  load("data/MAYO_TCX_RNA_workspace.RData")
  interests <- c("NC_001716.2_region_1_153080__ID=id0", "NC_001664.2_region_1_159322__ID=id0")
  
  # Viral level summary -----------------------------------------------------
  
  exprMat <- MAYO_RNA_workspace$virusLevelCounts[,MAYO_RNA_workspace$metadata$Sample_ID]
  
  otherCovariates <- MAYO_RNA_workspace$metadata[match(colnames(exprMat), MAYO_RNA_workspace$metadata$Sample_ID),c("RIN", "AOD", "Sex", "Source", "Flowcell", "Diagnosis")]
  
  otherCovariates$AOD <- as.numeric(gsub(x = otherCovariates$AOD, pattern = "_or_above", replacement = "", fixed = TRUE))
  otherCovariates$RIN <- as.numeric(otherCovariates$RIN)
  
  otherCovariates$Sex <- as.factor(otherCovariates$Sex)
  otherCovariates$Diagnosis <- as.factor(gsub(otherCovariates$Diagnosis, pattern = " ", replacement = "_"))
  otherCovariates$Source <- as.factor(otherCovariates$Source)
  otherCovariates$Flowcell <- as.factor(otherCovariates$Flowcell)
  
  
  # DE ---------------------------------------------------------------------
  
  virMat <- exprMat
  fullReadsPerSample <- as.numeric(MAYO_RNA_workspace$metadata$TotalReads[match(colnames(virMat), MAYO_RNA_workspace$metadata$Sample_ID)])
  
  path_level <- c
  
  df_count <- data.frame()
  #log2CPM <- cpm(virMat, lib.size = fullReadsPerSample, log=T, prior.count = 0.1)
  
  hhv7 <- data.frame(
    cpm = unlist(virMat[interests[1],]),
    name = "HHV7",
    status = otherCovariates$Diagnosis
  )
  hhv6a <- data.frame(
    cpm = unlist(virMat[interests[2],]),
    name = "HHV6A",
    status = otherCovariates$Diagnosis
  )
  print(table(hhv6a$cpm, hhv6a$status)) 
  print(table(hhv7$cpm, hhv7$status)) 
  df_count <- rbind(df_count, hhv7, hhv6a)

  ggplot(df_count, aes(x=cpm)) +
    geom_histogram() + facet_grid(status~name) + ggtitle("MAYO") + theme_bw() +
    ylab("Count") + xlab("Raw read count")
}

draw_cpm_hist()
