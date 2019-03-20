options(stringsAsFactors = FALSE)

library(openxlsx)
library(edgeR)
library(limma)
library(plyr)
library(pheatmap)
library(reshape2)
library(tidyverse)

# Load up workspace containing viral counts and metadata ------------------
load("data/MAYO_TCX_RNA_workspace.RData")


run_voom <- function(virMat, fullReadsPerSample, otherCovariates, MAYO_RNA_workspace, filter = F) {
  #otherCovariates$Diagnosis <- factor(as.character(otherCovariates$Diagnosis), levels = c("AD", "Control"))
  designMat = model.matrix( ~ 0 + Diagnosis + AOD + Sex + Flowcell + Source + RIN  ,data = otherCovariates)
  contrast.matrix<-makeContrasts(DiagnosisAD-DiagnosisControl, levels=designMat)
  
  sign_threshold_population <- 10
  
  virMat <- virMat[rowSums(sign(virMat)) > 0,]
  retainVir <- rowSums(virMat >= 2) >= sign_threshold_population

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
  viral_tT <- data.frame(sequence = MAYO_RNA_workspace$virus_name_2_accession_ID$VirusName[match(sapply(strsplit(rownames(viral_tT), split = "_"), function(x) paste(x[1:2], collapse = "_")), MAYO_RNA_workspace$virus_name_2_accession_ID$Accession)], name = rownames(viral_tT),viral_tT, row.names = NULL)
  viral_tT
}

run_edgeR <- function(virMat, fullReadsPerSample, otherCovariates, ROSMAP_RNA_workspace, filter = F) {
  #otherCovariates$Diagnosis <- factor(as.character(otherCovariates$Diagnosis), levels = c("AD", "Control"))
  designMat = model.matrix( ~ 0 + Diagnosis + AOD + Sex + Flowcell + Source + RIN  ,data = otherCovariates)
  contrast.matrix<-makeContrasts(DiagnosisAD-DiagnosisControl, levels=designMat)
  
  sign_threshold_population <- 10
  
  virMat <- virMat[rowSums(sign(virMat)) > 0,]
  retainVir <- rowSums(virMat >= 2) >= sign_threshold_population
  
  contrast.matrix<-makeContrasts(DiagnosisAD-DiagnosisControl, levels=designMat)
  
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
  viral_tT
}
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


#will_keep <- otherCovariates$Diagnosis == "AD" | otherCovariates$Diagnosis == "Control"
#virMat <- virMat[,will_keep]
#fullReadsPerSample <- fullReadsPerSample[will_keep]
#otherCovariates <- otherCovariates[will_keep,]

run_voom(virMat, fullReadsPerSample, otherCovariates, MAYO_RNA_workspace, filter = F) %>% arrange(P.Value)
run_voom(virMat, fullReadsPerSample, otherCovariates, MAYO_RNA_workspace, filter = T)  %>% arrange(P.Value)
run_edgeR(virMat, fullReadsPerSample, otherCovariates, MAYO_RNA_workspace, filter = F) %>% arrange(PValue)
run_edgeR(virMat, fullReadsPerSample, otherCovariates, MAYO_RNA_workspace, filter = T) %>% arrange(PValue)


deseq2_coldata <- otherCovariates
rownames(deseq2_coldata) <- colnames(virMat)
sign_threshold_population <- 10
virMat <- virMat[rowSums(sign(virMat)) > 0,]
retainVir <- rowSums(virMat >= 2) >= sign_threshold_population

# perform DESeq before the filteration
dds <- DESeq2::DESeqDataSetFromMatrix(countData = virMat, 
                              colData = deseq2_coldata, 
                              design = ~ 0 +  AOD + Sex + Flowcell + Source + RIN + Diagnosis)
dds <- DESeq2::DESeq(dds[retainVir,])
as.data.frame(DESeq2::results(dds))

viral_tT <- as.data.frame(DESeq2::results(dds))
viral_tT <- data.frame(sequence = MAYO_RNA_workspace$virus_name_2_accession_ID$VirusName[match(sapply(strsplit(rownames(viral_tT), split = "_"), function(x) paste(x[1:2], collapse = "_")), MAYO_RNA_workspace$virus_name_2_accession_ID$Accession)], name = rownames(viral_tT),viral_tT, row.names = NULL)



