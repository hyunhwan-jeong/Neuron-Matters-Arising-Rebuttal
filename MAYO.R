options(stringsAsFactors = FALSE)

library(openxlsx)
library(edgeR)
library(limma)
library(plyr)
library(pheatmap)
library(reshape2)


# Load up workspace containing viral counts and metadata ------------------
load("data/MAYO_TCX_RNA_workspace.RData")


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


will_keep <- otherCovariates$Diagnosis == "AD" | otherCovariates$Diagnosis == "Control"
virMat <- virMat[,will_keep]
fullReadsPerSample <- fullReadsPerSample[will_keep]
otherCovariates <- otherCovariates[will_keep,]

otherCovariates$Diagnosis <- droplevels(otherCovariates$Diagnosis)
designMat = model.matrix( ~ 0 + Diagnosis + AOD + Sex + Flowcell + Source + RIN  ,data = otherCovariates)
contrast.matrix<-makeContrasts(DiagnosisAD-DiagnosisControl)

sign_threshold_population <- 10

virMat <- virMat[rowSums(sign(virMat)) > 0,]
retainVir <- rowSums(virMat >= 2) >= sign_threshold_population
dge <- DGEList(counts=virMat, lib.size = fullReadsPerSample)


v <- voom(dge,designMat, lib.size = fullReadsPerSample, normalize.method = "quantile")
fit <- lmFit(object = v[retainVir,], design = designMat)

contrast.matrix<-makeContrasts(
  AD_vs_Controls = DiagnosisAD-DiagnosisControl,
  levels=designMat)##

fit<- contrasts.fit(fit,contrast.matrix)
fit <- eBayes(fit, robust = TRUE)

viral_level_tT_list <- lapply(colnames(contrast.matrix), function(x) {
  
  viral_tT <- topTable(fit,number=Inf, coef = x)
  lh <- limma.one.sided.use.coef(fit, lower = TRUE, coef = x)
  rh = 1 - lh
  viral_tT <- data.frame(sequence = MAYO_RNA_workspace$virus_name_2_accession_ID$VirusName[match(sapply(strsplit(rownames(viral_tT), split = "_"), function(x) paste(x[1:2], collapse = "_")), MAYO_RNA_workspace$virus_name_2_accession_ID$Accession)], name = rownames(viral_tT),viral_tT, row.names = NULL)
  viral_tT$downregulated_pvalue <- unname(lh)[match(viral_tT$name, names(lh))]
  viral_tT$upregulated_pvalue <- unname(rh)[match(viral_tT$name, names(rh))]
  
  viral_tT
  
})

names(viral_level_tT_list) <- colnames(contrast.matrix)

viral_level_DE <- ldply(viral_level_tT_list, rbind, .id = "Comparison") 
viral_level_DE <- viral_level_DE[order(viral_level_DE$P.Value, decreasing = FALSE),]



# Viral genomic feature level summary -----------------------------------------------------

exprMat <- MAYO_RNA_workspace$virusGenomicFeatureCounts[,MAYO_RNA_workspace$metadata$Sample_ID]

otherCovariates <- MAYO_RNA_workspace$metadata[match(colnames(exprMat), MAYO_RNA_workspace$metadata$Sample_ID),c("RIN", "AOD", "Sex", "Source", "Flowcell", "Diagnosis")]

otherCovariates$AOD <- as.numeric(gsub(x = otherCovariates$AOD, pattern = "_or_above", replacement = "", fixed = TRUE))
otherCovariates$RIN <- as.numeric(otherCovariates$RIN)

otherCovariates$Sex <- as.factor(otherCovariates$Sex)
otherCovariates$Diagnosis <- as.factor(gsub(otherCovariates$Diagnosis, pattern = " ", replacement = "_"))
otherCovariates$Source <- as.factor(otherCovariates$Source)
otherCovariates$Flowcell <- as.factor(otherCovariates$Flowcell)

designMat = model.matrix( ~ 0 + Diagnosis + AOD + Sex + Flowcell + Source + RIN  ,data = otherCovariates)

# DE ---------------------------------------------------------------------

virMat <- exprMat

fullReadsPerSample <- as.numeric(MAYO_RNA_workspace$metadata$TotalReads[match(colnames(virMat), MAYO_RNA_workspace$metadata$Sample_ID)])

sign_threshold_population <- 10

virMat <- virMat[rowSums(sign(virMat)) > 0,]


retainVir <- rowSums(virMat >= 2) >= sign_threshold_population

dge <- DGEList(counts=virMat, lib.size = fullReadsPerSample)


v <- voom(dge,designMat, lib.size = fullReadsPerSample, normalize.method = "quantile")
fit <- lmFit(object = v[retainVir,], design = designMat)

contrast.matrix<-makeContrasts(
  AD_vs_Controls = DiagnosisAD-DiagnosisControl,
  PSP_vs_Controls = DiagnosisPSP-DiagnosisControl,
  PathologicalAgeing_vs_Controls = DiagnosisPathologic_Aging-DiagnosisControl,
  AD_vs_PSP = DiagnosisAD-DiagnosisPSP,
  AD_vs_PathologicalAgeing = DiagnosisAD-DiagnosisPathologic_Aging,
  PSP_vs_PathologicalAgeing = DiagnosisPSP - DiagnosisPathologic_Aging,
  levels=designMat)##

fit<- contrasts.fit(fit,contrast.matrix)
fit <- eBayes(fit, robust = TRUE)

viral_genomic_feature_level_tT_list <- lapply(colnames(contrast.matrix), function(x) {
  
  viral_tT <- topTable(fit,number=Inf, coef = x)
  lh <- limma.one.sided.use.coef(fit, lower = TRUE, coef = x)
  rh = 1 - lh
  viral_tT <- data.frame(sequence = MAYO_RNA_workspace$virus_name_2_accession_ID$VirusName[match(sapply(strsplit(rownames(viral_tT), split = "_"), function(x) paste(x[1:2], collapse = "_")), MAYO_RNA_workspace$virus_name_2_accession_ID$Accession)], name = rownames(viral_tT),viral_tT, row.names = NULL)
  viral_tT$downregulated_pvalue <- unname(lh)[match(viral_tT$name, names(lh))]
  viral_tT$upregulated_pvalue <- unname(rh)[match(viral_tT$name, names(rh))]
  
  viral_tT
  
})

names(viral_genomic_feature_level_tT_list) <- colnames(contrast.matrix)