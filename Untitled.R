options(stringsAsFactors = FALSE)

library(openxlsx)
library(edgeR)
library(limma)
library(plyr)
library(pheatmap)
library(reshape2)

# Load up workspace containing viral counts and metadata ------------------
load("data/ROSMAP_RNA_workspace.RData")

ROSMAP_RNA_workspace$metadata <- ROSMAP_RNA_workspace$metadata[which(ROSMAP_RNA_workspace$metadata$Study == "ROS"),]

# Viral level -------------------------------------------------------------

expr <- ROSMAP_RNA_workspace$virusLevelCounts[,ROSMAP_RNA_workspace$metadata$Sample_ID]

otherCovariates <- ROSMAP_RNA_workspace$metadata[match(colnames(expr), ROSMAP_RNA_workspace$metadata$Sample_ID),
                                                 c("Study", "AOD", "PMI","MSex", "Race", "RIN", "Batch", "CeradScore", "Education")]

otherCovariates$MSex <- as.factor(otherCovariates$MSex)
otherCovariates$Race <- as.factor(otherCovariates$Race)
otherCovariates$Study <- as.factor(otherCovariates$Study)
otherCovariates$Batch <- as.factor(otherCovariates$Batch)
otherCovariates$CeradScore <- as.factor(otherCovariates$CeradScore)

otherCovariates$AOD <- as.numeric(gsub(x = otherCovariates$AOD, pattern = "+", replacement = "", fixed = TRUE))
otherCovariates$PMI <- as.numeric(otherCovariates$PMI)
otherCovariates$RIN <- as.numeric(otherCovariates$RIN)
otherCovariates$Education <- as.numeric(otherCovariates$Education)


designMat = model.matrix( ~ 0 + CeradScore + Education + AOD  + MSex + Race + PMI + RIN + Batch,
                          data = droplevels.data.frame(otherCovariates))

# DE ---------------------------------------------------------------------

virMat <- expr

fullReadsPerSample <- as.numeric(ROSMAP_RNA_workspace$metadata$TotalReads[match(colnames(virMat), ROSMAP_RNA_workspace$metadata$Sample_ID)])# 
sign_threshold_population <- 10

virMat <- virMat[rowSums(sign(virMat)) > 0,]


retainVir <- rowSums(virMat >= 2) >= sign_threshold_population

# CERAD AD
# 1 Definite
# 2 Probable
# 3 Possible
# 4 No AD
# 9 Missing

contrast.matrix <- makeContrasts(cerad_1_minus_4 = CeradScore1 - CeradScore4,
                                 cerad_1_2_minus_4 = (CeradScore1+CeradScore2)/2 - CeradScore4,
                                 cerad_1_2_3_minus_4 = (CeradScore1+CeradScore2+CeradScore3)/3 - CeradScore4,
                                 levels=designMat)



dge <- DGEList(counts=virMat, lib.size = fullReadsPerSample)
v <- voom(dge,designMat, lib.size = fullReadsPerSample, normalize.method = "quantile")
fit <- lmFit(object = v[retainVir,], design = designMat)
fit<- contrasts.fit(fit,contrast.matrix)
fit <- eBayes(fit, robust = TRUE)


viral_level_tT_list <- lapply(colnames(contrast.matrix), function(x) {
  viral_tT <- topTable(fit,number=Inf, coef = x)
  viral_tT <- data.frame(sequence = ROSMAP_RNA_workspace$virus_name_2_accession_ID$VirusName[match(sapply(strsplit(rownames(viral_tT), split = "_"), function(x) paste(x[1:2], collapse = "_")), ROSMAP_RNA_workspace$virus_name_2_accession_ID$Accession)], name = rownames(viral_tT),viral_tT, row.names = NULL)
  viral_tT
})


names(viral_level_tT_list) <- colnames(contrast.matrix)

viral_level_DE <- ldply(viral_level_tT_list, rbind, .id = "Comparison") 
viral_level_DE <- viral_level_DE[order(viral_level_DE$P.Value, decreasing = FALSE),]


