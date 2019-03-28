library(openxlsx)
library(edgeR)
library(limma)
library(plyr)
library(pheatmap)
library(reshape2)
library(tidyverse)
load("data/MAYO_TCX_RNA_workspace.RData")

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
designMat = model.matrix( ~ 0 + Diagnosis + AOD + Sex + Flowcell + Source + RIN  ,data = otherCovariates)
sign_threshold_population <- 10

virMat <- virMat[rowSums(sign(virMat)) > 0,]
retainVir <- rowSums(virMat >= 2) >= sign_threshold_population

dge <- DGEList(counts=virMat, lib.size = fullReadsPerSample) 
v <- voom(dge,designMat, lib.size = fullReadsPerSample, normalize.method = "quantile")
cont<-makeContrasts(DiagnosisAD-DiagnosisControl, levels=designMat)
sign_threshold_population <- 10

fit <- lmFit(object = v, design = designMat)
fit <- contrasts.fit(fit,cont)
fit <- eBayes(fit, robust = TRUE)


virMat <- exprMat
fullReadsPerSample <- as.numeric(MAYO_RNA_workspace$metadata$TotalReads[match(colnames(virMat), MAYO_RNA_workspace$metadata$Sample_ID)])
will_keep <- otherCovariates$Diagnosis == "AD" | otherCovariates$Diagnosis == "Control"
virMat <- virMat[,will_keep]
fullReadsPerSample <- fullReadsPerSample[will_keep]
otherCovariates <- otherCovariates[will_keep,]
otherCovariates$Diagnosis <- as.factor(droplevels(otherCovariates$Diagnosis))
designMat = model.matrix( ~ 0 + Diagnosis + AOD + Sex + Flowcell + Source + RIN, data = otherCovariates)
sign_threshold_population <- 10

virMat <- virMat[rowSums(sign(virMat)) > 0,]
retainVir <- rowSums(virMat >= 2) >= sign_threshold_population

dge2 <- DGEList(counts=virMat, lib.size = fullReadsPerSample) 
v2 <- voom(dge2, designMat, lib.size = fullReadsPerSample, normalize.method = "quantile")
cont2<-makeContrasts(DiagnosisAD-DiagnosisControl, levels=designMat)

fit2 <- lmFit(object = v2, design = designMat)
fit2 <- contrasts.fit(fit2 , cont2)
fit2 <- eBayes(fit2, robust = TRUE)

df_v <- v$E %>% as.data.frame %>% rownames_to_column("id") %>%  gather("sample", "voom", -id)
df_v2 <- v2$E %>% as.data.frame %>% rownames_to_column("id") %>%  gather("sample", "voom", -id)

df_v %>% full_join(df_v2, by = c("id", "sample")) %>% 
  ggplot(aes(x=voom.x, y=voom.y)) +
  geom_point() +
  geom_abline()


df1 <- topTable(fit, number = Inf) %>% rownames_to_column("id")
df2 <- topTable(fit2, number = Inf) %>% rownames_to_column("id")

left_join(df1, df2, by = "id") %>% 
  ggplot(aes(x=logFC.x, y=logFC.y)) +
  geom_point()

interests <- c("NC_001716.2_region_1_153080__ID=id0", "NC_001664.2_region_1_159322__ID=id0")


left_join(df1, df2, by = "id") %>% 
  mutate(id = ifelse(id %in% interests, id, "")) %>% 
  ggplot(aes(x=logFC.x, y=logFC.y)) +
  geom_point() +
  geom_abline() +
  geom_text(aes(label=id)) 

left_join(df1, df2, by = "id") %>% 
  mutate(id = ifelse(id %in% interests, id, "")) %>% 
  ggplot(aes(x=-log10(P.Value.x), y=-log10(P.Value.y))) +
  geom_point() +
  geom_abline() +
  geom_text(aes(label=id)) 

left_join(df1, df2, by = "id") %>% 
  mutate(id = ifelse(id %in% interests, id, "")) %>% 
  ggplot(aes(x=-log10(adj.P.Val.x), y=-log10(adj.P.Val.y))) +
  geom_point() +
  geom_abline() +
  geom_hline(yintercept = -log10(0.1)) +
  geom_vline(xintercept = -log10(0.1)) +
  geom_text(aes(label=id)) 
