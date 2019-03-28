library(edgeR)
library(DESeq2)

run_voom <- function(virMat, 
                     fullReadsPerSample, 
                     comparison_annot, design, 
                     filter = "prefilter") {
  # create a design matrix and a contrast matrix for the DE analysis
  designMat = model.matrix(design, data = comparison_annot)
  contrast.matrix<-makeContrasts(statusCase-statusControl,levels=designMat)
  
  # identify viruses will be retained in further steps (by read counts)
  sign_threshold_population <- 10
  virMat <- virMat[rowSums(sign(virMat)) > 0,]
  retainVir <- rowSums(virMat >= 2) >= sign_threshold_population
  
  if(filter != "prefilter") {
    dge <- DGEList(counts=virMat, lib.size = fullReadsPerSample) 
    v <- voom(dge,designMat, lib.size = fullReadsPerSample, 
              normalize.method = "quantile")
    fit <- lmFit(object = v[retainVir,], design = designMat)
    fit<- contrasts.fit(fit,contrast.matrix)
    fit <- eBayes(fit, robust = TRUE)
    as.data.frame(topTable(fit,number=Inf))
  }
  else {
    dge <- DGEList(counts=virMat, lib.size = fullReadsPerSample) 
    v <- voom(dge[retainVir,],designMat, lib.size = fullReadsPerSample, 
              normalize.method = "quantile")
    fit <- lmFit(object = v, design = designMat)
    fit<- contrasts.fit(fit,contrast.matrix)
    fit <- eBayes(fit, robust = TRUE)
    as.data.frame(topTable(fit,number=Inf))
  }
}

run_edgeR <- function(virMat, fullReadsPerSample, comparison_annot, design, 
                      filter = "prefilter") {
  # create a design matrix and a contrast matrix for the DE analysis
  designMat = model.matrix(design, data = comparison_annot)
  contrast.matrix<-makeContrasts(statusCase-statusControl,levels=designMat)
  # identify viruses will be retained in further steps (by read counts)
  sign_threshold_population <- 10
  virMat <- virMat[rowSums(sign(virMat)) > 0,]
  retainVir <- rowSums(virMat >= 2) >= sign_threshold_population
  
  # perform edgeR (QL F-test) and return a data.frame
  if(filter != "prefilter") {
    dge <- DGEList(counts=virMat, lib.size = fullReadsPerSample) 
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, designMat)
    fit <- glmQLFit(dge[retainVir,], designMat, robust = T)
    glf <- glmQLFTest(fit,contrast = contrast.matrix)
    as.data.frame(topTags(glf, n = Inf))
  }
  else {
    dge <- DGEList(counts=virMat[retainVir,], lib.size = fullReadsPerSample) 
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, designMat)
    fit <- glmQLFit(dge, designMat, robust = T)
    glf <- glmQLFTest(fit,contrast = contrast.matrix)
    as.data.frame(topTags(glf, n = Inf))
  }
}

run_deseq2 <- function(virMat, fullReadsPerSample, comparison_annot, design, 
                       filter = "independent") {
  deseq2_coldata <- comparison_annot
  rownames(deseq2_coldata) <- colnames(virMat)
  sign_threshold_population <- 10
  virMat <- virMat[rowSums(sign(virMat)) > 0,]
  retainVir <- rowSums(virMat >= 2) >= sign_threshold_population
  
  # perform DESeq before the filteration
  if(filter != "independent") {
    dds <- DESeqDataSetFromMatrix(countData = virMat[retainVir,], 
                                  colData = deseq2_coldata, 
                                  design = design)
    dds <- DESeq(dds)
    as.data.frame(results(dds, independentFiltering = F))
  } else {
    dds <- DESeqDataSetFromMatrix(countData = virMat,
                                  colData = deseq2_coldata, 
                                  design = design)
    
    dds <- DESeq(dds)
    as.data.frame(results(dds[retainVir,], independentFiltering = T))
  }
}

run_all_DEs <- function(dataset,
                        design_others,
                        design_deseq2) {
  
  virMat <- dataset$virMat
  fullReadsPerSample <- dataset$fullReadsPerSample
  covariates <- dataset$covariates
  
  ret <- list()
  ret[["voom+limma no prefilter"]] <- 
    run_voom(
      virMat,
      fullReadsPerSample,
      covariates,
      design_others,
      "no_prefilter")
  
  ret[["voom+limma prefilter"]] <- 
    run_voom(
      virMat,
      fullReadsPerSample,
      covariates,
      design_others,
      "prefilter")
  
  ret[["edgeR prefilter"]] <- 
    run_edgeR(
      virMat,
      fullReadsPerSample,
      covariates,
      design_others,
      "no_prefilter")
  
  ret[["edgeR no prefilter"]] <- 
    run_edgeR(
      virMat,
      fullReadsPerSample,
      covariates,
      design_others,
      "prefilter")
  
  ret[["DESeq2 indepfilter"]] <-
    run_deseq2(
      virMat, 
      fullReadsPerSample,
      covariates,
      design_deseq2,
      "independent")
  
  ret[["DESeq2 prefilter"]] <- 
    run_deseq2(
      virMat,
      fullReadsPerSample,
      covariates,
      design_deseq2,
      "prefilter")
  
  ret
}


