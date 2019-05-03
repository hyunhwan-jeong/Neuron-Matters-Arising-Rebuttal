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
  interests <- c("NC_001716.2_region_1_153080__ID=id0", 
                 "NC_001664.2_region_1_159322__ID=id0")
  
  retainVir <- (rowSums(virMat >= 2) >= sign_threshold_population) | (rownames(virMat) %in% interests)
  
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
  interests <- c("NC_001716.2_region_1_153080__ID=id0", 
                 "NC_001664.2_region_1_159322__ID=id0")
  
  retainVir <- (rowSums(virMat >= 2) >= sign_threshold_population) | (rownames(virMat) %in% interests)
  
  
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
  interests <- c("NC_001716.2_region_1_153080__ID=id0", 
                 "NC_001664.2_region_1_159322__ID=id0")
  
  retainVir <- (rowSums(virMat >= 2) >= sign_threshold_population) | (rownames(virMat) %in% interests)
  
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

run_glm <- function(virMat, 
                    fullReadsPerSample, 
                    comparison_annot, 
                    design) {
  
  
  glm_func <- function(mat, cov, f = NULL) {
    require(magrittr)
    get_coeff <- function(dat, f, i) {
      glm(f, data = dat, family = binomial(), na.action = na.exclude) %>% 
        summary %>% coefficients -> coef
      
      logFC <- log2(mean(dat$cpm[dat$status=="Case"]) / mean(dat$cpm[dat$status=="Control"]))
      #logFC <- mean(log2(dat$cpm[dat$status=="Case"]+0.01)) - mean(log2(dat$cpm[dat$status=="Control"]+0.01))
      
      list(P.Value = coef[i, "Pr(>|z|)"],
           z = coef[i, "z value"],
           logFC = logFC)
    }
    
    if(is.null(f)) {
      f <- status ~ 0 + cpm + AOD + Race + RIN + Sex + Batch + PMI
    }
    df_all_glm <- data.frame()
    for(i in 1:nrow(mat)) {
      cpm <- unlist(mat[i,]) 
      df_all <- data.frame(cov, cpm = cpm)
      
      df_all_glm <- rbind(df_all_glm, get_coeff(df_all, f, "cpm"))
    }
    df_all_glm$adj.P.Val <- p.adjust(df_all_glm$P.Val, method = "fdr")
    rownames(df_all_glm) <- rownames(mat)
    df_all_glm
  }
  
  sign_threshold_population <- 10
  virMat <- virMat[rowSums(sign(virMat)) > 0,]
  interests <- c("NC_001716.2_region_1_153080__ID=id0", 
                 "NC_001664.2_region_1_159322__ID=id0")
  
  retainVir <- (rowSums(virMat >= 2) >= sign_threshold_population) | (rownames(virMat) %in% interests)
  
  
  cpm <- edgeR::cpm(virMat, lib.size = fullReadsPerSample)
  
  glm_func(cpm[retainVir,], comparison_annot, design)
  
}



run_all_DEs <- function(dataset,
                        design_others,
                        design_deseq2,
                        design_glm) {
  
  virMat <- dataset$virMat
  fullReadsPerSample <- dataset$fullReadsPerSample
  covariates <- dataset$covariates
  
  ret <- list()
  ret[["voom+limma"]] <- 
    run_voom(
      virMat,
      fullReadsPerSample,
      covariates,
      design_others,
      "no_prefilter")
  
  # ret[["voom+limma prefilter"]] <- 
  #   run_voom(
  #     virMat,
  #     fullReadsPerSample,
  #     covariates,
  #     design_others,
  #     "prefilter")
  
  ret[["edgeR"]] <- 
    run_edgeR(
      virMat,
      fullReadsPerSample,
      covariates,
      design_others,
      "no_prefilter")
  
  # ret[["edgeR no prefilter"]] <- 
  #   run_edgeR(
  #     virMat,
  #     fullReadsPerSample,
  #     covariates,
  #     design_others,
  #     "prefilter")
  
  ret[["DESeq2"]] <-
    run_deseq2(
      virMat, 
      fullReadsPerSample,
      covariates,
      design_deseq2,
      "independent")
  
  # ret[["DESeq2 prefilter"]] <- 
  #   run_deseq2(
  #     virMat,
  #     fullReadsPerSample,
  #     covariates,
  #     design_deseq2,
  #     "prefilter")
  
  ret[["GLM"]] <-
    run_glm(
      virMat,
      fullReadsPerSample,
      covariates,
      design_glm
    )
  
  ret
}


