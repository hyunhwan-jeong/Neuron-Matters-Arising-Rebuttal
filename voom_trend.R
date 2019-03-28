voom_trend <- function (counts, 
                        design = NULL, 
                        lib.size = NULL, 
                        normalize.method = "none", 
                        span = 0.5, 
                        save.plot = FALSE, ...) {
  out <- list()
  if (is(counts, "DGEList")) {
    out$genes <- counts$genes
    out$targets <- counts$samples
    if (is.null(design) && diff(range(as.numeric(counts$sample$group))) > 
        0) 
      design <- model.matrix(~group, data = counts$samples)
    if (is.null(lib.size)) 
      lib.size <- counts$samples$lib.size * counts$samples$norm.factors
    counts <- counts$counts
  }
  else {
    isExpressionSet <- suppressPackageStartupMessages(is(counts, 
                                                         "ExpressionSet"))
    if (isExpressionSet) {
      if (length(Biobase::fData(counts))) 
        out$genes <- Biobase::fData(counts)
      if (length(Biobase::pData(counts))) 
        out$targets <- Biobase::pData(counts)
      counts <- Biobase::exprs(counts)
    }
    else {
      counts <- as.matrix(counts)
    }
  }
  n <- nrow(counts)
  if (n < 2L) 
    stop("Need at least two genes to fit a mean-variance trend")
  if (is.null(design)) {
    design <- matrix(1, ncol(counts), 1)
    rownames(design) <- colnames(counts)
    colnames(design) <- "GrandMean"
  }
  if (is.null(lib.size)) 
    lib.size <- colSums(counts)
  y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
  y <- normalizeBetweenArrays(y, method = normalize.method)
  fit <- lmFit(y, design, ...)
  if (is.null(fit$Amean)) 
    fit$Amean <- rowMeans(y, na.rm = TRUE)
  sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)
  sy <- sqrt(fit$sigma)
  allzero <- rowSums(counts) == 0
  if (any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  l <- lowess(sx, sy, f = span)
  data.frame(
    sx = sx,
    sy = sy,
    lx = l$x,
    ly = l$y
  )
}

run_voom_trend<- function(virMat, fullReadsPerSample, comparison_annot, design, 
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
    voom_trend(dge,designMat, lib.size = fullReadsPerSample, 
              normalize.method = "quantile")
  }
  else {
    dge <- DGEList(counts=virMat, lib.size = fullReadsPerSample) 
    voom_trend(dge[retainVir,],designMat, lib.size = fullReadsPerSample, 
              normalize.method = "quantile")
  }
}

run_trend_analysis <- function(dataset, design) {
  virMat <- dataset$virMat
  fullReadsPerSample <- dataset$fullReadsPerSample
  covariates <- dataset$covariates
  
  ret <- list()
  ret[["no prefilter"]] <- 
    run_voom_trend(
      virMat,
      fullReadsPerSample,
      covariates,
      design,
      "no_prefilter")
  
  ret[["prefilter"]] <- 
    run_voom_trend(
      virMat,
      fullReadsPerSample,
      covariates,
      design,
      "prefilter")
  ret
}