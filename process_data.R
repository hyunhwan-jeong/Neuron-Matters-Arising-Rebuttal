process_ROSMAP <- function(study_id = "ROS") {
  load("data/ROSMAP_RNA_workspace.RData")
  ROSMAP_RNA_workspace$metadata <- 
    ROSMAP_RNA_workspace$metadata[
      which(ROSMAP_RNA_workspace$metadata$Study == study_id),
      ]
  selectedCovariates <- c("Study", "AOD", "PMI","MSex", "Race", "RIN", 
                          "Batch", "CeradScore", "Education")
 
  
  ret <- list()
  CERAD <- list(
    AD_definite = c(1, 4),
    AD_likely = c(1, 2, 4),
    AD_possible = c(1, 2, 3, 4)
  )
  CERAD_CTRL <- 4
 
  for(cerad_cur in names(CERAD)) {
    sample_id <- ROSMAP_RNA_workspace$metadata$Sample_ID
    expr <- 
      ROSMAP_RNA_workspace$virusLevelCounts[, sample_id]
 
    otherCovariates <- 
      ROSMAP_RNA_workspace$metadata[
        match(colnames(expr),
              ROSMAP_RNA_workspace$metadata$Sample_ID),
        selectedCovariates
        ]
    
    selected <- otherCovariates$CeradScore %in% CERAD[[cerad_cur]]
    
    otherCovariates <- otherCovariates[selected,]
  
    otherCovariates$MSex <- as.factor(otherCovariates$MSex)
    otherCovariates$Race <- as.factor(otherCovariates$Race)
    otherCovariates$Study <- as.factor(otherCovariates$Study)
    otherCovariates$Batch <- as.factor(otherCovariates$Batch)
    otherCovariates$status <- ifelse(
      otherCovariates$CeradScore == CERAD_CTRL, 
      "Control", 
      "Case"
      )
    otherCovariates$status <- factor(otherCovariates$status, 
                                     levels = c("Control", "Case"))
    otherCovariates$CeradScore <- as.factor(otherCovariates$CeradScore)
    
    otherCovariates$AOD <- as.numeric(
      gsub(x = otherCovariates$AOD, 
           pattern = "+", 
           replacement = "", 
           fixed = TRUE))
    
    otherCovariates$PMI <- as.numeric(otherCovariates$PMI)
    otherCovariates$RIN <- as.numeric(otherCovariates$RIN)
    otherCovariates$Education <- as.numeric(otherCovariates$Education)  

    virMat <- expr[,selected] 
    fullReadsPerSample <- as.numeric(
      ROSMAP_RNA_workspace$metadata$TotalReads[
        match(colnames(virMat), ROSMAP_RNA_workspace$metadata$Sample_ID)]
    )
    
    ret[[cerad_cur]] <- list(
      virMat = virMat,
      covariates = droplevels.data.frame(otherCovariates),
      fullReadsPerSample = fullReadsPerSample
    )
  } 
  ret
}
process_MSBB <- function(tissue) {
  load("data/MSBB_RNA_workspace.RData")
  
  ret <- list()
  for(subset_i in 1:length(MSBB_RNA_workspace$pathologySampleSets)){
   
    tissue_samples <-  
      MSBB_RNA_workspace$metadata$Sample_ID[
        which(MSBB_RNA_workspace$metadata$Region == tissue)]
    
    case_selected <-
      intersect(
        MSBB_RNA_workspace$pathologySampleSets[[subset_i]], 
        tissue_samples
      )
    caseMat <- MSBB_RNA_workspace$virusLevelCounts[,case_selected]
    
    control_selected <- 
      intersect(
        MSBB_RNA_workspace$controlSampleSets, 
        tissue_samples)
    
    controlMat <- MSBB_RNA_workspace$virusLevelCounts[,control_selected]
    
   
    selected_samples <- 
      match(c(colnames(controlMat),colnames(caseMat)), 
            MSBB_RNA_workspace$metadata$Sample_ID) 
    selected_covariates <- 
      c("RIN", "AOD", "Race", "Sex", "Batch", "PMI", "NP.1")
    otherCovariates <- MSBB_RNA_workspace$metadata[selected_samples,
                                                   selected_covariates]
    otherCovariates$AOD <- gsub(x = otherCovariates$AOD, 
                                pattern = "+", 
                                replacement = "", 
                                fixed = TRUE)
    
    otherCovariates$AOD <- as.numeric(otherCovariates$AOD)
    otherCovariates$RIN <- as.numeric(otherCovariates$RIN)
    otherCovariates$PMI <- as.numeric(otherCovariates$PMI)
    
    otherCovariates$Race <- as.factor(otherCovariates$Race)
    otherCovariates$Sex <- as.factor(otherCovariates$Sex)
    otherCovariates$Batch <- as.factor(otherCovariates$Batch)

    otherCovariates$status <- factor(x = c(rep("Control", ncol(controlMat)), 
                                           rep("Case", ncol(caseMat))), 
                                     levels = c("Control","Case"))
    
    
    virMat <- cbind(controlMat,caseMat)
    
    fullReadsPerSample <- 
      as.numeric(
        MSBB_RNA_workspace$metadata$TotalReads[
          match(colnames(virMat), 
                MSBB_RNA_workspace$metadata$Sample_ID)]
        )
  
    cur_AD_level <- names(MSBB_RNA_workspace$pathologySampleSets)[[subset_i]]
    ret[[cur_AD_level]] <- list(
      virMat = virMat,
      covariates = otherCovariates,
      fullReadsPerSample = fullReadsPerSample
    )
  } 
  ret
}
process_MAYO <- function() {
  load("data/MAYO_TCX_RNA_workspace.RData")
  
  
  # Viral level summary -----------------------------------------------------
  
  exprMat <- 
    MAYO_RNA_workspace$virusLevelCounts[,
                                        MAYO_RNA_workspace$metadata$Sample_ID]
 
  selected_samples <- 
    match(
      colnames(exprMat), 
      MAYO_RNA_workspace$metadata$Sample_ID)
 
  selected_covariates <- c("RIN", "AOD", "Sex", 
                           "Source", "Flowcell", "Diagnosis") 
  otherCovariates <- 
    droplevels.data.frame(
      MAYO_RNA_workspace$metadata[selected_covariates])
  
  selected_samples <- 
    otherCovariates$Diagnosis %in% c("AD", "Control")
  
  
  otherCovariates <- otherCovariates[selected_samples,] 
  
  otherCovariates$status <-
    ifelse(otherCovariates$Diagnosis == "AD", "Case", "Control")
  otherCovariates$status <-
    factor(otherCovariates$status, levels = c("Control", "Case"))
  
  otherCovariates$AOD <- 
    as.numeric(
      gsub(x = otherCovariates$AOD, 
           pattern = "_or_above", 
           replacement = "", 
           fixed = TRUE))
  otherCovariates$RIN <- as.numeric(otherCovariates$RIN)
  
  otherCovariates$Sex <- as.factor(otherCovariates$Sex)
  otherCovariates$Diagnosis <- as.factor(
    gsub(otherCovariates$Diagnosis, 
         pattern = " ", 
         replacement = "_"))
  
  otherCovariates$Source <- as.factor(otherCovariates$Source)
  otherCovariates$Flowcell <- as.factor(otherCovariates$Flowcell)
  
  virMat <- exprMat
  
  fullReadsPerSample <- 
    as.numeric(
      MAYO_RNA_workspace$metadata$TotalReads[
        match(
          colnames(virMat), 
          MAYO_RNA_workspace$metadata$Sample_ID)
        ]
      )
  list(
    AD_definite = 
      list(
        virMat = virMat[,selected_samples],
        covariates = otherCovariates,
        fullReadsPerSample = fullReadsPerSample[selected_samples]
      )
    )  
}
