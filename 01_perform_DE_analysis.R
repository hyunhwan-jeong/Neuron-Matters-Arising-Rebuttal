source("process_data.R")
source("DE_utils.R")

ret <- list()

AD_level <- c("AD_definite", "AD_likely", "AD_possible")
for(tissue in c("BM_22", "BM_36", "BM_10", "BM_44")) {
  tissue_id <- list(
    "BM_22" = "MSBB_BM22",
    "BM_36" = "MSBB_BM36",
    "BM_10" = "MSBB_BM10",
    "BM_44" = "MSBB_BM44"
  )
  
  dataset <- process_MSBB(tissue)
  cur <- tissue_id[[tissue]] 
  ret[[cur]] <- list()
  for(cur_ad_lv in AD_level) {
    ret[[cur]][[cur_ad_lv]] <-
      run_all_DEs(
        dataset[[cur_ad_lv]],
        ~ 0 + status + AOD + Race + RIN + Sex + Batch + PMI,
        ~ 0 + AOD + Race + RIN + Sex + Batch + PMI + status
      )
  }
  
}

dataset <- process_ROSMAP("ROS")
ret[["ROS_DFC"]] <- list()
for(cur_ad_lv in AD_level) {
  ret[["ROS_DFC"]][[cur_ad_lv]] <-
    run_all_DEs(
      dataset[[cur_ad_lv]],
      ~ 0 + status + Education + AOD  + MSex + Race + PMI + RIN + Batch, 
      ~ 0 + Education + AOD + MSex + Race + PMI + RIN + status 
    )
}

dataset <- process_ROSMAP("MAP")
ret[["MAP_DFC"]] <- list()
for(cur_ad_lv in AD_level) {
  ret[["MAP_DFC"]][[cur_ad_lv]] <-
    run_all_DEs(
      dataset[[cur_ad_lv]],
      ~ 0 + status + Education + AOD  + MSex + Race + PMI + RIN + Batch,
      ~ 0 + Education + AOD + MSex + Race + PMI + RIN + status 
    )
}


ret[["MAYO_TCX"]] <-
  list(
    AD_definite = 
      run_all_DEs(process_MAYO()$AD_definite, 
                  ~ 0 + status + AOD + Sex + Flowcell + Source + RIN,
                  ~ 0 + AOD + Sex + Flowcell + Source + RIN + status)
  )


DE_results <- ret

save(DE_results, file = "results/DE_results.RData")
