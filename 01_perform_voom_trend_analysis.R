source("process_data.R")
source("DE_utils.R")
source("voom_trend.R")
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
      run_trend_analysis(
        dataset[[cur_ad_lv]],
        ~ 0 + status + AOD + Race + RIN + Sex + Batch + PMI
      )
  }
  
}

dataset <- process_ROSMAP("ROS")
ret[["ROS_DFC"]] <- list()
for(cur_ad_lv in AD_level) {
  ret[["ROS_DFC"]][[cur_ad_lv]] <-
    run_trend_analysis(
      dataset[[cur_ad_lv]],
      ~ 0 + status + Education + AOD  + MSex + Race + PMI + RIN + Batch
    )
}

dataset <- process_ROSMAP("MAP")
ret[["MAP_DFC"]] <- list()
for(cur_ad_lv in AD_level) {
  ret[["MAP_DFC"]][[cur_ad_lv]] <-
    run_trend_analysis(
      dataset[[cur_ad_lv]],
      ~ 0 + status + Education + AOD  + MSex + Race + PMI + RIN + Batch
    )
}


ret[["MAYO_TCX"]] <-
  list(
    AD_definite = 
      run_trend_analysis(process_MAYO()$AD_definite, 
                  ~ 0 + status + AOD + Sex + Flowcell + Source + RIN)
  )

library(tidyverse)
df_trend <- tibble()
for(tissue in names(ret)) {
  for(lv in names(ret[[tissue]])) {
    cur_trend <- ret[[tissue]][[lv]]
    for(met in names(cur_trend)) {
      cur_df <- cur_trend[[met]] %>% 
        rownames_to_column("name") %>% 
        mutate(tissue = tissue,
               AD_level = lv,
               method = met) %>% 
        select(tissue, AD_level, method, everything())
      df_trend <-
        bind_rows(
          df_trend,
          cur_df
        )
    }
  }
}

write_csv(df_trend, "results/voom_trend_results.csv")

