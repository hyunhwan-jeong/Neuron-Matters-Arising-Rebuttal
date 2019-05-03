source("process_data.R")
library(tidyverse)

get_count_df <- function(dataset, trait) {
  require(edgeR)
  viral_cpm <- cpm(dataset$virMat, 
                   lib.size = dataset$fullReadsPerSample, 
                   log = T, 
                   prior.count = 0.1)
  
  vir_name <- 
    c("HHV6A", "HHV7")

  vir_id <- 
    c("NC_001664.2_region_1_159322__ID=id0", 
    "NC_001716.2_region_1_153080__ID=id0")
  df_expr <- tibble()
  for(i in 1:2) {
    df_expr <- 
      bind_rows(
        df_expr,
        tibble(
          name = vir_name[i],
          CPM = unlist(viral_cpm[vir_id[i],]),
          rawcount = unlist(dataset$virMat[vir_id[i],]),
          AD_level = as.factor(dataset$covariates[,trait]))
      )
  }
  df_expr
}

df_expr <- tibble()
for(tissue in c("BM_22", "BM_36", "BM_10", "BM_44")) {
  tissue_id <- list(
    "BM_22" = "MSBB_BM22",
    "BM_36" = "MSBB_BM36",
    "BM_10" = "MSBB_BM10",
    "BM_44" = "MSBB_BM44"
  )
  
  dataset <- process_MSBB(tissue)
  cur <- tissue_id[[tissue]] 
  df_expr <-
    bind_rows(
      df_expr,
      get_count_df(dataset$AD_possible, "NP.1") %>% 
        mutate(tissue = cur) %>% 
        mutate(AD_level =
                 fct_recode(AD_level, 
                            Control = "1",
                            AD_definite = "2",
                            AD_likely = "3",
                            AD_possible = "4")) %>% 
        select(tissue, everything()))
}

for(tissue in c("ROS", "MAP")) {
  tissue_id <- list(
    "ROS" = "ROS_DFC",
    "MAP" = "MAP_DFC"
  )
  
  dataset <- process_ROSMAP(tissue)
  cur <- tissue_id[[tissue]] 
  df_expr <-
    bind_rows(
      df_expr,
      get_count_df(dataset$AD_possible, "CeradScore") %>% 
        mutate(tissue = cur) %>% 
        mutate(AD_level =
                 fct_recode(AD_level, 
                            Control = "4",
                            AD_definite = "1",
                            AD_likely = "2",
                            AD_possible = "3")) %>%         
        select(tissue, everything()))
}

dataset <- process_MAYO()
df_expr <-
  bind_rows(
    df_expr,
    get_count_df(dataset$AD_definite, "Diagnosis") %>% 
      mutate(tissue = "MAYO_TCX") %>% 
      mutate(AD_level = str_replace(AD_level, "AD", "AD_definite")) %>% 
      select(tissue, everything()))

df_expr <-
  df_expr %>% 
  mutate(
    tissue = 
      fct_relevel(
        tissue,
        "MSBB_BM10", "MSBB_BM22", "MSBB_BM36", "MSBB_BM44", 
        "ROS_DFC", "MAP_DFC", "MAYO_TCX"
      )
  )

library(cowplot)
p1 <- df_expr %>% 
  filter(name == "HHV6A") %>% 
  ggplot(aes(rawcount)) +
  geom_histogram(aes(y=..density..*..width..)) +
  facet_grid(AD_level ~ tissue, scales = "free_x") +
  ggtitle("HHV6A rawcount") +
  ylab("Frequency") +
  theme_bw()
  
p2 <- df_expr %>% 
  filter(name == "HHV7") %>% 
  ggplot(aes(rawcount)) +
  geom_histogram(aes(y=..density..*..width..)) +
  facet_grid(AD_level ~ tissue, scales = "free_x") +
  ggtitle("HHV7 rawcount") +
  ylab("Frequency") +  
  theme_bw()

p3 <- df_expr %>% 
  filter(name == "HHV6A") %>% 
  ggplot(aes(CPM)) +
  geom_histogram(aes(y=..density..*..width..)) +
  facet_grid(AD_level ~ tissue, scales = "free_x") +
  ggtitle("HHV6A") +
  ylab("Frequency") + xlab("log2CPM") +
  theme_bw()

p4 <- df_expr %>% 
  filter(name == "HHV7") %>% 
  ggplot(aes(CPM)) +
  geom_histogram(aes(y=..density..*..width..)) +
  facet_grid(AD_level ~ tissue, scales = "free_x") +
  ggtitle("HHV7") +  
  ylab("Frequency") + xlab("log2CPM") +
  theme_bw()

plot_grid(p3, p4, ncol = 2)
save_plot(last_plot(), 
          filename = "results/viral_expression.pdf",
          base_height = 6,
          base_width = 16)
