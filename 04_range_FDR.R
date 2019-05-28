library(tidyverse)

df <- read_csv("results/DE_results.csv") 

unique(df$method)

interests <- c("NC_001716.2_region_1_153080__ID=id0", 
               "NC_001664.2_region_1_159322__ID=id0")

df %>% 
  filter(method != "voom+limma no prefilter", method != "GLM", method != "voom+limma prefilter") %>% 
  filter(name %in% interests) %>% pull(FDR) %>% summary
