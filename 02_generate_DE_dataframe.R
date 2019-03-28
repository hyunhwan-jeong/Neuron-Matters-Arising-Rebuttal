load("results/DE_results.RData")
library(tidyverse)

df_de <- tibble()
for(tissue in names(DE_results)) {
  for(level in names(DE_results[[tissue]])) {
    cur_de_results <- DE_results[[tissue]][[level]]
    for(met in names(cur_de_results)) {
      cur_df <- cur_de_results[[met]] %>% 
        rownames_to_column("name")
      if(startsWith(met, "voom")) {
        cur_df <- 
          cur_df %>% 
          select(name, log2FC = logFC, FDR = adj.P.Val)
      } else if(startsWith(met, "edgeR")) {
        cur_df <- 
          cur_df %>% 
          select(name, log2FC = logFC, FDR = FDR)
      } else {
        cur_df <- 
          cur_df %>% 
          select(name, log2FC = log2FoldChange, FDR = padj)
      }
      df_de <-
        bind_rows(
          df_de,
          cur_df %>% 
            mutate(tissue = tissue) %>% 
            mutate(AD_level = level) %>% 
            mutate(method = met) %>% 
            select(tissue, AD_level, method, everything())
        )
    }
  }
}

df_de %>% write_csv("results/DE_results.csv")
