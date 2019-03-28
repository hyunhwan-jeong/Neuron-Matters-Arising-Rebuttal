library(tidyverse)

interests <- c("NC_001716.2_region_1_153080__ID=id0", 
               "NC_001664.2_region_1_159322__ID=id0")

df_de <- read_csv("results/DE_results.csv") %>% 
  filter(name %in% interests) %>% 
  mutate(name = ifelse(name == interests[1], "HHV7", "HHV6A")) %>% 
  mutate(
    tissue = 
      fct_relevel(
        tissue,
        "MSBB_BM10", "MSBB_BM22", "MSBB_BM36", "MSBB_BM44", "ROS_DFC", "MAP_DFC"
      ),
    method = 
      fct_relevel(
        method, 
        unique(method)
      )
  )
  
p1 <- df_de %>%
  ggplot(aes(x=AD_level, y=-log10(FDR))) 

p2 <- df_de %>% 
  ggplot(aes(x=AD_level, y=log2FC))


draw_bars <- function(gg_obj, yintercept = 0) {
  gg_obj +
    geom_bar(aes(fill = method), 
             stat = "identity", 
             position = "dodge") +
    facet_grid(name~tissue) +
    geom_hline(yintercept = yintercept, color = "red") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90)
    )
}

library(cowplot)
plot_grid(
  draw_bars(p1, -log10(0.1)) + ggtitle("FDR"),
  draw_bars(p2) + ggtitle("Effect size"),
  ncol = 1
)
save_plot(last_plot(), 
          filename = "results/viral_stats.pdf",
          base_height = 10,
          base_width = 14)

