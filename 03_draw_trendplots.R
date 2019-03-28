vir_id <- 
  c("NC_001664.2_region_1_159322__ID=id0", 
    "NC_001716.2_region_1_153080__ID=id0")
df_trend <- read_csv("results/voom_trend_results.csv") %>% 
  filter(
    #tissue %in% c("MSBB_BM10", "MSBB_BM22"),
    AD_level == "AD_possible" |
      tissue == "MAYO_TCX"
  ) %>% 
  mutate(
    tissue = fct_relevel(
      tissue,
      "MSBB_BM10", "MSBB_BM22", "MSBB_BM36", "MSBB_BM44", 
      "ROS_DFC", "MAP_DFC", "MAYO_TCX"
    )
  )

library(cowplot)  
df_trend %>% 
  ggplot(aes(x=sx, y=sy)) +
  geom_point(alpha=0.5, size=0.3) +
  geom_line(aes(x=lx, y=ly), alpha=0.5, color="red") +
  geom_point(
    data = df_trend %>% filter(name %in% vir_id),
    aes(x=sx, y=sy), 
    size=0.8,
    color = "red"
  )  + 
  facet_grid(method~tissue) +
  theme_bw() + xlab("log2(count+0.5)") + ylab("sqrt(SD)")

save_plot(last_plot(), 
          filename = "results/viral_trends.pdf",
          base_height = 8,
          base_width = 16)
# x <- read_csv("results/voom_trend_results.csv") %>% 
#   filter(
#     tissue == "ROS_DFC",
#     AD_level == "AD_possible",
#     method == "prefilter"
#   )

# plot(x=x$sx,
#      y=x$sy)
# 
# lines(lowess(x$sx, x$sy, f=0.5))
