rm(list=ls())
source("dotplot_util.R")
load("data/MSBB_RNA_workspace.RData")
HHV6A <- get_count(MSBB_RNA_workspace, "NC_001664.2")
HHV7 <- get_count(MSBB_RNA_workspace, "NC_001716.2")
cpm <- edgeR::cpm(rbind(HHV6A, HHV7), MSBB_RNA_workspace$metadata$TotalReads, log = T, prior.count = 0.1)
HHV6A <- cpm[1,]
HHV7 <- cpm[2,]
join_counts_and_meta(HHV6A, HHV7, MSBB_RNA_workspace) %>% treat_MSBB(MSBB_RNA_workspace) %>% mutate(Study = "MSBB") -> df_MSBB


load("data/ROSMAP_RNA_workspace.RData")
HHV6A <- get_count(ROSMAP_RNA_workspace, "NC_001664.2")
HHV7 <- get_count(ROSMAP_RNA_workspace, "NC_001716.2")
cpm <- edgeR::cpm(rbind(HHV6A, HHV7), ROSMAP_RNA_workspace$metadata$TotalReads, log = T, prior.count = 0.1)
HHV6A <- cpm[1,]
HHV7 <- cpm[2,]

join_counts_and_meta(HHV6A, HHV7, ROSMAP_RNA_workspace) -> df_ROSMAP
df_ROSMAP %>% filter(Study == "ROS") %>% treat_rosmap() %>% mutate(Study = "ROS", Region = "DPFC") -> df_ROS
df_ROSMAP %>% filter(Study == "MAP") %>% treat_rosmap() %>% mutate(Study = "MAP", Region = "DPFC") -> df_MAP


load("data/MAYO_TCX_RNA_workspace.RData")
HHV6A <- get_count(MAYO_RNA_workspace, "NC_001664.2")
HHV7 <- get_count(MAYO_RNA_workspace, "NC_001716.2")
cpm <- edgeR::cpm(rbind(HHV6A, HHV7), MAYO_RNA_workspace$metadata$TotalReads, log = T, prior.count = 0.1)
HHV6A <- cpm[1,]
HHV7 <- cpm[2,]
join_counts_and_meta(HHV6A, HHV7, MAYO_RNA_workspace) %>% treat_mayo() %>% mutate(Study = "MAYO", Region = "TCX") -> df_MAYO 


bind_rows(
    df_MSBB %>% select(HHV6A, HHV7, Study, Region, status),
    df_ROS %>% select(HHV6A, HHV7, Study, Region, status), 
    df_MAP %>% select(HHV6A, HHV7, Study, Region, status),
    df_MAYO %>% select(HHV6A, HHV7, Study, Region, status)
) %>% gather("feature", "log2CPM", -Study, -Region, -status) %>% 
    #filter(Study == "MSBB") %>% 
    mutate(status = factor(status, levels = c("AD", "Ctrl"))) %>% 
    unite("dataset", c("Study", "Region"), sep = " ") %>% 
    mutate(dataset = str_remove(dataset, "_")) %>% 
    mutate(dataset = factor(dataset, 
                            levels = c("MSBB BM10",
                                       "MSBB BM22",
                                       "MSBB BM36",
                                       "MSBB BM44",
                                       "ROS DPFC",
                                       "MAP DPFC",
                                       "MAYO TCX"))) %>% 
    ggplot(aes(x=status,y=log2CPM)) +
    geom_jitter(width=0.2, alpha=0.5, size = 2) +
    #geom_boxplot(width=0.5, outlier.size = 0, coef = 0, outlier.shape = NA, fill = NA, color = "black", alpha=0.5) +
    geom_violin(width=0.5, fill=NA) +
    facet_grid(feature ~ dataset) +
    scale_y_continuous(breaks = seq(-10,1,1)) +
    #scale_color_manual(values = c("#0571b0", "#ca0020")) +
    theme_minimal() +
    theme(legend.position = "none") +
    xlab(NULL) +
    ylab("log2(CPM)") +
    theme(strip.text.x = element_text(size = 16)) +
    theme(strip.text.y = element_text(size = 16))  +
    theme(axis.title=element_text(size=16)) +
    theme(axis.text.x = element_text(size=16)) +
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                 geom = "crossbar", width = 0.5, color = "red")
ggsave("results/logCPMplot.pdf")
