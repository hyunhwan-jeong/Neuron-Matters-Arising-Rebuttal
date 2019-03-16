library(tidyverse)
load("data/MSBB_RNA_workspace.RData")
interests <- c("NC_001716.2_region_1_153080__ID=id0", "NC_001664.2_region_1_159322__ID=id0")


head(MSBB_RNA_workspace$virusLevelCounts)

head(MSBB_RNA_workspace$metadata)

head(MSBB_RNA_workspace$controlSampleSets)

is_bm22 <- MSBB_RNA_workspace$metadata$Region == "BM_22"

new_meta <- MSBB_RNA_workspace$metadata %>% 
  mutate(ad_possible = Sample_ID %in% MSBB_RNA_workspace$pathologySampleSets$AD_possible) %>% 
  mutate(ad_likely = Sample_ID %in% MSBB_RNA_workspace$pathologySampleSets$AD_likely) %>% 
  mutate(ad_definite = Sample_ID %in% MSBB_RNA_workspace$pathologySampleSets$AD_definite) %>% 
  mutate(is_control = Sample_ID %in% MSBB_RNA_workspace$controlSampleSets) %>% 
  mutate(ad_level = ad_possible + ad_likely + ad_definite)
  
all(new_meta$ad_possible >= new_meta$ad_likely)
all(new_meta$ad_likely >= new_meta$ad_definite)
all(new_meta$ad_possible != new_meta$is_control)

new_meta %>% 
#  mutate(HHV6A = MSBB_RNA_workspace$virusLevelCounts[interests[2], Sample_ID] %>% unlist) %>% 
  mutate(read_count = MSBB_RNA_workspace$virusLevelCounts[interests[1], Sample_ID] %>% unlist) %>% 
  filter(Region == "BM_22") %>% 
  filter(ad_definite | is_control) %>% 
  mutate(trait = ifelse(is_control, "Control", "AD")) %>% 
  ggplot(aes(x=read_count)) +
  geom_histogram(aes(fill=trait)) +
  facet_wrap(~trait, ncol=1, scales = "free_y")
  
new_meta %>% 
  #  mutate(HHV6A = MSBB_RNA_workspace$virusLevelCounts[interests[2], Sample_ID] %>% unlist) %>% 
  mutate(read_count = MSBB_RNA_workspace$virusLevelCounts[interests[1], Sample_ID] %>% unlist) %>% 
  filter(Region == "BM_22") %>% 
  mutate(trait = ifelse(ad_level == 0, "Control", "AD")) %>% 
  ggplot(aes(x=read_count)) +
  geom_histogram(aes(fill=trait)) +
  facet_wrap(~trait, ncol=1, scales = "free_y")

new_meta %>% 
  #  mutate(HHV6A = MSBB_RNA_workspace$virusLevelCounts[interests[2], Sample_ID] %>% unlist) %>% 
  mutate(read_count = MSBB_RNA_workspace$virusLevelCounts[interests[1], Sample_ID] %>% unlist) %>% 
  filter(Region == "BM_22") %>% 
  ggplot(aes(x=read_count)) +
  geom_histogram() +
  facet_wrap(~as.factor(ad_level), ncol=1, scales = "free_y") + 
  ggtitle("histogram of HHV7 by AD level (higher is severe), no overlap")

library(cowplot)
df_count <- new_meta %>% 
  mutate(read_count = MSBB_RNA_workspace$virusLevelCounts[interests[1], Sample_ID] %>% unlist) %>% 
  filter(Region == "BM_22") %>% 
  filter(ad_level == 0) %>% select(read_count) %>% 
  mutate(ad_level = "Control")

for(i in 1:3) {
  df_count <- 
    bind_rows(df_count,
      new_meta %>% 
      mutate(read_count = MSBB_RNA_workspace$virusLevelCounts[interests[1], Sample_ID] %>% unlist) %>% 
      filter(Region == "BM_22") %>% 
      filter(ad_level >= i) %>% select(read_count) %>% 
      mutate(ad_level = sprintf("Level :%d", i))
    )
}

ggplot(df_count, aes(x=read_count)) + geom_histogram() +
  facet_wrap(~ad_level, ncol = 1, scales="free_y")
