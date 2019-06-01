library(tidyverse)
library(broman)
library(glue)

join_counts_and_meta <- function(HHV6A, HHV7, workspace) {
    stopifnot(all(names(HHV6A) == workspace$metadata$Sample_ID))
    stopifnot(all(names(HHV7) == workspace$metadata$Sample_ID))
    data.frame(
        workspace$metadata,
        HHV6A = HHV6A %>% unlist,
        HHV7 = HHV7 %>% unlist
    ) 
}

treat_MSBB <- function(df_MSBB, ws) {
    df_MSBB$status <-
        ifelse(
            df_MSBB$Sample_ID %in% ws$pathologySampleSets$AD_possible,
            "AD",
            ""
        )
    df_MSBB$status <-
        ifelse(
            df_MSBB$Sample_ID %in% ws$controlSampleSets,
            "Ctrl",
            df_MSBB$status
        ) %>% factor(level = c("Ctrl", "AD"))
    df_MSBB
}

treat_rosmap <- function(df) {
    df$AOD <-  as.numeric(gsub(x = df$AOD, pattern = "+", replacement = "", fixed = TRUE))
    df$MSex <- as.factor(df$MSex)
    df$Race <- as.factor(df$Race)
    df$Study <- as.factor(df$Study)
    df$Batch <- as.factor(df$Batch)
    df %>% mutate(status = ifelse(CeradScore>3, "Ctrl", "AD")) %>%
        mutate(status = factor(status, level = c("Ctrl", "AD")))
}

treat_mayo <- function(df_MAYO) {
    df_MAYO$AOD <-
        as.numeric(gsub(
            x = df_MAYO$AOD,
            pattern = "_or_above",
            replacement = "",
            fixed = TRUE
        ))
    df_MAYO$RIN <- as.numeric(df_MAYO$RIN)
    
    df_MAYO$Sex <- as.factor(df_MAYO$Sex)
    df_MAYO$status <-
        gsub(df_MAYO$Diagnosis,
             pattern = " ",
             replacement = "_") %>%
        gsub(pattern = "Control", replacement = "Ctrl")
    df_MAYO$status <-
        ifelse(df_MAYO$status == "AD" |
                   df_MAYO$status == "Ctrl",
               df_MAYO$status,
               NA)
    df_MAYO$status <-
        factor(df_MAYO$status, level = c("Ctrl", "AD"))
    df_MAYO %>% filter(!is.na(status))
}

get_cpm <- function(workspace, viral_name, offset1 = 0.1, offset2 = 0, trans_func) {
    (offset1 + workspace$virusLevelCounts[str_detect(rownames(workspace$virusLevelCounts),
                                                     viral_name),]) /
        (offset2 + workspace$metadata$TotalReads) * 10 ^ 6 -> cpm
    if(trans_func == "CPM") {
        cpm
    } else {
        log2(cpm)
    }
}

get_count <- function(workspace, viral_name) {
    workspace$virusLevelCounts[str_detect(rownames(workspace$virusLevelCounts),
                                                     viral_name),]
}


draw_dotplot <-
    function(df,
             vir_name,
             title,
             AD = "AD",
             CT = "Ctrl",
             YLIM = c(-10, 0)) {
        cnt <- df[[vir_name]]
        status <- df$status
        paste0(title, " p=", format(ks.test(cnt[status == AD], cnt[status == CT])$p.value, digits = 2)) -> title
        # dotplot(status,
        #         cnt,
        #         jiggle = "random",
        #         main = title,
        #         #ylim = YLIM,
        #         bg = ifelse(status == AD, rgb(0,0, 1, 0.1), rgb(1,0,0,0.1)), pch=20)
        # 
        x_jit <- jitter(ifelse(status == AD, 0.25, 0.75) , factor = 0.4)
        color <- ifelse(status == AD, rgb(0,0, 1, 0.1), rgb(1,0,0,0.1))
        print(cnt)
        print(x_jit)
        plot(y = cnt, x = x_jit, xaxt = "n", xlim=c(0,1), main = title, pch = 20, col = color, cex=2)
        axis(1, at=c(0.25,0.75), labels=c(AD, CT))
    }
