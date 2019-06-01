source("MSBB_Viral_differential_RNA_abundance_analysis.R")

pvalMat_label_df <- format(Pval_expansive_matrix, scientific = TRUE, digits = 2)
#pvalMat_label_df[data.matrix(FDR_expansive_matrix) > 0.1] <- ""
pvalMat_label_df[is.na(data.matrix(FDR_expansive_matrix))] <- ""

for(i in 1:nrow(pvalMat_label_df)) {
  for(j in 1:ncol(pvalMat_label_df)) {
    if( is.na(Pval_expansive_matrix[i,j]) ) next()
    if( Pval_expansive_matrix[i,j] >= 0.01 ) {
      pvalMat_label_df[i,j] <- sprintf("%.2f", Pval_expansive_matrix[i,j])
    }
  }
}
logp <- -log10(Pval_expansive_matrix) * sign(DE_expansive_matrix)
logp[is.na(logp)] <- 0
dataValues <- sort(unique(c(unlist(logp))),decreasing=FALSE)
outerLim <- ceiling(max(abs(dataValues)))
bkPrezero <- seq((-1 * outerLim)-0.01,to=-0.01,0.01)
bkPostzero <- seq(0.01,(outerLim+0.01),0.01)

bks <- c(bkPrezero,bkPostzero)
colPrezero <- colorRampPalette(c("darkblue", "aliceblue"))(length(bkPrezero)-1)
colPostzero <- colorRampPalette(c("salmon", "darkred"))(length(bkPostzero)-1)
colSeq <- c(colPrezero, "white",colPostzero)

orderCols <- unlist(lapply(c("BM_10", "BM_22", "BM_36", "BM_44"), function(x) paste(x, c("AD_definite", "AD_likely","AD_possible"), sep = "_")))

colGaps <- c(3,6,9)

#Collectively Significant
retainVir <- rowSums(FDR_expansive_matrix <= 0.1, na.rm = TRUE) > 0


pheatmap(gaps_col = colGaps, mat = logp[retainVir,orderCols], cellwidth = 45, cellheight = 20, 
               filename =  "results/heatmap_asym.pdf",
               color=colSeq,breaks=bks,
               show_rownames = TRUE, display_numbers = pvalMat_label_df[retainVir,orderCols],
               border_color="white",scale="none", cluster_rows = T, cluster_cols = FALSE,
               number_color = "white", fontsize_number = 10)

colPostzero <- colorRampPalette(c("#fff8f0", "darkred"))(length(bkPostzero)-1)
colSeq <- c(colPrezero, "white",colPostzero)

pheatmap(gaps_col = colGaps, mat = logp[retainVir,orderCols], cellwidth = 45, cellheight = 20, 
         filename =  "results/heatmap_sym.pdf",
         color=colSeq,breaks=bks,
         show_rownames = TRUE, display_numbers = pvalMat_label_df[retainVir,orderCols],
         border_color="white",scale="none", cluster_rows = T, cluster_cols = FALSE,
         number_color = "white", fontsize_number = 10)

