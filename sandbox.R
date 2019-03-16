# source: https://github.com/csoneson/conquer_comparison
library(MAST)
library(edgeR)

load("data.Rdata")
sign_threshold_population <- 10
mat <- virMat[rowSums(sign(virMat)) > 0,]
retainVir <- rowSums(mat >= 2) >= sign_threshold_population
dge <- DGEList(counts = mat[retainVir,], lib.size = fullReadsPerSample)
dge <- calcNormFactors(dge)
cpms <- cpm(dge)
cdr <- scale(colMeans(mat[retainVir,]> 0))
comparison_annot <- data.frame(comparison_annot, wellKey = colnames(cpms), cdr = cdr)
sca <- FromMatrix(exprsArray = log2(cpms+1),
                  cData = comparison_annot)

zlmdata <- zlm(~ 0 +  AOD + Race + RIN + Sex + Batch + PMI + status, sca)
mast <- lrTest(zlmdata, "status")

df <- data.frame(pval = mast[, "hurdle", "Pr(>Chisq)"], row.names = names(mast[, "hurdle", "Pr(>Chisq)"])) 
df$fdr <- p.adjust(df$pval, method = "fdr")

