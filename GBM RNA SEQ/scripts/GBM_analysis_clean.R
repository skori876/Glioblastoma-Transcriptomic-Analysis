library(GEOquery)
library(limma)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(hgu133plus2.db)
library(AnnotationDbi)
gset <- getGEO("GSE4290", GSEMatrix = TRUE)
eset <- gset[[1]]

pheno <- pData(eset)
expr_matrix <- exprs(eset)
group <- pheno$`Histopathological diagnostic:ch1`

keep <- group %in% c(
  "glioblastoma, grade 4",
  "non-tumor"
)

pheno_gbm <- pheno[keep, ]
expr_gbm <- expr_matrix[, keep]
condition <- factor(
  ifelse(
    pheno_gbm$`Histopathological diagnostic:ch1` ==
      "glioblastoma, grade 4",
    "GBM",
    "Normal"
  )
)
expr_log <- log2(expr_gbm + 1)

design <- model.matrix(~0 + condition)
colnames(design) <- levels(condition)

fit <- lmFit(expr_log, design)

contrast <- makeContrasts(
  GBM - Normal,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

de_table <- topTable(
  fit2,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)
sig_genes <- subset(
  de_table,
  adj.P.Val < 0.05 & abs(logFC) > 1
)

top50 <- head(
  de_table[order(de_table$adj.P.Val), ],
  50
)
write.csv(
  de_table,
  "GSE4290_GBM_vs_Normal_DEGs.csv"
)

write.csv(
  sig_genes,
  "GSE4290_Significant_DEGs.csv"
)

write.csv(
  top50,
  "Top50_GBM_Genes.csv"
)
EnhancedVolcano(
  de_table,
  lab = rownames(de_table),
  x = "logFC",
  y = "adj.P.Val",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "GBM vs Normal"
)
top_genes <- rownames(top50)

heatmap_data <- expr_log[top_genes, ]

annotation_col <- data.frame(
  Condition = condition
)

rownames(annotation_col) <- colnames(heatmap_data)

pheatmap(
  heatmap_data,
  scale = "row",
  annotation_col = annotation_col,
  show_rownames = FALSE,
  main = "Top 50 Differentially Expressed Genes"
)
probe_ids <- rownames(sig_genes)

gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

gene_symbols <- na.omit(gene_symbols)

gene_ids <- bitr(
  unique(gene_symbols),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
ego <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH"
)

ekegg <- enrichKEGG(
  gene = gene_ids$ENTREZID,
  organism = "hsa"
)
sessionInfo()
