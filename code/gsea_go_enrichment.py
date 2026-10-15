library(Seurat)
library(msigdbr)
library(fgsea)
library(pheatmap)
library(Matrix)
library(Seurat)
library(fgsea)
library(tidyverse)
library(pheatmap)
library(dplyr) 
library(tibble) 
counts <- read.csv("C:/Users/86158/Desktop/论文复现1/DATA/AfterQC_count.csv", row.names = 1)
metadata <- read.csv("C:/Users/86158/Desktop/论文复现1/DATA/metadata_with_super_clusters.csv", row.names = 1)
scvi_latent <- read.csv("C:/Users/86158/Desktop/论文复现1/scvi_latent.csv", row.names = 1)

seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

scvi_latent <- as.matrix(scvi_latent)
colnames(scvi_latent) <- paste0("SCVI_", 1:ncol(scvi_latent))

clean_ids <- function(x) {
  x <- gsub("[^[:alnum:]]", ".", x)
  x <- trimws(x)
  x <- gsub("\\.+", ".", x)
  x
}

rownames(scvi_latent) <- clean_ids(rownames(scvi_latent))
colnames(seurat_obj) <- clean_ids(colnames(seurat_obj))

common_cells <- intersect(rownames(scvi_latent), colnames(seurat_obj))
seurat_obj <- subset(seurat_obj, cells = common_cells)
scvi_latent <- scvi_latent[common_cells, ]

seurat_obj[["scvi"]] <- CreateDimReducObject(
  embeddings = scvi_latent,
  key = "SCVI_",
  assay = "RNA"
)

hallmark <- fgsea::gmtPathways("C:/Users/86158/Desktop/论文复现1/DATA/h.all.v2024.1.Hs.symbols.gmt")

# 差异基因分析
Idents(seurat_obj) <- "seurat_clusters"
all_clusters <- unique(seurat_obj$seurat_clusters)
diff_genes <- list()

for (cluster in all_clusters) {
  markers <- FindMarkers(
    object = seurat_obj,
    ident.1 = cluster,
    ident.2 = NULL, 
    logfc.threshold = 0.25,
    test.use = "wilcox",
    min.pct = 0.1
  )
  
  ranked_genes <- markers %>%
    arrange(desc(avg_log2FC)) %>%
    rownames_to_column("gene") %>%
    select(gene, avg_log2FC)
  
  diff_genes[[as.character(cluster)]] <- ranked_genes
}

# GSEA分析
gsea_results <- list()

for (cluster in names(diff_genes)) {
  ranked_list <- setNames(
    diff_genes[[cluster]]$avg_log2FC,
    diff_genes[[cluster]]$gene
  )

  fgseaRes <- fgsea(
    pathways = hallmark, 
    stats = ranked_list,
    minSize = 15,
    maxSize = 500,
    nperm = 10000
  )
  
  gsea_results[[cluster]] <- fgseaRes %>%
    select(pathway, NES) %>%
    mutate(cluster = cluster)
}

combined_nes <- bind_rows(gsea_results) %>%
  pivot_wider(
    names_from = cluster,
    values_from = NES,
    values_fill = 0 
  ) %>%
  column_to_rownames("pathway")

combined_nes <- combined_nes[1:25, ]

pheatmap(
  combined_nes,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  scale = "row", 
  clustering_method = "complete",
  show_colnames = TRUE,
  fontsize_row = 8,
  main = "GSEA NES Heatmap",
  border_color = NA 
)

super_cluster_mapping <- c(
  "0" = "A", "2" = "A", "4" = "A", "5" = "A", "8" = "A",
  "3" = "B", "6" = "B", "10" = "B", "11" = "B",
  "1" = "C", "7" = "C", "9" = "C"
)

seurat_obj$super_cluster <- ifelse(
  seurat_obj$seurat_clusters %in% names(super_cluster_mapping),
  super_cluster_mapping[as.character(seurat_obj$seurat_clusters)],
  "Other"  
)

updated_metadata <- seurat_obj@meta.data
write.csv(updated_metadata, "C:/Users/86158/Desktop/论文复现1/DATA/metadata_with_super_clusters.csv", row.names = TRUE)



library(clusterProfiler)
library(org.Hs.eg.db) 
library(ggplot2)

Idents(seurat_obj) <- "super_cluster"

# 提取A类vs其他类的差异基因
markers_A <- FindMarkers(
  seurat_obj,
  ident.1 = "A",
  ident.2 = NULL, 
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

sig_genes_A <- rownames(subset(markers_A, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5))

# 基因ID转换（Symbol转Entrez ID）
entrez_ids_A <- bitr(
  sig_genes_A,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

go_A <- enrichGO(
  gene = entrez_ids_A$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",  # 生物过程
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

go_A_simplified <- clusterProfiler::simplify(go_A, cutoff = 0.7)
barplot(go_A_simplified, 
        showCategory = 10, 
        title = "GO Enrichment (Cluster A)",
        font.size = 8)



# 提取B类vs其他类的差异基因
markers_B <- FindMarkers(
  seurat_obj,
  ident.1 = "B",
  ident.2 = NULL, 
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

sig_genes_B <- rownames(subset(markers_B, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5))

entrez_ids_B <- bitr(
  sig_genes_B,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

go_B <- enrichGO(
  gene = entrez_ids_B$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",  # 生物过程
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

go_B_simplified <- clusterProfiler::simplify(go_B, cutoff = 0.7)


barplot(go_B_simplified, 
        showCategory = 10, 
        title = "GO Enrichment (Cluster B)",
        font.size = 8)


# 提取C类vs其他类的差异基因
markers_C <- FindMarkers(
  seurat_obj,
  ident.1 = "C",
  ident.2 = NULL,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# 筛选显著基因（按p值排序）
sig_genes_C <- rownames(subset(markers_C, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5))

entrez_ids_C <- bitr(
  sig_genes_C,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# GO富集（生物过程BP）
go_C <- enrichGO(
  gene = entrez_ids_C$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",  # 生物过程
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

go_C_simplified <- clusterProfiler::simplify(go_C, cutoff = 0.7)

barplot(go_C_simplified, 
        showCategory = 10, 
        title = "GO Enrichment (Cluster C)",
        font.size = 8)
