library(monocle3)
library(Seurat)

cds <- new_cell_data_set(
  expression_data = GetAssayData(seurat_obj, assay = "RNA", slot = "counts"),
  cell_metadata = seurat_obj@meta.data,
  gene_metadata = data.frame(
    gene_short_name = rownames(seurat_obj),
    row.names = rownames(seurat_obj)
  )
)

reducedDims(cds)$UMAP <- seurat_obj@reductions$scvi@cell.embeddings[,1:2]

cds@preprocess_aux <- list(
  gene_loadings = matrix(0, nrow = nrow(cds), 
  prop_var_expl = rep(1, ncol(reducedDims(cds)$UMAP)
)

colData(cds)$cluster <- seurat_obj$super_cluster

cds <- learn_graph(cds, use_partition = FALSE)

root_cells <- colnames(cds)[which(
  colData(cds)$cluster == "A" & 
  colData(cds)$nCount_RNA <= quantile(
    colData(cds)$nCount_RNA[colData(cds)$cluster == "A"], 
    probs = 0.1
  )
)]
cds <- order_cells(cds, root_cells = root_cells)

plot_cells(cds, 
          color_cells_by = "cluster",
          label_cell_groups = FALSE)