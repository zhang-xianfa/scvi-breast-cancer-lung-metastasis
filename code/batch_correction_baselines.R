library(Seurat)
library(sva)
library(harmony)
library(batchelor)
library(kBET)
library(lisi)
library(aricode)
library(cluster)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(patchwork)
library(fpc)
library(SingleCellExperiment)
library(scuttle)

set.seed(123)

counts <- read.csv("C:/Users/86158/Desktop/论文复现1/DATA/AfterQC_count.csv", row.names = 1)
metadata <- read.csv("C:/Users/86158/Desktop/论文复现1/metadata_with_superclusters.csv", row.names = 1)


common_cells <- intersect(colnames(counts), rownames(metadata))
seurat_obj <- CreateSeuratObject(
  counts = counts[, common_cells],
  meta.data = metadata[common_cells, ]
)

VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), group.by = "model")

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("nFeature_RNA", "percent.mt"))
seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)

#' ComBat校正
run_combat <- function(seurat){
  tryCatch({
    # 强制锁定操作对象为RNA assay
    DefaultAssay(seurat) <- "RNA"
    
    # 创建组合批次变量（增加错误处理）
    if(!all(c("model", "Phase") %in% colnames(seurat@meta.data))){
      stop("Metadata missing 'model' or 'Phase' columns")
    }
    batch <- interaction(seurat$model, seurat$Phase, drop = TRUE)
    
    exp_matrix <- GetAssayData(
      object = seurat,
      assay = "RNA",     
      layer = "data"      
    )
    
    if(!inherits(exp_matrix, "matrix")){
      exp_matrix <- as.matrix(exp_matrix)
    }
    
    combat_data <- sva::ComBat(
      dat = exp_matrix,
      batch = batch,
      mod = model.matrix(~1, data = seurat@meta.data),
      par.prior = TRUE
    )
    
    corrected <- CreateSeuratObject(
      counts = GetAssayData(seurat, assay = "RNA", layer = "counts"),
      assay = "RNA",             
      meta.data = seurat@meta.data  
    )
    
    corrected <- SetAssayData(
      object = corrected,
      assay = "RNA",        
      layer = "data",     
      new.data = combat_data
    )
    
    corrected <- NormalizeData(corrected)  
    corrected <- ScaleData(corrected)      
    
    corrected <- FindVariableFeatures(corrected, nfeatures = 2000)
    corrected <- RunPCA(
      corrected,
      npcs = 30,
      features = VariableFeatures(corrected), 
      verbose = FALSE
    )
    
    return(corrected)
  }, error = function(e){
    message("ComBat Error: ", e$message)
    return(NULL)
  })
}
#' Harmony校正
run_harmony <- function(seurat){
  seurat %>% 
    RunHarmony(
      group.by.vars = c("model", "Phase"),
      theta = c(1.2, 0.8), 
      lambda = 0.8,          
      max.iter = 50,
      project.dim = FALSE  
    )
}

#' Seurat CCA整合
run_seurat_cca <- function(seurat){
  tryCatch({
    seurat_list <- SplitObject(seurat, split.by = "model")
    seurat_list <- lapply(seurat_list, function(x){
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x)
      return(x)
    })
    
    anchors <- FindIntegrationAnchors(
      object.list = seurat_list,
      dims = 1:30,
      k.anchor = 10
    )
    
    integrated <- IntegrateData(
      anchorset = anchors,
      dims = 1:30,
      new.assay.name = "integrated"
    )
    
    DefaultAssay(integrated) <- "integrated"
    integrated <- ScaleData(
      integrated,
      vars.to.regress = c("Phase")
    )
    integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)
    return(integrated)
  }, error = function(e){
    message("Seurat CCA Error: ", e)
    return(NULL)
  })
}



#' fastMNN校正
run_mnn <- function(seurat){
  tryCatch({
    library(SingleCellExperiment)
    library(batchelor)
    library(scuttle)
    
    batches <- SplitObject(seurat, split.by = "model")
    batch_names <- names(batches)
    
    sce_list <- lapply(batch_names, function(x){
      sce <- SingleCellExperiment(
        list(counts = GetAssayData(batches[[x]], assay = "RNA", layer = "counts")),
        colData = batches[[x]]@meta.data
      )
      sce$batch <- x
      logcounts(sce) <- log1p(calculateCPM(sce))
      return(sce)
    })
    names(sce_list) <- batch_names
    
    mnn_out <- fastMNN(
      sce_list,
      d = 30,
      k = 15,
      BSPARAM = BiocSingular::RandomParam(),
      auto.merge = TRUE
    )
    
    corrected_embedding <- reducedDim(mnn_out, "corrected")
    rownames(corrected_embedding) <- colnames(mnn_out)
    
    corrected <- CreateSeuratObject(
      counts = GetAssayData(seurat, assay = "RNA", layer = "counts"),
      assay = "RNA",
      meta.data = seurat@meta.data
    )

    corrected[["mnn"]] <- CreateDimReducObject(
      embeddings = corrected_embedding,
      key = "mnn_",
      assay = "RNA"
    )
    

    corrected <- NormalizeData(corrected)
    corrected <- ScaleData(corrected)  # 生成scale.data
    
    corrected <- FindVariableFeatures(corrected, selection.method = "vst", nfeatures = 2000)
    
    corrected <- RunUMAP(
      corrected,
      reduction = "mnn",
      dims = 1:30,
      reduction.name = "umap"
    )
    
    return(corrected)
    
  }, error = function(e){
    message("MNN Error Details:\n", conditionMessage(e))
    return(NULL)
  })
}



#统一评估函数

evaluate_correction <- function(corrected, method_name, true_labels){
  kbet_res <- kBET(
    corrected@reductions$pca@cell.embeddings[,1:30],
    batch = corrected$model,
    k0 = 15,
    plot = FALSE
  )
  
  lisi_res <- compute_lisi(
    corrected@reductions$pca@cell.embeddings[,1:30],
    corrected@meta.data[, c("model", "Phase")],
    c("model", "Phase")
  )
  
  corrected <- FindNeighbors(corrected, reduction = "pca", dims = 1:30)
  corrected <- FindClusters(corrected, resolution = 0.5)
  clusters <- corrected$seurat_clusters
  
  dist_mat <- dist(corrected@reductions$pca@cell.embeddings[,1:30])
  silhouette <- mean(silhouette(as.numeric(clusters), dist_mat)[,3])
  
  calinski <- calinhara(
    corrected@reductions$pca@cell.embeddings[,1:30],
    as.numeric(clusters))
    
    data.frame(
      Method = method_name,
      kBET = mean(kbet_res$stats$kBET.observed),
      LISI_model = mean(lisi_res$model),
      LISI_Phase = mean(lisi_res$Phase),
      Silhouette = silhouette,
      Calinski = calinski,
      N_Clusters = length(unique(clusters))
      )
}

methods <- list(
  "ComBat" = run_combat,
  "Harmony" = run_harmony,
  "Seurat_CCA" = run_seurat_cca,
  "fastMNN" = run_mnn
)
results <- purrr::map_dfr(names(methods), function(method_name){
  message("\nRunning ", method_name, "...")
  start_time <- Sys.time()
  
  corrected <- methods[[method_name]](seurat_obj)
  if(is.null(corrected)) return(NULL)
  
  if(!"pca" %in% names(corrected@reductions)){
    corrected <- RunPCA(corrected, npcs = 30)
  }
  
  eval_res <- evaluate_correction(corrected, method_name, seurat_obj$cell_type)
  eval_res$Runtime <- difftime(Sys.time(), start_time, units = "mins") %>% as.numeric()
  p <- DimPlot(corrected, group.by = "model") + 
    ggtitle(paste(method_name, "Corrected"))
  print(p)
  
