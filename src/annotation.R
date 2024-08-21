#! /hsfscqjf1/ST_CQ/Reference/hxh/miniconda3/envs/R4.2.1/bin/Rscript

cat("init\n")
# setwd("/hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/test/lungcancer")
setwd("/data/work/lungcancer")

library(data.table)
library(Matrix)
library(ggplot2)
library(reticulate)
library(Seurat)
library(SPOTlight)
library(scater)
library(scran)
library(SingleCellExperiment)

config <- import_from_path("config", "src/utils")

fi (config$use != "stereonote"){
    use_condaenv(condaenv = "/hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/tools/anaconda/envs/stereopy-rapids/bin/python", require = TRUE)
}
ad <- import("anndata")


cat("1. loading data ...\n")
adata <- ad$read_h5ad(config$tissue_seurat_temp_path)
seurat_ref <- readRDS(config$reference_rds_path)

expr_matrix <- t(adata$X)
colnames(expr_matrix) <- as.character(adata$obs_names$tolist())
rownames(expr_matrix) <- as.character(adata$var_names$tolist())
# colnames(adata$X) <- adata$var$features
# 
# var_names <- as.character(adata$var_names$tolist())
# obs_names <- as.character(adata$obs_names$tolist())
# 
# expr_matrix <- as(adata$X, "CsparseMatrix")
# rownames(expr_matrix) <- obs_names
# #colnames(expr_matrix) <- var_names
# seurat_obj <- CreateSeuratObject(counts = expr_matrix)





cat("2. preprocessing ...\n")
sce <- as.SingleCellExperiment(seurat_ref)
sce <- logNormCounts(sce)

# # 获取sce和seurat_obj中的共有基因集
# common_gene <- intersect(rownames(sce), rownames(seurat_obj))
# 
# # 过滤sce和seurat_obj以仅保留共有基因
# sce <- sce[common_gene, ]
# expr_matrix <- expr_matrix[common_gene, ]
# seurat_obj <- CreateSeuratObject(counts = expr_matrix)
# sce <- as.SingleCellExperiment(seurat_ref)

# 去掉核糖体和线粒体基因
gene <- !grepl(
    pattern = "^RP[L|S]|MT",
    x = rownames(sce)
)
dec <- modelGeneVar(sce , subset.row = gene)

# 计算高变基因
hvg <- getTopHVGs(dec, n = config$n_top_genes)

# 加上细胞注释信息
colLabels(sce) <- colData(sce)$celltype

# Compute marker genes
mgs <- scoreMarkers(sce, subset.row = gene)

# 保留最相关的marker基因
mgs_fil <- lapply(names(mgs), function(i) {
  x <- mgs[[i]]
  # Filter and keep relevant marker genes, those with AUC > 0.8
  x <- x[x$mean.AUC > 0.8, ]
  # Sort the genes from highest to lowest weight
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  # Add gene and cluster id to the dataframe
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)




cat("3. spotlight runing ...\n")
x_matrix <- GetAssayData(seurat_ref, assay = "RNA", layer = "data")
# x_matrix <- as(x_matrix, "RsparseMatrix")
x <- SingleCellExperiment(list(counts = x_matrix))
y <- SingleCellExperiment(list(counts = expr_matrix))
groups <- sce$celltype
print(class(x_matrix))
print(class(expr_matrix))
rm(seurat_ref,adata,expr_matrix,x_matrix,sce,dec)

res <- SPOTlight(
  x = x,
  y = y,
  groups = groups,
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene"
)

cat("saving ... \n")
saveRDS(res, file = config$spotlight_res_path)

