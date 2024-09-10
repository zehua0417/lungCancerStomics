#! /hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/tools/anaconda/envs/r_env/bin/Rscript

cat("> 0. init\n")
setwd("/hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/test/lungcancer")
# setwd("/data/work/lungcancerstomics")

options(warn = -1)
library(data.table)
library(Matrix)
library(ggplot2)
library(Seurat)
library(SingleCellExperiment)
library(SPOTlight)
library(scater)
library(scran)
library(scatterpie)
options(warn = 0)


cat("> 3. loading data ...\n")
tissue_data <- readRDS("temp/temp_bin20.rds")
seurat_ref <- readRDS("data/12T_harmony_celltype.rds")
x <- SingleCellExperiment(seurat_ref)
y <- SingleCellExperiment(tissue_data)
y_matrix <- GetAssayData(tissue_data, assay = "RNA", layer = "data")
y <- SingleCellExperiment(list(counts = y_matrix))
rm(y_matrix)



cat("> 2. preprocessing ...\n")
sce <- as.SingleCellExperiment(seurat_ref)
sce <- logNormCounts(sce)

cat("> 2.1. high variable genes ...\n")
# 去掉核糖体和线粒体基因
gene <- !grepl(
  pattern = "^RP[L|S]|MT",
  x = rownames(sce)
)
dec <- modelGeneVar(sce , subset.row = gene)

cat("> 2.2. marker genes ...\n")
# 计算高变基因
hvg <- getTopHVGs(dec, n = 1500)
# 加上细胞注释信息
colLabels(sce) <- colData(sce)$celltype
# Compute marker genes
mgs <- scoreMarkers(sce, subset.row = gene)
# 保留最相关的marker基因
mgs_fil <- lapply(names(mgs),
  function(i) {
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.8
    x <- x[x$mean.AUC > 0.8, ]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
  }
)
mgs_df <- do.call(rbind, mgs_fil)


cat("> 2.3. Downsample ...\n")
idx <- split(seq(ncol(sce)), sce$celltype)
n_cells <- 2000
cs_keep <- lapply(idx,
  function(i) {
    n <- length(i)
    if (n < n_cells)
      n_cells <- n
    sample(i, n_cells)
  }
)
sce <- sce[, unlist(cs_keep)]

cat("> 2.4. data celltype ...\n")
groups <- sce$celltype

# rm(seurat_ref, sce, dec)
cat("> 3. save/load temp data\n")
# save(sce, hvg, mgs_df, groups, file = "temp/spotlight_data.rda")
load("temp/spotlight_data.rda")
tissue_data <- readRDS("temp/temp_bin20.rds")
y <- SingleCellExperiment(tissue_data)
y_matrix <- GetAssayData(tissue_data, assay = "RNA", layer = "data")
y <- SingleCellExperiment(list(counts = y_matrix))
rm(y_matrix)
res <- readRDS("temp/spotlight_res.rds")


cat("> 4. running spotlight ...\n")
res <- SPOTlight(
  x = sce,
  y = y,
  groups = groups,
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene"
)
cat("> 5. saving ... \n")
rm(sce, y, groups, hvg)
# saveRDS(res, file = "temp/spotlight_res.rds")
df_merged <- merge(res$mat, tissue_data@reductions$spatial@cell.embeddings, all.x = TRUE)
write.csv(df_merged, "out/annotation/tissue_anno.csv", row.names = TRUE)



# cat("> 6. plot \n")
# mat <- res$mat
# ct <- colnames(mat)
# mat[mat < 0.1] <- 0
# palettemartin <- c(
#   "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db",
#   "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff",
#   "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d"
# )
# pal <- colorRampPalette(palettemartin)(length(ct))
# names(pal) <- ct
# x_coords <- as.matrix(tissue_data@reductions$spatial@cell.embeddings)
# plot <- plotSpatialScatterpie(
#   x = x_coords,
#   y = mat,
#   cell_types = colnames(mat),
#   img = FALSE,
#   scatterpie_alpha = 1,
#   pie_scale = 0.4
# )
# # + scale_fill_manual(
# #   values = pal,
# #   breas = names(pal)
# # )
# ggsave(filename = "temp/spatial_scatterpie_plot.png", plot = plot, width = 10, height = 8, dpi = 300)


