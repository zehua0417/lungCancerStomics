#! /hsfscqjf1/ST_CQ/Reference/hxh/miniconda3/envs/R4.2.1/bin/Rscript

setwd("/hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/test/lungcancer")

library(reticulate)
library(Matrix)
library(Seurat)
use_condaenv(condaenv = "/hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/tools/anaconda/envs/stereopy-rapids/bin/python", require = TRUE)
ad <- import("anndata")

adata <- ad$read_h5ad("temp/temp.h5ad")

colnames(adata$X) <- adata$var$features

var_names <- as.character(adata$var_names$tolist())
obs_names <- as.character(adata$obs_names$tolist())

expr_matrix <- as(adata$X, "CsparseMatrix")
rownames(expr_matrix) <- obs_names
colnames(expr_matrix) <- var_names

seurat_obj <- CreateSeuratObject(counts = expr_matrix)

