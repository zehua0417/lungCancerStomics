#! /hsfscqjf1/ST_CQ/Reference/hxh/miniconda3/envs/R4.2.1/bin/Rscript

.libPaths("/hsfscqjf1/ST_CQ/Reference/hxh/miniconda3/envs/R4.2.1/lib/R/library")
library(Seurat)
library(SeuratDisk)
setwd("/hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/test/lungcancer/temp/")

Convert("temp_bin20.h5ad", dest = "h5Seurat", overwrite = TRUE)

seurat_object <- LoadH5Seurat("temp_bin20.h5seurat")

saveRDS(seurat_object, file = "temp_bin20.rds")
