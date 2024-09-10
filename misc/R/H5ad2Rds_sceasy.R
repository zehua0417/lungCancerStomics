#! /hsfscqjf1/ST_CQ/Reference/hxh/miniconda3/envs/R4.2.1/bin/R

library(anndata)
library(sceasy)
setwd("/hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/test/lungcancer")

#ad <- anndata::read_h5ad('temp/temp.h5ad')
sceasy::convertFormat("temp/temp.h6ad", from="anndata", to="seurat", outFile='temp/temp.rds')
