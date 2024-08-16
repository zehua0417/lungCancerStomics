library(Seurat)
library(sceasy)
library(reticulate)
#loompy <- reticulate::import('loompy')
setwd("/hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/test/lungcancer/")

data <- readRDS(file = "data/reference_gene.rds")
print(class(data))

sceasy::convertFormat(
	data,
	from="seurat",
	to="anndata",
	outFile='data/reference_gene.h5ad'
)
