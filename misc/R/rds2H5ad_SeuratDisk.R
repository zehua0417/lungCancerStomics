library(Seurat)
library(SeuratDisk)

data <- readRDS(file = "data/reference_gene.rds")
SaveH5Seurat(data, filename = "../../data/reference_gene.h5Seurat")
Convert("../../data/reference_gene.h5Seurat", dest = "h5ad")
