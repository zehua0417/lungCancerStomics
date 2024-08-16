#!/hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/tools/anaconda/envs/stereo/bin/python

import anndata
import scanpy as sc

adata = anndata.read_h5ad('/hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/test/lungcancer/data/reference_gene.h5ad')

print(data.X[:5, :])
