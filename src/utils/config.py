import os
import sys
import numpy as np

import warnings

warnings.filterwarnings('ignore')

#* stream control parameters *#
#status = "dev"	# 开发阶段, bin_size = 100
status = "test"  # 测试阶段, bin_size = 50
debug = False  # 调试, 使用小数据集

#* data path *#
if debug:
    tissue_path = "data/A03982E1_microbiome.genus.label.gem"
else:
    tissue_path = "data/A03982E1.tissue.gem"
genus_path = "data/A03982E1_microbiome.genus.label.gem"
species_path = "data/A03982E1_microbiome.species.label.gem"

if not os.path.isfile(tissue_path):
    print(f"文件{tissue_path}不存在")
    sys.exit(1)
if not os.path.isfile(genus_path):
    print(f"文件{genus_path}不存在")
    sys.exit(1)
if not os.path.isfile(species_path):
    print(f"文件{species_path}不存在")
    sys.exit(1)

############## parameters ##############
#* bin_size *#
if status == "dev":
    bin_size = 100
else:
    bin_size = 50

#* filter cells and genes *#
# minimum/maximum number of counts required for a cell to pass fitlering.
# 细胞通过过滤所需的最小/大表达计数
min_cell_counts = 300
max_cell_counts = 3000
# minimum/maximum number of counts expressed required for a gene to pass filtering.
# 基因通过过滤所需的最小/大表达计数
min_gene_counts = 2500
max_gene_counts = 12500
# minimum/maximum number of genes expressed required for a cell to pass filtering.
# 细胞通过过滤所需的最小/大表达基因数
min_genes = 200
max_genes = 4300
# minimum/maximum number of cells expressed required for a gene to pass filering:w
# .
# 基因通过过滤所需的最小/大表达细胞数
min_cells = None
max_cells = None
# maximum number of pct_counts_mt required for a cell to pass filtering.
# 细胞能通过过滤的最大pct_counts_mt
pct_counts_mt = 1.5
# mean counts greater than this value for a gene to pass filtering.
# 基因通过过滤所需的平均表达计数
mean_umi_gt = 0.1

#* high variable genes *#
method = "seurat"  # 'seurat' or 'cell_ranger' or 'seurat_v3'
# 要保留的顶级高度变异基因的数量
n_top_genes = None
# 基因表达值的最小/大正态化离散度
min_disp = 0.5
max_disp = np.inf
# 基因表达值的最小/大平均值
min_mean = 0.0125
max_mean = 3
# Loess模型拟合中使用的数据比例
span = 0.3
n_bins = bin_size