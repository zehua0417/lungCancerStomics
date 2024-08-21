import os
import sys
import numpy as np

import warnings

warnings.filterwarnings('ignore')

user_name = os.getlogin()

# * stream control parameters *#
# status = "dev"	# 开发阶段, bin_size = 100
status = "test"  # 测试阶段, bin_size = 20/50
debug = False  # 调试, 使用小数据集
use_GPU = False

if use_GPU:
    method = 'rapids'
else:
    method = None

# * data path *#
if user_name == "stereonote":  # which means on DCS cloud
    tissue_path = "/data/input/Files/A03982E1/A03982E1.tissue.gem.gz"
    genus_path = "/data/input/Files/A03982E1/A03982E1.microbiome.genus.label.gem"
    species_path = "/data/input/Files/A03982E1/A03982E1.microbiome.species.label.gem"
    tissue_temp_path = "/data/input/Files/lizehua/temp.h5ad"
    tissue_seurat_temp_path = "/data/input/Files/lizehua/temp_seurat.h5ad"
    ref_path = "/data/input/Files/lizehua/ref.xlsx"  # auto annotation ref file
    ref_gene_path = "/data/input/Files/lizehua/reference_gene.h5ad"
    reference_path = "/data/input/Files/lizehua/reference_gene.rds"
    reference_h5_path = "/data/input/Files/lizehua/12T_harmony_celltype.h5ad"
    reference_rds_path = "/data/input/Files/lizehua/12T_harmony_celltype.rds"
    reference_dict_path = "/data/input/Files/lizehua/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    spotlight_res_path = "/data/input/Files/lizehua/spotlight_results.rds"
else:
    if debug:
        tissue_path = "data/A03982E1_microbiome.genus.label.gem"
    else:
        tissue_path = "data/A03982E1.tissue.gem"
    genus_path = "data/A03982E1_microbiome.genus.label.gem"
    species_path = "data/A03982E1_microbiome.species.label.gem"
    tissue_temp_path = "temp/temp.h5ad"
    tissue_seurat_temp_path = "temp/temp_seurat.h5ad"
    ref_path = "data/ref.xlsx"  # auto annotation ref file
    ref_gene_path = "data/reference_gene.h5ad"
    reference_path = "data/reference_gene.rds"
    reference_h5_path = "data/12T_harmony_celltype.h5ad"
    reference_rds_path = "data/12T_harmony_celltype.rds"
    reference_dict_path = "data/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    spotlight_res_path = "temp/spotlight_results.rds"


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
# * bin_size *#
if status == "dev":
    bin_size = 100
else:
    bin_size = 20
chip_resolution = 500

# * filter cells and genes *#
# # minimum/maximum number of counts required for a cell to pass fitlering.
# # 细胞通过过滤所需的最小/大表达计数
# min_cell_counts = 8     # 300
# max_cell_counts = None  # 3000
# # minimum/maximum number of counts expressed required for a gene to pass filtering.
# # 基因通过过滤所需的最小/大表达计数
# min_gene_counts = None  # 2500
# max_gene_counts = None  # 12500
# # minimum/maximum number of genes expressed required for a cell to pass filtering.
# # 细胞通过过滤所需的最小/大表达基因数
# min_genes = None  # 200
# max_genes = None  # 4300
# # minimum/maximum number of cells expressed required for a gene to pass filering:w
# # .
# # 基因通过过滤所需的最小/大表达细胞数
# min_cells = None
# max_cells = None
# # maximum number of pct_counts_mt required for a cell to pass filtering.
# # 细胞能通过过滤的最大pct_counts_mt
# pct_counts_mt = None  # 1.5
# # mean counts greater than this value for a gene to pass filtering.
# # 基因通过过滤所需的平均表达计数
# mean_umi_gt = None  # 0.1
counts_n_metric = 5
mt_n_metric = 3
highest_mt_pct = 8

# * high variable genes *#
method = "seurat"  # 'seurat' or 'cell_ranger' or 'seurat_v3'
# 要保留的顶级高度变异基因的数量
n_top_genes = 2000
# 基因表达值的最小/大正态化离散度
min_disp = 0.5
max_disp = np.inf
# 基因表达值的最小/大平均值
min_mean = 0.0125
max_mean = 3
# Loess模型拟合中使用的数据比例
span = 0.3
n_bins = bin_size

# * Embedding *#
# PCA
# 要计算的主成分数量
n_pcs = 3
# Leiden
# 分辨率参数
# resolution = 1  # 14 clusters
# resolution = 0.5  # 3 clusters
resolution = 0.2   # 6 clusters

# * annotation *#
select_keys = [0, 1, 2, 6, 7, 9, 10, 11, 13, 14, 15, 16, 17, 18, 21, 23],
marker_num_of_cluster = 10

fine_tune_times = 5
singleR_n_jobs = 16
