#!/bin/bash

# R
R_ENV_NAME="r_env"
R_PACKAGES="r-base=4.3.3 r-essentials"
conda create -n $R_ENV_NAME $R_PACKAGES -y
source activate $R_ENV_NAME
Rscript -e "install.packages(c('data.table', 'Matrix', 'ggplot2'))"
Rscript -e "install.packages('BiocManager')"
Rscript -e "BiocManager::install(c('reticulate', 'Seurat', 'scater', 'scran', 'SingleCellExperiment', 'SPOTlight'))"
conda deactivate

# python
PYTHON_ENV_NAME="stereo"
conda create -n $PYTHON_ENV_NAME python=3.8 -y
source activate $PYTHON_ENV_NAME
pip install matplotlib 
pip install seaborn
pip install anndata
pip install stereopy
pip install pyreadr
conda deactivate

# 公司网不好, 见笑了>_
# 函数：尝试创建环境并安装包
# create_env() {
# 	echo "尝试创建 Conda 环境并安装包..."
# 	conda create -n $R_ENV_NAME $PACKAGES -y
# }
# 无限循环，直到成功创建环境并安装包
# while true; do
# 	create_env
# 	if [ $? -eq 0 ]; then
# 		echo "成功创建 Conda 环境并安装包。"
# 		break
# 	else
# 		echo "创建 Conda 环境失败，重试..."
# 	fi
# done
