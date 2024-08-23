# 基于官方的带有Conda的基础镜像
FROM continuumio/miniconda3

# 配置conda使用清华镜像源
RUN conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/ && \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/ && \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/ && \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r/ && \
    conda config --set show_channel_urls yes

# 设置R环境变量
ENV R_ENV_NAME="r_env"
ENV R_PACKAGES="r-base=4.3.3 r-essentials"

# 创建R环境并安装基本包，同时设置R的CRAN镜像
RUN conda create -n $R_ENV_NAME $R_PACKAGES -y && \
    echo 'options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))' > /opt/conda/envs/$R_ENV_NAME/lib/R/etc/Rprofile.site && \
    conda run -n $R_ENV_NAME Rscript -e "install.packages(c('data.table', 'ggplot2'))" && \
    conda run -n $R_ENV_NAME Rscript -e "install.packages('BiocManager')"
RUN conda run -n $R_ENV_NAME Rscript -e "BiocManager::install(c('reticulate', 'Seurat', 'scater', 'scran', 'SingleCellExperiment', 'SPOTlight'))"

# 创建Python环境并安装指定的Python包
ENV PYTHON_ENV_NAME="stereo"
ENV PIP_INDEX_URL=https://pypi.tuna.tsinghua.edu.cn/simple

RUN conda create -n $PYTHON_ENV_NAME python=3.8 -y && \
    conda run -n $PYTHON_ENV_NAME pip install stereopy -i https://pypi.tuna.tsinghua.edu.cn/simple && \
    # conda install -n $PYTHON_ENV_NAME -c conda-forge -c stereopy stereopy
    conda install -n $PYTHON_ENV_NAME -c conda-forge matplotlib seaborn anndata pyreadr

# 配置基本工具
RUN apt-get update && apt-get install -y make && apt-get install -y vim

# 创建一个新的用户，并设置其工作目录
RUN useradd -m -s /bin/bash stereonote

# 将用户切换为新创建的用户
USER stereonote

# 设置工作目录
WORKDIR /home/stereonote/workspace

# 复制脚本到镜像中
COPY . .

# 初始化
RUN make init

# 使用bash作为默认shell
SHELL ["/bin/bash", "-c"]

CMD ["bash"]
