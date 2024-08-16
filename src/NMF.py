import src.utils.config as config
import src.utils.loader as loader
import src.utils.imgcatcher as imgc
import pandas as pd
from sklearn.decomposition import NMF
from matplotlib import pyplot as plt
import seaborn as sns

###  1. load data
print("loading data ...")
tissue_data = loader.load_data("tissue_temp")
cell_data = tissue_data.tl.result['annotation']

# init NMF model
component = config.nmf_component
model = NMF(
		n_components = component,
		init = 'nndsvd',  # 初始化方法
		random_state = 123,  # 随机种子
		max_iter = 1000,  # 最大迭代次数
		verbose = 0  # 是否打印输出
		)

# 拟合模型并转换数据
W = model.fit_transform(cell_data)
H = model.components_  # 模型的成分矩阵

# 转换为DataFrame
W_df = pd.DataFrame(
		W,
		index = cell_data.index,
		columns = [f"component_{i}" for i in range(component)]
		)
H_df = pd.DataFrame(
		H,
		index = [f"component_{i}" for i in range(component)],
		columns = cell_data.columns
		)

# 保存结果
W_df.to_csv('out/nmf/W.csv')
H_df.to_csv('out/nmf/H.csv')

# 可视化
imgcatcher_heatmap = imgc.ImgCatcher('out/nmf/heatmap.pdf')
plt.figure(figsize=(10, 10))
sns.heatmap(H_df, cmap='viridis')
plt.title('NMF Heatmap')
plt.xlabel('Genes')
plt.ylabel('Components')
imgcatcher_heatmap.save_and_close()

