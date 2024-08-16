#!F:/Program\ Data/condaEnvs/stereo/python
#!/hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/tools/anaconda/envs/stereopy-rapids/bin/python

print("init ...")
import sys
import os

os.chdir("/hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/test/lungcancer")
# os.chdir("f:/onedrive/study/biology/stomics/lungcancer")
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, project_root)

import src.utils.annotation as anno
import src.utils.filter as ftr
import src.utils.imgcatcher as imgc
import src.utils.loader as loader
import src.utils.converter as cvter
import src.utils.config as config


print("loading data ...")
tissue_data = loader.load_data("tissue", "gem")

################ I preprocess ################
# * 1. calc qc and visualize *#
print("1.1. calc qc and visualize")

tissue_data.tl.cal_qc()

imgcatcher_scatter = imgc.ImgCatcher('out/preprocess/scatter_qc1.pdf')
tissue_data.plt.spatial_scatter()
imgcatcher_scatter.save_and_close()
del imgcatcher_scatter

imgcatcher_violin = imgc.ImgCatcher('out/preprocess/violin_qc1.pdf')
tissue_data.plt.violin(
    keys=['n_genes_by_counts', 'total_counts'],
    show_stripplot=True
)
imgcatcher_violin.save_and_close()
del imgcatcher_violin

imgcatcher_violin_mt = imgc.ImgCatcher('out/preprocess/violin_qc_mt1.pdf')
tissue_data.plt.violin(
    keys=['pct_counts_mt'],
    show_stripplot=True
)
imgcatcher_violin_mt.save_and_close()
del imgcatcher_violin_mt

# * 2. filter cells and genes *#
print("1.2. filtering ...")
# observe cell distribution
imgcatcher_scatter = imgc.ImgCatcher('out/preprocess/scatter_gene_count1.pdf')
tissue_data.plt.genes_count()
imgcatcher_scatter.save_and_close()
del imgcatcher_scatter

# filter by cells
tissue_data.tl.filter_cells(
    min_counts=config.min_cell_counts,
    max_counts=config.max_cell_counts,
    min_genes=config.min_genes,
    max_genes=config.max_genes,
    pct_counts_mt=config.pct_counts_mt,
    inplace=True
)

# filter zero mt cells
ftr.filter_zero_mt_cells(tissue_data)

# filter by genes
# tissue_data.tl.filter_genes(
#     min_counts=config.min_gene_counts,
#     max_counts=config.max_gene_counts,
#     min_cells=config.min_cells,
#     max_cells=config.max_cells,
#     mean_umi_gt=config.mean_umi_gt,
#     filter_mt_genes=True,
#     inplace=True
# )

# observe cell distribution after filter
imgcatcher_scatter = imgc.ImgCatcher('out/preprocess/scatter_gene_count2.pdf')
tissue_data.plt.genes_count()
imgcatcher_scatter.save_and_close()
del imgcatcher_scatter

imgcatcher_violin = imgc.ImgCatcher('out/preprocess/violin_qc2.pdf')
tissue_data.plt.violin(
    keys=['n_genes_by_counts', 'total_counts'],
    show_stripplot=True
)
imgcatcher_violin.save_and_close()
del imgcatcher_violin

imgcatcher_violin_mt = imgc.ImgCatcher('out/preprocess/violin_qc_mt2.pdf')
tissue_data.plt.violin(
    keys=['pct_counts_mt'],
    show_stripplot=True
)
imgcatcher_violin_mt.save_and_close()
del imgcatcher_violin_mt

imgcatcher_scatter = imgc.ImgCatcher('out/preprocess/scatter_qc2.pdf')
tissue_data.plt.spatial_scatter()
imgcatcher_scatter.save_and_close()
del imgcatcher_scatter

sys.exit(0)

# save data to self.raw
print("archive data ...")
tissue_data.tl.raw_checkpoint()
print("==tissue_data.tl.raw==")
print(tissue_data.tl.raw)

# * Normalization *#
print("1.3. normalization ...")
tissue_data.tl.normalize_total()
tissue_data.tl.log1p()

################ II High Variability Genes ################
print("2.1. high variability genes ...")
tissue_data.tl.highly_variable_genes(
    n_top_genes=config.n_top_genes,
    min_disp=config.min_disp,
    max_disp=config.max_disp,
    min_mean=config.min_mean,
    max_mean=config.max_mean,
    res_key="highly_variable_genes"
)

imgcatcher_scatter = imgc.ImgCatcher('out/hvg/scatter_hvg1.pdf')
tissue_data.plt.highly_variable_genes(res_key="highly_variable_genes")
imgcatcher_scatter.save_and_close()
del imgcatcher_scatter

tissue_data.tl.scale(zero_center=False)

################ III Embedding ################
print("3.1. PCA ...")
# from sklearn.decomposition import PCA
# hvgs=tissue_data.tl.result["highly_variable_genes"]['highly_variable']
# exp_matrix=tissue_data.tl.data.exp_matrix[:, hvgs]
# pca=PCA(n_components=config.n_pcs)
# exp_dense=exp_matrix.toarray()
# fitted=pca.fit(exp_dense)
# fitted.n_components_
tissue_data.tl.pca(
    use_highly_genes=True,
    n_pcs=config.n_pcs,
    res_key="pca"
)
tissue_data.tl.neighbors(
    pca_res_key="pca",
    res_key="neighbors",
    method=config.method
)

print("3.2. UMAP ...")
tissue_data.tl.umap(
    pca_res_key="pca",
    neighbors_res_key="neighbors",
    res_key="umap",
    # method=config.method
)

# imgcatcher_umap=imgc.ImgCatcher('out/umap/umap_gene.pdf')
# tissue_data.plt.umap(
#     gene_names=['AASS'],
#     res_key="umap"
# )
# imgcatcher_umap.save_and_close()
# del imgcatcher_umap

print("3.3. Leiden ...")
tissue_data.tl.leiden(
    resolution=config.resolution,
    neighbors_res_key="neighbors",
    res_key="leiden",
    method=config.method
)
imgcatcher_leiden = imgc.ImgCatcher('out/leiden/leiden.pdf')
tissue_data.plt.cluster_scatter(
    res_key="leiden"
)
imgcatcher_leiden.save_and_close()
del imgcatcher_leiden
imgcatcher_umap = imgc.ImgCatcher('out/umap/umap.pdf')
tissue_data.plt.umap(
    res_key="umap",
    cluster_key="leiden"
)
imgcatcher_umap.save_and_close()
del imgcatcher_umap

################ IV Find Marker genes ################
print("4.1. marker genes ...")
tissue_data.tl.find_marker_genes(
    cluster_res_key="leiden",
    method="wilcoxon_test",
    use_highly_genes=True
)

# print top 10 marker genes score
imgcatcher_marker = imgc.ImgCatcher('out/marker/genes_text.pdf')
tissue_data.plt.marker_genes_text(
    res_key='marker_genes',
    markers_num=10,
    sort_key='scores'
)
imgcatcher_marker.save_and_close()
del imgcatcher_marker
# scat top 10 marker genes of each cluster
imgcatcher_marker = imgc.ImgCatcher('out/marker/genes_scatter.pdf')
tissue_data.plt.marker_genes_scatter(
    res_key='marker_genes',
    markers_num=10
)
imgcatcher_marker.save_and_close()
del imgcatcher_marker

cvter.save_ster2h5ad(tissue_data, 'temp/temp.h5ad', flavor='anndata')
cvter.save_ster2h5ad(tissue_data, 'temp/temp.h5ad', flavor='seurat')
# del tissue_data
# reload data
# tissue_data=loader.load_data("tissue_temp", "h5ad")

################ V Annotation ################
# print("5.1. manual annotation ...")
# ma_obj=ma.ManualAnnotation(
#     st_data=tissue_data,
#     ref_path=config.ref_path,
#     select_keys=config.select_keys,
#     marker_num_of_cluster=config.marker_num_of_cluster
# )
# ma_obj.run()
# del ma_obj

# print("5.2. auto annotation ...")
# init reference
# ref=loader.load_data("ref_h5", "h5ad")
# ref_dict=loader.load_data("ref_dict", "txt.gz")
# ref_dict=ref_dict[['Index', 'Cell_type.refined']]
# ref_dict=ref_dict.rename(columns={'Index': 'index', 'Cell_type.refined': 'cell_type'})
# cell_index=ref.adata.obs.index
# ref_dict.set_index('index', inplace=True)
# cell_type_list=[ref_dict.loc[idx, 'cell_type'] if idx in ref_dict.index else None for idx in cell_index]
# ref.adata.obs['cell_type']=cell_type_list
#
# tissue_data.tl.single_r(
# 	ref,
# 	ref_use_col='cell_type',
# 	method=config.method,
# 	n_jobs=config.singleR_n_jobs,
# 	fine_tune_times=config.fine_tune_times,
# 	res_key='annotation'
# )
#
# imgcatcher_annotation =imgc.ImgCatcher('out/annotation/annotation.pdf')
# tissue_data.plt.cluster_scatter(
# 	res_key='annotation'
# )
# imgcatcher_annotation.save_and_close()
#
#
# 微生物在淋巴结的空间分布
