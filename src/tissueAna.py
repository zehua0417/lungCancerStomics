#!F:/Program\ Data/condaEnvs/stereo/python
#!/hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/tools/anaconda/envs/stereopy-rapids/bin/python

import os
import sys
print("init ...")

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
tissue_data.tl.cal_qc()
ftr.filter_zero_mt_cells(tissue_data)

filter_obj = ftr.Filter(tissue_data)
filter_obj.calc_QC_metric()
filter_obj.plot_QC_metric("before")
filter_obj.filter_low_quality_cells(config.counts_n_metric, config.mt_n_metric, config.highest_mt_pct)
filter_obj.plot_QC_metric("after")
tissue_data = filter_obj.return_st_data()

del filter_obj

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
tissue_data.tl.pca(
    use_highly_genes=True,
    n_pcs=config.n_pcs,
    res_key="pca"
)
img_catcher_elbow = imgc.ImgCatcher('out/umap/elbow.pdf')
tissue_data.plt.elbow(pca_res_key='pca')
img_catcher_elbow.save_and_close()
del img_catcher_elbow

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
cvter.save_ster2h5ad(tissue_data, 'temp/temp_seurat.h5ad', flavor='seurat')
del tissue_data
# reload data
# tissue_data = loader.load_data("tissue_temp", "h5ad")

################ V Annotation ################
# print("5.1. manual annotation ...")
# ma_obj=anno.ManualAnnotation(
#     st_data=tissue_data,
#     ref_path=config.ref_path,
#     select_keys=config.select_keys,
#     marker_num_of_cluster=config.marker_num_of_cluster
# )
# ma_obj.run()
# del ma_obj

print("5.2. auto annotation ...")
# init reference
ref = loader.load_data("ref_h5", "h5ad")
# ref_dict = loader.load_data("ref_dict", "txt.gz")
# ref_dict = ref_dict[['Index', 'Cell_type.refined']]
# ref_dict = ref_dict.rename(columns={'Index': 'index', 'Cell_type.refined': 'cell_type'})
# cell_index = ref.adata.obs.index
# ref_dict.set_index('index', inplace=True)
# cell_type_list = [ref_dict.loc[idx, 'cell_type'] if idx in ref_dict.index else None for idx in cell_index]
# ref.adata.obs['cell_type'] = cell_type_list

# tissue_data.tl.single_r(
#     ref,
#     ref_use_col='celltype',
#     method=config.method,
#     n_jobs=config.singleR_n_jobs,
#     fine_tune_times=config.fine_tune_times,
#     res_key='annotation'
# )
# 
# imgcatcher_annotation = imgc.ImgCatcher('out/annotation/annotation.pdf')
# tissue_data.plt.cluster_scatter(
#     res_key='annotation'
# )
# imgcatcher_annotation.save_and_close()

