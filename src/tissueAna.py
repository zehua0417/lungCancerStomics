#!F:/Program\ Data/condaEnvs/stereo/python
#!/hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/tools/anaconda/envs/stereo/bin/python

print("init ...")
import os
import sys
os.chdir("/hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/test/lungcancer")
#os.chdir("f:/onedrive/study/biology/stomics/lungcancer_new")
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, project_root)

import src.utils.config as config
import src.utils.loader as loader
import src.utils.imgcatcher as imgc
import src.utils.filter as ftr

print("loading data ...")
tissue_data = loader.load_data("tissue")

################ I preprocess ################
#* 1. calc qc and visualize *#
print("1.1. calc qc and visualize")

tissue_data.tl.cal_qc()

imgcatcher_scatter = imgc.ImgCatcher('out/preprocess/scatter_qc1.pdf')
tissue_data.plt.spatial_scatter()
imgcatcher_scatter.save_and_close()

imgcatcher_violin = imgc.ImgCatcher('out/preprocess/violin_qc1.pdf')
tissue_data.plt.violin(
    keys = ['n_genes_by_counts', 'total_counts'],
    show_stripplot=True
)
imgcatcher_violin.save_and_close()

imgcatcher_violin_mt = imgc.ImgCatcher('out/preprocess/violin_qc_mt1.pdf')
tissue_data.plt.violin(
    keys = ['pct_counts_mt'],
    show_stripplot=True
)
imgcatcher_violin_mt.save_and_close()

#* 2. filter cells and genes *#
print("1.2. filtering ...")
# observe cell distribution
imgcatcher_scatter = imgc.ImgCatcher('out/preprocess/scatter_gene_count1.pdf')
tissue_data.plt.genes_count()
imgcatcher_scatter.save_and_close()

# filter by cells
tissue_data.tl.filter_cells(
    min_counts = config.min_cell_counts,
    max_counts = config.max_cell_counts,
    min_genes = config.min_genes,
    max_genes = config.max_genes,
    pct_counts_mt = config.pct_counts_mt,
    inplace = True
)

# filter zero mt cells
ftr.filter_zero_mt_cells(tissue_data)

# filter by genes
tissue_data.tl.filter_genes(
    min_counts = config.min_gene_counts,
    max_counts = config.max_gene_counts,
    min_cells = config.min_cells,
    max_cells = config.max_cells,
    mean_umi_gt = config.mean_umi_gt,
    filter_mt_genes = True,
    inplace = True
)

# observe cell distribution after filter
imgcatcher_scatter = imgc.ImgCatcher('out/preprocess/scatter_gene_count2.pdf')
tissue_data.plt.genes_count()
imgcatcher_scatter.save_and_close()

imgcatcher_violin = imgc.ImgCatcher('out/preprocess/violin_qc2.pdf')
tissue_data.plt.violin(
    keys = ['n_genes_by_counts', 'total_counts'],
    show_stripplot=True
)

imgcatcher_violin.save_and_close()

imgcatcher_violin_mt = imgc.ImgCatcher('out/preprocess/violin_qc_mt2.pdf')
tissue_data.plt.violin(
    keys = ['pct_counts_mt'],
    show_stripplot=True
)
imgcatcher_violin_mt.save_and_close()

imgcatcher_scatter = imgc.ImgCatcher('out/preprocess/scatter_qc2.pdf')
tissue_data.plt.spatial_scatter()
imgcatcher_scatter.save_and_close()

# save data to self.raw
print("archive data ...")
tissue_data.tl.raw_checkpoint()
print("==tissue_data.tl.raw==")
print(tissue_data.tl.raw)


#* Normalization *#
print("2.1. normalization ...")
tissue_data.tl.normalize_total()
tissue_data.tl.log1p()
