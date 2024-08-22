import numpy as np
import scanpy as sc
import seaborn as sns
from stereo.preprocess.qc import cal_cells_indicators
from stereo.core import StereoExpData
from scipy.stats import median_abs_deviation
from src.utils.converter import save_ster2h5ad
from src.utils.converter import conv_adata2ster


def filter_zero_mt_cells(data):
    """
    filter cells with zero mitochondrial gene counts
    please forgive me for not writing the function in the class
    """
    cal_cells_indicators(data)
    cell_subset = np.ones(data.cells.size, dtype=np.bool8)
    cell_subset &= data.cells.pct_counts_mt > 0
    data.sub_by_index(cell_index=cell_subset)


class Filter:
    """
    Filter class
    adata: AnnData object
    """

    def __init__(self, data):
        """
        data: StereoExpData or AnnData object
        """
        self.adata = None
        if isinstance(data, StereoExpData):
            self.adata = save_ster2h5ad(data, flavor = "anndata", file = None)
            self.adata.obsm["spatial"] = data.position
        elif isinstance(data, sc.AnnData):
            self.adata = data
        else:
            raise TypeError("The input data type is not supported")

    def calc_QC_metric(self):
        """
        calculate the QC metric 
        search 3 special genes:
        1. mitochondrial genes
        2. ribosomal genes
        3. hemoglobin genes
        then calculate the Qc metric
        """
        self.adata.var["mt"] = self.adata.var_names.str.startswith("MT-")
        self.adata.var["ribo"] = self.adata.var_names.str.startswith(("RPS", "RPL"))
        self.adata.var["hb"] = self.adata.var_names.str.contains(("^HB[^(P)]"))
        sc.pp.calculate_qc_metrics(
            self.adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
        )

    def plot_QC_metric(self, name):
        """
        plot the QC metric
        """
        p1 = sns.displot(
            self.adata.obs["total_counts"], bins=50, kde=True
        )
        p1.savefig(f"figures/{name}_total_counts.png")
        sc.pl.violin(
            self.adata, "pct_counts_mt", jitter=0.4, multi_panel=True,
            save=f"_{name}_pct_counts_mt.png"
        )
        sc.pl.scatter(
            self.adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",
            save=f"_{name}_qc_scatter.png"
        )

    def is_outlier(self, metric: str, nmads: int):
        """
        detect outliers
        通过计算数据集中某个指标与中位数之间的偏差，来判断该指标的值是否是异常值。
        """
        M = self.adata.obs[metric]
        outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
            np.median(M) + nmads * median_abs_deviation(M) < M
        )
        return outlier

    def filter_low_quality_cells(self, counts_n_mad, mt_n_mad, highest_mt_pct):
        """
        remove low quality cells
        according to:
        - the number of counts per barcode (count depth)
        - the number of genes per barcode
        - the fraction of counts from mitochondrial genes per barcode
        """
        # filter cells with low quality
        self.adata.obs["outlier"] = (
            self.is_outlier("log1p_total_counts", counts_n_mad)
            | self.is_outlier("log1p_n_genes_by_counts", counts_n_mad)
            | self.is_outlier("pct_counts_in_top_20_genes", counts_n_mad)
        )
        # filter cells with high mitochondrial gene expression
        self.adata.obs["mt_outlier"] = self.is_outlier("pct_counts_mt", mt_n_mad) | (
            self.adata.obs["pct_counts_mt"] > highest_mt_pct
        )

        print(self.adata.obs.mt_outlier.value_counts())
        self.adata = self.adata[(~self.adata.obs.outlier) & (~self.adata.obs.mt_outlier)].copy()

        print(f"Number of cells after filtering of low quality cells: {self.adata.n_obs}")

    def return_data(self):
        """
        return the filtered data
        """
        return self.adata

    def return_st_data(self):
        """
        return the filtered data as StereoExpData object
        """
        return conv_adata2ster(self.adata, 'spatial')
