import os
import sys
import pandas as pd
import scanpy as sc
# import stereo as st

sys.stdout = None


class ManualAnnotation:
    """
    Class for manual annotation of cell types
    st_data: StereopyData
    marker_gene_key: str
    marker_df: pd.DataFrame
    marker_num_of_cluster: int
    cell_dict: dict
    ref: pd.DataFrame
    ref_dict: dict
    marker_gene_in_data: dict
    select_keys: list
    """

    def __init__(
        self,
        st_data,
        ref_path,
        select_keys=None,
        marker_gene_key="marker_genes",
        marker_num_of_cluster=10
    ):
        self.st_data = st_data
        self.an_data = st.io.stereo_to_anndata(
            st_data,
            flavor="scanpy"
        )
        if os.path.isfile(ref_path):
            self.ref = pd.read_excel(ref_path)
        else:
            print(f"Error: {ref_path} not found")
            sys.exit(1)
        self.marker_gene_key = marker_gene_key
        if isinstance(select_keys, tuple) and len(select_keys) == 1:
            self.select_keys = select_keys[0]
        else:
            self.select_keys = select_keys
        self.marker_num_of_cluster = marker_num_of_cluster
        self.ref_dict = {}
        self.cell_dict = {}
        self.marker_gene_in_data = {}
        self.selected_items = {}
        self.marker_df = st_data.tl.result[marker_gene_key]['pct']
        self.figure = None

    def run(self):
        self.load_ref_dict()
        self.load_cell_dict()
        self.init_marker_gene_in_data()
        self.select_genes()
        self.draw()

    def load_cell_dict(self):
        for col in range(1, self.marker_df.shape[1]):
            temp_df = self.marker_df.iloc[:, [0, col]]
            # sort by value of col
            temp_df = temp_df.sort_values(
                by=temp_df.columns[1],
                ascending=False
            )
            # get top genes
            temp_df = temp_df.head(self.marker_num_of_cluster)
            top_genes = temp_df.iloc[:, 0].tolist()
            self.cell_dict[col] = top_genes

    def load_ref_dict(self):
        for i in range(self.ref.shape[0]):
            # get cell type
            cell_type = self.ref.iloc[i, 0]
            # if cell type not in cell_dict, add it
            if cell_type not in self.ref_dict:
                self.ref_dict[cell_type] = []
            # get gene
            gene = self.ref.iloc[i, 1]
            # add gene to ref dict
            self.ref_dict[cell_type].append(gene)

    def init_marker_gene_in_data(self):
        for ct, markers in self.ref_dict.items():
            markers_found = list()
            for marker in markers:
                if marker in self.an_data.var.index:
                    markers_found.append(marker)
            self.marker_gene_in_data[ct] = markers_found

    def select_genes(self):
        if self.select_keys is not None:
            tmp_select_keys = [
                list(self.marker_gene_in_data)[i] for i in self.select_keys
            ]
        else:
            tmp_select_keys = list(self.marker_gene_in_data.keys())
        self.selected_items = {
            key: self.marker_gene_in_data[key] for key in tmp_select_keys
        }

    def draw(self):
        sc.tl.dendrogram(
            self.an_data,
            groupby="leiden"
        )
        return sc.pl.DotPlot(
            self.an_data,
            var_names=self.selected_items,
            groupby="leiden",
            standard_scale="var",
        ).savefig('./out/annotation/marker_plot.pdf')


sys.stdout = sys.__stdout__

# cell2location annotation
class Cell2locAnnotation:
    """
    Class for cell2location annotation
    sapatial_data: annData
    ref_path: str
    """
    def __init__(
        self,
        spatial_data,
        ref_path
    ):
        self.spatial_data = spatial_data
        self.ref_data = sc.read_h5ad(ref_path)

    def train(self):
        from cell2location.models import RegressionModel
        from scvi.model import SCVI

        # 初始化和训练单细胞模型
        sc_model = SCVI(self.ref_data)
        sc_model.train()

        # 创建参考模型
        self.reference_model = RegressionModel.from_rna_model(sc_model)
        self.reference_model.train(self.spatial_data)

    def predict(self):
        import cell2location
        # 利用参考模型进行空间组注释
        self.results = cell2location.run_cell2location(
            spatial_data=self.spatial_data,
            reference_model=self.reference_model,
        )
    
    def get_results(self):
        return self.results
    
    def save_results(self, path):
        self.results.to_csv(path)
    
    def plot_heatmap(self, path):
        import seaborn as sns
        import matplotlib.pyplot as plt
        sns.heatmap(
            self.results.uns["spatial_coefs"],
            cmap="coolwarm",
            center=0,
            cbar_kws={"label": "Coefficient"},
        )
        plt.savefig(path)

    def plot_spacial(self, path):
        sc.pl.spatial(
            self.spatial_data,
            color="cell2location",
            spot_size=10,
            alpha_img=0.5,
            save=path
        )
