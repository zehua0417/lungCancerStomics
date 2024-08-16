import os
import sys
import pandas as pd
import scanpy as sc
import stereo as st

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


