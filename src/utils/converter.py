import src.utils.config as config


def save_ster2h5ad(data, file, flavor):
    """
    save stereo data to h5ad file and return anndata object
    """
    import stereo as st
    return st.io.stereo_to_anndata(
        data,
        flavor=flavor,
        output=file
    )


def conv_adata2ster(data, keys=None):
    """
    convert AnnData to stereo data
    """
    import stereo as st
    return st.io.anndata_to_stereo(
        data,
        spatial_key=keys,
        resolution=config.chip_resolution
    )


def conv_rds2h5ad(rds_filepath, adata_filepath):
    """
    convert rds file to h5ad file
    """
    import pyreadr
    import anndata
    import sys
    import pandas as pd

    print('loading data ...')
    data = pyreadr.read_r('data/reference_gene.rds')
    data = data[None].transpose()
    print('data loaded')
    # 如果是一个df
    if isinstance(data, pd.DataFrame):
        print('data is a DataFrame')
        # try to convert it to anndata
        try:
            data = anndata.AnnData(data)
        except Exception:
            print('data is a DataFrame, but can not be converted to anndata')
            sys.exit(1)
        # save it to h5ad
        print('save it to h5ad')
        data.write('data/reference_gene.h5ad')
        print('data is saved')
    else:
        print('data is not a df')
