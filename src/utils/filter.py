import numpy as np
import stereo as st
from stereo.preprocess.qc import cal_cells_indicators

def filter_zero_mt_cells(data):
    """ 
    Filter zero mt cells.
    Refer to the code of stereopy preprocess.filter.py
    """
    cal_cells_indicators(data)

    cell_subset = np.ones(data.cells.size, dtype=np.bool8)
    cell_subset &= data.cells.pct_counts_mt > 0
    data.sub_by_index(cell_index=cell_subset)
    return data
