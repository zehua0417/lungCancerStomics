import os

os.chdir("f:/onedrive/study/biology/stomics/lungcancer")
# os.chdir("/hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/test/lungcancer")

import src.utils.loader as loader
import src.utils.config as config
import src.utils.Moran as Moran

import pandas as pd
import numpy as np
import scanpy as sc
import stereo as st

# Load the data
micro_data = loader.load_data("genus", "gem")
exp = micro_data._exp_matrix
# convert to np matrix
exp = exp.toarray()
row = exp.shape[0]
col = exp.shape[1]
# convert to list
exp = exp.tolist()

# Calculate the Moran's I
moran = Moran.Moran(exp, row, col)
moran_I = moran.I()
moran_p = moran.p_sim()
print(f"Moran's I: {moran_I}, p-value: {moran_p}")
mpran.save("out/moran/moran.csv")

